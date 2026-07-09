# Installing on an air-gapped cluster

The cluster has no internet. Every conda environment, every pip wheel and every
reference database is therefore built or downloaded elsewhere and shipped in. No
package may come from the `defaults` or `anaconda` channels.

## Why the channels matter

`defaults` and `anaconda` are not neutral mirrors of conda-forge. They carry
Anaconda's commercial terms of service, and they ship different builds of the
same version number. An environment that silently resolves one package from
`defaults` is not the environment a reader reproducing this analysis will get,
even from the identical YAML file. bioconda is built and tested against
conda-forge at **strict** channel priority; other configurations solve by
accident.

So: every file in `envs/` lists `conda-forge`, `bioconda`, `nodefaults` in that
order, `config/condarc` pins strict priority and empties `default_channels`, and
`scripts/sh/build_envs.sh` passes `--override-channels` so that a stray
`~/.condarc` on the build machine cannot reintroduce them. That is three
independent guards, because each alone is bypassable.

`python scripts/python/check_env_channels.py envs/` refuses to build if any of
them is violated, and `bash tests/run_tests.sh` runs it. After the build,
`scripts/sh/verify_envs.sh` reads `conda list --explicit` for every environment
and checks the URL each package actually came from — the difference between "we
asked for conda-forge" and "we got conda-forge".

## What makes offline deployment work

Snakemake stores each environment in `<conda-prefix>/<hash>`, where `<hash>` is
the MD5 of the environment **file's content**. Not of the prefix path, not of
the machine. So a directory built anywhere is already named what the cluster
will look for.

What does *not* travel is the absolute path baked inside a conda environment:
shebang lines, `pkg-config` files, R's `Makeconf`. Two ways around it.

**Same path (simplest).** Pick a prefix you can create on both machines, build
there, `rsync` it across. Nothing is rewritten because nothing moved.

**conda-pack (when the paths cannot match).** `build_envs.sh --pack` produces
one relocatable tarball per environment; `deploy_envs.sh --packed` unpacks each
and runs `conda-unpack`, which rewrites the embedded paths.

The build machine must be **Linux x86_64**. `qhost` reports `lx-amd64` for every
node; a conda environment built on macOS or aarch64 will not run there, and
conda-pack will not rescue it.

## Build (machine with internet)

```bash
# bootstrap: mamba + snakemake, conda-forge only
mamba create -p ./bootstrap -c conda-forge --override-channels snakemake=8.20 mamba
conda activate ./bootstrap

git clone <this repo> && cd human-gut-niche-pipeline
export CONDARC=$PWD/config/condarc

# choose a prefix that will exist on the cluster too
bash scripts/sh/build_envs.sh /shared/hgn/conda-envs
#   ... or, if the paths cannot match:
bash scripts/sh/build_envs.sh /shared/hgn/conda-envs --pack
```

That script, in order: rejects any env file that could reach `defaults`; creates
the Snakemake launcher environment; `pip download`s the wheels the cluster cannot
fetch (`snakemake-executor-plugin-sge`); runs `snakemake --sdm conda
--conda-create-envs-only` for **every** target, since rules reachable only from
`transition_all` or `synthesis_all` are otherwise missed and discovered by the
first cluster job; records the env-file → hash mapping; verifies no package came
from a banned channel; and fetches the KEGG module and gene-set tables, which
need the network.

Reference databases are large and separate. Download them on the same machine and
rsync them to the paths named in `config/config.yaml`:

| Data | Config key | Notes |
|---|---|---|
| GTDB bac120/ar53 trees + taxonomy | `references.bac120_tree`, `ar53_tree`, `*_taxonomy` | the taxonomy stage needs only these |
| GTDB-Tk database | `references.gtdbtk_db` | only if you graft unplaced species |
| dbCAN database | `references.dbcan_db` | `run_dbcan database` on the build machine |
| antiSMASH databases | `references.antismash_db` | `download-antismash-databases` |
| AMRFinderPlus database | `references.amrfinder_db` | `amrfinder -u` |

## Deploy (cluster)

```bash
# same-path build
rsync -a /shared/hgn/conda-envs/  node519:/shared/hgn/conda-envs/
rsync -a deploy/                  node519:/path/to/repo/deploy/
ssh node519 'cd /path/to/repo && bash scripts/sh/deploy_envs.sh /shared/hgn/conda-envs'

# packed build
rsync -a deploy/ node519:/path/to/repo/deploy/
ssh node519 'cd /path/to/repo && bash scripts/sh/deploy_envs.sh /shared/hgn/conda-envs --packed deploy/packed'
```

`deploy_envs.sh` checks the tarball checksums, unpacks and `conda-unpack`s where
needed, installs the pip wheels with `--no-index`, runs `verify_envs.sh`, and
then does the one test that matters: `snakemake --sdm conda
--conda-create-envs-only` on the air-gapped host. If any environment is missing
or its hash does not match, that fails *there*, in seconds, rather than three
hours into the first `run_dbcan` job.

## Run

```bash
export PATH="/shared/hgn/conda-envs/launcher/bin:$PATH"
export CONDARC="$PWD/config/condarc"

bash tests/run_tests.sh
python scripts/python/check_sge_profile.py config/sge

snakemake --workflow-profile config/sge --sdm conda \
          --conda-prefix /shared/hgn/conda-envs -j 200
```

`--conda-prefix` is what points Snakemake at the deployed environments; without
it, it looks in `.snakemake/conda` and tries to build them, which on an
air-gapped host means the run dies at the first rule.

## Freezing the solve

Once the environments build, freeze them so the next rebuild is bit-identical:

```bash
conda list -p /shared/hgn/conda-envs/<hash> --explicit > envs/r.linux-64.pin.txt
```

Snakemake uses a `<env>.<platform>.pin.txt` sitting next to the env file in
preference to solving the YAML. Update the pins whenever the YAML changes, or the
two silently diverge.

## Known sharp edges

`envs/scoary.yaml` installs `scoary-2` from PyPI, so that environment cannot be
built without internet. It is built on the build machine like the others; nothing
extra is needed, but it does mean the scoary environment cannot be rebuilt on the
cluster from the YAML alone.

`snakemake-executor-plugin-sge` is not on conda-forge. `build_envs.sh` runs `pip
download` for it and `deploy_envs.sh` installs it with `--no-index`. If that
fails, use `config/sge-generic`, which needs only `qsub`/`qstat`/`qacct`/`qdel`
and the conda-forge `snakemake-executor-plugin-cluster-generic`.

Version pins in `envs/*.yaml` have not been solved anywhere yet. The first
`build_envs.sh` run is where an unsatisfiable pin surfaces; relax the offending
pin and re-run rather than dropping `nodefaults` to make the solve succeed.

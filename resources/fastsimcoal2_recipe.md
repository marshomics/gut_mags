# Optional: coalescent split/migration modelling with fastsimcoal2

The transition module computes the 2D folded joint SFS per niche pair
(`results/transition/work/<species>/demography/joint_sfs_<X>_<Y>.tsv`) and a
directionality call from polymorphism nestedness. For a full split-time, Ne and
directional-migration estimate, fit a coalescent model to that joint SFS with
fastsimcoal2. This is left as a documented manual step because the model
templates must be checked against your sampling, and the binary is not on conda.

This step is optional. The recent-acquisition verdict does not depend on it; it
rests on the rooted ancestral-state polarization, the within-niche diversity and
Tajima's D, the nestedness directionality, and the accessory gene gains. The IM
model adds a calendar/coalescent split time and an explicit migration direction.

## Steps

1. Convert the joint SFS to fastsimcoal2 `.obs` format
   (`_jointMAFpop1_0.obs`); the matrix written by the pipeline is already
   folded, populations X (rows) by Y (columns).

2. Define three competing models as `.tpl` / `.est` pairs:
   - SI  strict isolation (split, two Ne, no migration);
   - IM  isolation with symmetric or asymmetric migration (split, two Ne, m12, m21);
   - AM  ancient migration (migration only before a time point).

3. Run each model with several independent starts:
   ```
   fsc28 -t model.tpl -e model.est -m -0 -C 10 -n 200000 -L 40 -M -q
   ```
   repeat `transition.demography.fsc_runs` times; keep the run with the highest
   likelihood per model.

4. Compare models by AIC (from the best `.bestlhoods` per model); the lowest AIC
   wins. Read off T (split time), N1/N2/NANC and m12/m21. Convert to years with
   `transition.demography.mutation_rate_per_site_per_year` and a generation time.

5. Interpretation for "recent acquisition": a small, recently expanded derived
   Ne, a young T, and migration biased source -> derived agree with the
   nestedness and diversity evidence. Report all models and the AIC table, not
   only the winner.

Set `transition.demography.fastsimcoal2_bin` in the config to record the binary
path you used, for provenance.

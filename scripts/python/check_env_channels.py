#!/usr/bin/env python3
"""
check_env_channels.py
---------------------
Refuse to build an environment that could resolve a package from `defaults` or
`anaconda`.

Those are not neutral mirrors. They carry Anaconda's commercial terms of service,
and they ship different builds of the same version number than conda-forge -- so
an environment that pulls one package from `defaults` is not the environment
anyone else reproducing this analysis will get, even from the same YAML. bioconda
is built and tested against conda-forge with strict priority; anything else
solves by accident.

Every envs/*.yaml must therefore:
  * list conda-forge (bioconda too, where it needs it),
  * list `nodefaults`, which strips the implicit `defaults` channel,
  * not name `defaults` / `anaconda` / `main` / `r` / `pkgs/*` anywhere,
  * not pin a package with a `channel::package` spec pointing at those.

Run standalone, or as step 1 of scripts/sh/build_envs.sh:
    python scripts/python/check_env_channels.py envs/
"""
import glob
import os
import sys

import yaml

BANNED = {"defaults", "anaconda", "main", "r", "free", "pro", "msys2",
          "pkgs/main", "pkgs/r", "pkgs/free", "anaconda/pkgs/main"}
BANNED_PREFIX = ("defaults::", "anaconda::", "main::", "pkgs/")
BANNED_URL = ("repo.anaconda.com", "repo.continuum.io", "anaconda.com/pkgs")


def check(path):
    problems = []
    cfg = yaml.safe_load(open(path)) or {}
    chans = cfg.get("channels") or []
    if isinstance(chans, str):
        chans = [chans]
    lower = [str(c).strip().lower() for c in chans]

    for c in lower:
        if c in BANNED or any(u in c for u in BANNED_URL):
            problems.append(f"channel {c!r} is banned")
    if "nodefaults" not in lower:
        problems.append("channels must include `nodefaults`, or conda appends "
                        "`defaults` implicitly")
    if "conda-forge" not in lower:
        problems.append("channels must include `conda-forge`")
    if lower and lower[0] != "conda-forge":
        problems.append(f"conda-forge must come first for strict priority, got {lower[0]!r}")

    def scan(dep, where):
        s = str(dep).strip().lower()
        if s.startswith(BANNED_PREFIX):
            problems.append(f"{where}: spec {dep!r} names a banned channel")
        if any(u in s for u in BANNED_URL):
            problems.append(f"{where}: spec {dep!r} points at an Anaconda repo URL")

    for dep in cfg.get("dependencies") or []:
        if isinstance(dep, dict):
            for extra in dep.get("pip") or []:
                scan(extra, "pip")
        else:
            scan(dep, "dependencies")
    return problems


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "envs"
    files = sorted(glob.glob(os.path.join(root, "*.yaml")))
    if not files:
        sys.exit(f"no environment files under {root}")
    bad = 0
    for f in files:
        probs = check(f)
        if probs:
            bad += 1
            print(f"{f}:", file=sys.stderr)
            for p in probs:
                print(f"  - {p}", file=sys.stderr)
    if bad:
        print(f"\n{bad}/{len(files)} environment file(s) could resolve from "
              f"defaults/anaconda.", file=sys.stderr)
        return 1
    print(f"OK: {len(files)} environment files are conda-forge/bioconda only, "
          f"with nodefaults.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

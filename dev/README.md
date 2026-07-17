# `dev/` — developer scratch area

## `dev/data/` is for local, private input data only (NOT tracked in Git)

The entire `dev/data/` directory is **deliberately excluded from version control**
(see the repository `.gitignore`: `dev/data/`). No real genotypes, pedigree
records, individual-level outputs, or derived haplotypes may ever be committed to
this repository. `git ls-files dev/data` must return nothing.

## Expected files

Place your local copies under `dev/data/` (or anywhere else on disk, then point the
environment variables below at them):

- a pedigree CSV  — e.g. `ped_HSMap.csv`
- a genotype CSV  — e.g. `geno_HSMap.csv`

These are read with `read_HSMap_data()`.

## How the analysis scripts find the data

The pilot script `analysis/real_data_pilot.R` takes **all** paths from the
environment (or command-line arguments); there are no hard-coded paths. Set, for
example:

```sh
export HSMAP_PEDIGREE=/path/to/ped_HSMap.csv
export HSMAP_GENOTYPES=/path/to/geno_HSMap.csv
export HSMAP_PILOT_OUTPUT=/path/to/output_dir     # created if missing; keep OUTSIDE Git
export HSMAP_THREADS=4                             # optional; safe auto-detect otherwise
export HSMAP_PILOT_N=1500                          # backbone size (300 = smoke test)
```

Keep `HSMAP_PILOT_OUTPUT` outside the repository (or under `analysis/output/`,
which is git-ignored). Do not move private data into a tracked path.

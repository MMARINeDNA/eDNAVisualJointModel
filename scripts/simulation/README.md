# Moved → HSGP4eDNA

The HSGP eDNA **simulation-and-recovery study** that lived in this folder has
moved to its own repository, with full commit history:

### https://github.com/MMARINeDNA/HSGP4eDNA
(local checkout: `Documents/MMARINeDNA/HSGP4eDNA`)

It contains everything that used to be here plus the shared Stan models that
moved with it:

- the **vS1–vS4** simulation scenarios (2-D and 3-D latent fields),
- the vS1 estimation pipeline (`01`–`05`) and diagnostics,
- the HSGP **basis-count sensitivity sweep** and the adaptive basis rule,
- the **3-D bottom-depth demos** (bottom depth as a GP axis vs a parametric covariate),
- the Stan models `whale_edna_hsgp_vS1.stan` and `whale_edna_hsgp_2D_bathycov.stan`,
- shared helpers (`vS1_functions.R`, `helper_functions.R`).

See that repository's `README.md` for the pipeline, how to run it, and results.

The older `stan/whale_edna_hsgp_v1..v4*.stan` models and
`scripts/older_simulations/` remain in this repository.

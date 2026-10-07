# eDNAVisualJointModel

Joint modelling of **visual line-transect surveys** and **environmental DNA
(eDNA)** for marine megafauna off the US West Coast. One latent density
surface per species, a 2-D Hilbert-space Gaussian process (HSGP) over position
plus a spline on bottom depth, feeds two observation models:

```
log λ_s(x) = μ_s + f_s(X, Y) + h_s(Z_bathy)                       # animals / km²
eDNA:   log λ_edna = log λ_s + log zsample_eff + log conv_s + log vol  → qPCR hurdle + ZI-Beta-Binomial MB
visual: λ_groups_s = λ_s / E[group size_s]                          → half-normal detection, Poisson segment counts
```

The repo is being restructured toward **one eDNA model, one visual model, and
one joint model**. See [`ROADMAP.md`](ROADMAP.md) for the plan and the
per-phase checks. Everything removed in the restructure is preserved at the git
tag [`pre-restructure`](https://github.com/MMARINeDNA/eDNAVisualJointModel/tree/pre-restructure)
(the v1–v4.1 eDNA models and pipelines, the `distance/` line-transect models,
their outputs and notebooks).

## Status

| Model | Where | Status |
|---|---|---|
| **eDNA** (2-D HSGP + bottom-depth spline) | [HSGP4eDNA](https://github.com/MMARINeDNA/HSGP4eDNA) submodule, `stan/hsgp_2d_bathysp.stan`; driven by `R/fit_edna.r` | ✅ Validated in HSGP4eDNA (see its README) |
| **Visual** (2-D HSGP + bottom-depth spline, line transects) | `stan/hsgp_visual.stan`; driven by `R/fit_visual.r` | ✅ Validated over 5 replicates (see [Results](#results)) |
| Joint simulator | `R/functions_joint.R` (`simulate_joint()`) | ✅ One shared field; each half validated with its single-source model |
| **Joint** eDNA + visual model | `stan/hsgp_joint.stan` | 🚧 Roadmap Phase 5 |

## Getting started

Clone with the submodule:

```bash
git clone --recurse-submodules https://github.com/MMARINeDNA/eDNAVisualJointModel.git
```

In an existing clone, run `git submodule update --init` after pulling.

Requirements: R with `cmdstanr` + CmdStan, `tidyverse`, `posterior`, `MASS`,
`splines`, `patchwork`, `viridis`, `knitr`. All scripts run from the
**repo root**.

### Smoke test

Run this before committing any change to `R/`, `stan/`, or a submodule bump. It takes a few minutes, and the exit status is non-zero on
failure:

```bash
Rscript tests/smoke.R
```

It fits each working model on a small simulation with short chains and checks
that sampling is clean and the truth is recovered. Logs and results go to
`outputs/smoke/`.

### eDNA model

```r
source("R/fit_edna.r")
res <- fit_edna("bathysp_surface")   # simulate + fit, HSGP4eDNA validation settings
res$field                            # latent-field R² per species
res$recovery                         # GP hyperparameter recovery
```

The model, simulator, and formatter live in `external/HSGP4eDNA`; this repo
does not copy them. Outputs land in `external/HSGP4eDNA/outputs/<scenario>/`.

### Visual model

```r
source("R/fit_visual.r")
sim <- simulate_visual(seed = 1)      # one field + sightings for humpback and PWSD
res <- fit_visual(sim, "humpback")    # 4 chains x 1000/1000
res$recovery                          # density, detection, group size, GP hyperparameters
res$field                             # latent-field R² over segments
```

The model has these parts:

- **Latent field.** It is the eDNA model's field, with the same units
  (animals/km²), GP hyperparameters, spline basis and true spline
  coefficients. It is simulated over 25 E–W transects of 10 km segments.
- **Detection.** Half-normal detection, with a group-size covariate for
  PWSD.
- **Group sizes.** Resampled from the empirical pools in `data/grpsz/`.
- **Model structure.** The detection, population-average ESW and size-bias
  machinery is the former `distance/distance_hn_dens_v4.1.stan`, now in
  `stan/include/visual_functions.stan`. The data switches `use_gp = 0` and
  `K_bathy = 0` reduce the model to that non-spatial version.

As scripts:

```bash
Rscript R/fit_visual.r
```

```bash
SEEDS=1:5 Rscript R/validate_visual.r
```

`R/fit_visual.r` fits one dataset and writes to `outputs/visual/<SCENARIO>/`.
`R/validate_visual.r` runs the multi-replicate recovery study. The development
history of the detection model is in `docs/history/distance_v4.1_notebook.html`.

### Joint simulator

```r
source("R/functions_joint.R")
sim <- simulate_joint(seed = 1)   # one field; eDNA stations + transect segments
sim$edna                          # shaped like HSGP4eDNA's simulate_bathysp() output
sim$visual                        # shaped like simulate_visual() output
```

The joint model is not built yet (roadmap Phase 5). Until then, each half can
be fitted with its single-source model. `Rscript R/validate_joint_sim.r` does
exactly that, and is the Phase 4 check.

## Results

**Visual model** (`R/validate_visual.r`, 5 replicates, 4 chains × 1000/1000,
design-based GP priors, basis 28×14 = 392, spline df 4):

| Species | Detections | Runtime | Divergences | max Rhat | Field R² | Coverage D / σ_det / E[s] | Coverage gp_sigma / lx / ly |
|---|---|---|---|---|---|---|---|
| Humpback | 549–917 | 2.1–4.8 min | 0 in 5 fits | ≤ 1.01 | 0.87–0.93 | 5/5, 5/5, 5/5 | 5/5, 5/5, 5/5 |
| PWSD | 41–100 | 5.2–7.9 min | 1 each in 2 fits | ≤ 1.00 | 0.60–0.79 | 5/5, 4/5, 5/5 | 5/5, 4/5, 5/5 |

Median relative bias:
- mean density: −3% for both species;
- `gp_sigma`: +3% (humpback), +5% (PWSD);
- `lx`: −6% (humpback), +48% (PWSD, which has only 41–100 detections).

The earlier HSGP4eDNA priors gave 2/5 to 5/5 coverage, `gp_sigma` +40% and
PWSD `lx` +99%. See `ROADMAP.md`.

### GP priors

The field priors come from the **survey design only**: sample locations,
domain extent and the basis budget. Observed outcomes never enter
(`R/gp_priors.R`, `gp_prior_from_design()`):

- **Length-scales:** inverse-gamma on each axis, with 1% of mass below ℓ_min
  and 1% above the domain extent. ℓ_min is the larger of two floors:
  - the median spacing of sample coordinates along the axis;
  - the shortest length-scale a basis of `GP_BASIS_BUDGET` = 400 functions
    can represent.

  The basis is sized to ℓ_min, so it can represent everything the prior
  allows.
- **Marginal SD:** `gp_sigma ~ gamma(4.05, 3.37)`, with 1% of mass below 0.25
  and 1% above 3.
- **Floor check:** `floor_check()` flags a fit whose length-scale posterior
  piles up against ℓ_min. This uses no truth, so it works on real data. A flag
  means the budget may be limiting the fit: refit with a larger `M_max` and
  see whether the length-scale moves.

## Layout

```
.
├── ROADMAP.md              restructure plan + per-phase results
├── external/HSGP4eDNA/     submodule: eDNA model, HSGP helpers, simulators
├── R/                      drivers: fit_edna.r, functions_visual.R, functions_joint.R, gp_priors.R,
│                           fit_visual.r, validate_visual.r, validate_joint_sim.r
├── stan/                   hsgp_visual.stan; stan/include/visual_functions.stan
├── tests/smoke.R           fast end-to-end check of every model
├── data/                   real survey data (+ data/grpsz/ group-size pools)
├── scripts/data_figures/   real-data maps and exploration (no model fitting)
├── figures/                outputs of scripts/data_figures/
├── presentations/          collaborator slide decks + presenter notes
├── docs/history/           write-ups of earlier model versions
├── notes/                  correspondence
└── outputs/                generated artifacts (mostly git-ignored)
```

`data/MV1_MURI_df.csv` (MARVER1 metabarcoding, 30 MB) is git-ignored. Get it
from the project share. It is needed only by the eDNA data-figure scripts.

## Data figures

The scripts in `scripts/data_figures/` make summary maps of the **real**
survey data and write them to `figures/`. Run them as
`Rscript scripts/data_figures/<name>.R`.

| Script | Output |
|---|---|
| `plotLTData.R` | `lt_data_pwsd_humpback.png`: LT effort + sightings, PWSD and humpback |
| `ploteDNAData.R` | `edna_data_4panel.png`: hake qPCR; hake / PWSD / humpback MARVER1 |
| `plotPWSDData.R` | `lt_edna_pwsd.png`: PWSD LT + eDNA side by side |
| `plotHumpbackData.R` | `lt_edna_humpback.png`: humpback LT + eDNA side by side |
| `plotHarborPorpoiseData.R` | `lt_edna_harbor_porpoise.png`: harbor porpoise LT + eDNA |
| `plotHarborPorpoiseeDNA.R` | `edna_harbor_porpoise.png`: harbor porpoise eDNA only |
| `getGrpSz.R` | `data/grpsz/{humpback,pwsd}.rds`: empirical group-size pools used by the visual simulator |

`FinWhales.R` and `DetectionsBySpecies.R` are ad-hoc exploration.

## Presentations and history

- `presentations/` holds Quarto / RevealJS decks with companion
  presenter-notes documents: an introduction to GPs and HSGPs, and the
  project goals. Render them with `quarto render presentations/<file>.qmd`.
  The presenter notes still refer to the pre-restructure file layout, which
  will be refreshed in roadmap Phase 6.
- `docs/history/` has two self-contained write-ups. `v3.2_notebook.html` is
  the eDNA sampler-pathology debugging case study, and
  `distance_v4.1_notebook.html` covers the distance-sampling development
  (PRs #34–#37).

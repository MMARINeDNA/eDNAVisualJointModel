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
(the v1–v4.1 eDNA models and pipelines, their outputs and notebooks).

## Status

| Model | Where | Status |
|---|---|---|
| **eDNA** (2-D HSGP + bottom-depth spline) | [HSGP4eDNA](https://github.com/MMARINeDNA/HSGP4eDNA) submodule, `stan/hsgp_2d_bathysp.stan`; driven by `R/fit_edna.r` | ✅ Validated in HSGP4eDNA (see its README) |
| **Visual**, non-spatial | `distance/00_distance_v4.1.R`, `distance/distance_hn_dens_v4.1.stan` | ✅ Validated; the detection-function reference |
| **Visual**, spatial (2-D HSGP + spline) | `stan/hsgp_visual.stan` | 🚧 Roadmap Phase 3 |
| Visual, spatial v4.1a (old 3-D HSGP) | `distance/00_distance_v4.1a.R`, `distance/distance_v4.1a.stan` | ❌ Broken (HSGP weights missing the coordinate-scale Jacobian; see `ROADMAP.md`). Kept only as a template until Phase 3 replaces it. |
| **Joint** eDNA + visual | `stan/hsgp_joint.stan` | 🚧 Roadmap Phases 4–5 |

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

Run this before committing any change to `R/`, `stan/`, `distance/`, or a
submodule bump. It takes a few minutes, and the exit status is non-zero on
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

### Visual model (non-spatial)

```bash
Rscript distance/00_distance_v4.1.R
```

For each cetacean species (humpback, Pacific white-sided dolphin) this script
does three things:

- It simulates line-transect sightings from a GP density surface, using a
  half-normal detection function. For PWSD the detection function has a
  group-size covariate.
- It fits the model with a population-average ESW and a size-bias correction.
- It writes the recovery table and plots to `outputs/distance_v4.1/`.

The env overrides are `OUTPUT_DIR`, `CHAINS`, `ITER_WARMUP` and
`ITER_SAMPLING`. The development history is in
`docs/history/distance_v4.1_notebook.html`.

## Layout

```
.
├── ROADMAP.md              restructure plan + per-phase results
├── external/HSGP4eDNA/     submodule: eDNA model, HSGP helpers, simulators
├── R/                      drivers (fit_edna.r; visual + joint to come)
├── tests/smoke.R           fast end-to-end check of every model
├── distance/               current visual models (replaced in Phase 3)
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
| `getGrpSz.R` | `data/grpsz/{humpback,pwsd}.rds`: empirical group-size pools used by `distance/` |

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

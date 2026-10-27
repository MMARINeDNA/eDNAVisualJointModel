# eDNAVisualJointModel

Joint modelling of **visual line-transect surveys** and **environmental DNA
(eDNA)** for marine megafauna off the US West Coast. One latent density
surface per species, a 2-D Hilbert-space Gaussian process (HSGP) over position
plus a spline on bottom depth, feeds two observation models:

```
log λ_s(x) = μ_s + f_s(X, Y) + h_s(Z_bathy)                       # animals / km²
eDNA:   log λ_edna = log λ_s + log zsample_eff + log conv_s + log vol  → qPCR hurdle + ZI-Beta-Binomial MB
visual: λ_groups_s = λ_s / E[group size_s]                          → half-normal or hazard-rate detection,
                                                                      Poisson or negative-binomial segment counts
```

The repo is being restructured toward **one eDNA model, one visual model, and
one joint model**. See [`ROADMAP_2026-10-21.md`](ROADMAP_2026-10-21.md) for the
planned work, and [`ROADMAP_2026-10-06.md`](ROADMAP_2026-10-06.md) (archived) for
the restructure plan and its per-phase checks. Everything removed in the restructure is preserved at the git
tag [`pre-restructure`](https://github.com/MMARINeDNA/eDNAVisualJointModel/tree/pre-restructure)
(the v1–v4.1 eDNA models and pipelines, the `distance/` line-transect models,
their outputs and notebooks).

## Status

| Model | Where | Status |
|---|---|---|
| **eDNA** (2-D HSGP + bottom-depth spline) | [HSGP4eDNA](https://github.com/MMARINeDNA/HSGP4eDNA) submodule, `stan/hsgp_2d_bathysp.stan`; driven by `R/fit_edna.r` | ✅ Validated in HSGP4eDNA (see its README) |
| **Visual** (2-D HSGP + bottom-depth spline, line transects) | `stan/hsgp_visual.stan`; driven by `R/fit_visual.r` | ✅ Validated over 5 replicates (see [Results](#results)) |
| Joint simulator | `R/functions_joint.R` (`simulate_joint()`) | ✅ One shared field; each half validated with its single-source model |
| **Joint** eDNA + visual model | `stan/hsgp_joint.stan`; driven by `R/fit_joint.r` | ✅ One shared field; beats the best single source for the cetaceans and recovers their eDNA conversion factors (see [Results](#results)) |
| Real data | `R/01_load_muri.r`; real fits `R/fit_visual_real.r`, `R/fit_joint_real.r` | Visual-only fits done (35°N, Beaufort ≤ 5); first real eDNA and joint fits running on stats3 (see [Real data](#real-data)) |

## Study domain

Every analysis, simulated or real, is restricted to the MURI study area
**north of 35°N**. This is HSGP4eDNA's `STUDY_DOMAIN`: a 690 × 1770 km box in
UTM zone 10N, with its south-west corner at (E, N) = (−10, 3870) km and its
southern edge on 35°N.

- **Real eDNA:** 513 bottles at 167 stations (37.9–48.4°N), all inside it.
- **Real visual effort:** only segments whose midpoint is at or north of 35°N
  and whose mean Beaufort is at most 5 are used. That is 1,617 of 2,779
  segments, or 7,441 of 12,967 km.
- **Simulations:** the eDNA, visual and joint simulators all use the same box,
  so simulated and real fits share coordinates and prior extents.

> **Runs in progress (2026-10-26).** The current set-up is the 35°N domain,
> 50 × 150 km length-scales, 35 transects, an 800-function budget, HSGP4eDNA
> v1.12, and NB segment counts with non-centred junk noise.
> - **On stats3 now:** the seed-1 demonstration of the idealised simulation
>   with sea state (its joint and free-conversion fits; the eDNA-only and
>   visual-only fits are done), and the first real-data fits
>   (`R/fit_joint_real.r`: eDNA-only, visual-only, joint with free conversion
>   factors).
> - **Later this week:** the replicates, seeds 2–5.
> - **The results below** are from the earlier runs (v1.5, Poisson counts;
>   2026-10-12 and 2026-10-13) and will be replaced. HSGP4eDNA's own
>   validation is still from its earlier set-up (500 × 1270 km box, 50 × 300 km
>   length-scales, 400-function budget).
>
> The runs use the stats3 server (`notes/stats3_guide.md`); ROADMAP_2026-10-21,
> section 3, lists them.

## Simulated species

The simulations use three species, chosen to span how informative each data
source is:

| Species | eDNA | Visual line transects | Baseline density exp(μ) | GP σ | Visual detection and group size |
|---|---|---|---|---|---|
| **Pacific hake** (*Merluccius productus*) | qPCR (detection + Ct) and metabarcoding | not surveyed | 15 /km² | 2.0 | — |
| **Humpback whale** (*Megaptera novaeangliae*) | metabarcoding only | yes | 0.01 /km² | 1.0 | half-normal σ = 2.5 km; groups of 1–8 (mean 1.7) |
| **Pacific white-sided dolphin** (*Lagenorhynchus obliquidens*, PWSD) | metabarcoding only | yes | 0.02 /km² | 1.3 | σ = 1.5 km at the mean group size, wider for larger groups; groups of 2–452 (mean 43) |

- **Roles.**
  - **Hake** is the abundant, eDNA-rich anchor. It is the only species with
    qPCR, and it is not observed visually.
  - **Humpback** is rare in eDNA but conspicuous visually: about 400–1,000
    detections per simulated survey.
  - **PWSD** is weak in both sources: about 40–120 visual detections. It is
    where the joint model should help most.
- **Shared field structure.** Every species has an anisotropic GP with
  length-scales of 50 km (east–west) and 150 km (along the coast), plus a
  smooth bottom-depth effect. The depth shapes (`BATHY_BETA_TRUE` in
  HSGP4eDNA) are arbitrary test curves, not ecological estimates. Density is in
  animals/km²; for hake the unit is nominal, since hake is seen only through
  eDNA.
- **eDNA.** Expected bottle copies = density × conversion factor (hake 10,
  humpback 200, PWSD 110) × volume filtered, with negative-binomial
  bottle-to-bottle noise (θ = 10). Reads also include a "junk" category for
  all non-target DNA, which varies between samples. The whales sit near the
  detection limit: about 0.04–0.08 expected copies in a typical aliquot, so
  only about 10% of aliquots contain any whale DNA.
- **Water column** (`bathysp_depthspref` only). Bottles are taken at 6
  depths (0–500 m). Hake and humpback eDNA peak below the surface (~150 m);
  PWSD eDNA is surface-associated.
- **Visual surveys.** 35 east–west transects (about 51 km apart) across the study domain, cut into 10 km segments, truncation
  distance 5 km. Group sizes are resampled from the empirical pools in
  `data/grpsz/`. Options: with `sea_state = TRUE`, detection also declines
  with Beaufort sea state (β_bft = −0.08 humpback, −0.15 PWSD per Beaufort
  unit); with `visibility = TRUE` it rises with visibility (β_vis = 0.15,
  0.2 per nmi); `nb_size` makes sightings clustered (negative-binomial groups
  per segment); `segments` simulates on a given design, such as the real
  effort.
- **eDNA options** (`simulate_joint(edna_obs_args = ...)`, HSGP4eDNA's
  simulator): per-plate qPCR curves with standards, metabarcoding dilutions,
  several sampling depths.

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

Run this before committing any change to `R/`, `stan/`, or a submodule bump. It takes about 1.5 hours (44 checks; the joint cases are the longest), and the exit status is non-zero on
failure:

```bash
Rscript tests/smoke.R
```

It fits each working model on a small simulation with short chains and checks
that sampling is clean and the truth is recovered. Logs and results go to
`outputs/smoke/`. The exactness tests (`tests/joint_exactness.R`,
`tests/visual_nb.R`) and the detection test (`tests/visual_detection.R`) check
that each model option is exact when switched off and recovers its truth.

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
  coefficients. It is simulated over 35 E–W transects of 10 km segments.
- **Detection.** Half-normal or hazard-rate (`det_key`), with σ on the log
  scale as in Distance: log σ = log σ₀ + β_size (s − s̄) + β_bft (Beaufort − 3)
  + β_vis (visibility − 6), each covariate switched per species (`use_size_covar`,
  `use_bft`, `use_vis`). Segment Beaufort and visibility are rounded to 0.5,
  and the size-averaged ESW is computed once per detection-condition cell
  (Beaufort × visibility), so ESW differs between segments.
- **Segment counts.** Poisson, or negative binomial with size φ (`use_nb`;
  1/√φ ~ half-normal(0, 1), bounded so φ ≤ 10⁴) for clustered sightings.
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

```bash
SEA_STATE=1 SEEDS=1:3 Rscript R/validate_visual.r
```

`R/fit_visual.r` fits one dataset and writes to `outputs/visual/<SCENARIO>/`.
`R/validate_visual.r` runs the multi-replicate recovery study; with
`SEA_STATE=1` it simulates sea state and fits each species with and without
the covariate. The development
history of the detection model is in `docs/history/distance_v4.1_notebook.html`.

### Joint simulator

```r
source("R/functions_joint.R")
sim <- simulate_joint(seed = 1)   # one field; eDNA stations + transect segments
sim$edna                          # shaped like HSGP4eDNA's simulate_bathysp() output
sim$visual                        # shaped like simulate_visual() output
```

Each half can also be fitted with its single-source model.
`Rscript R/validate_joint_sim.r` does exactly that, and is the Phase 4 check.

### Joint model

```r
source("R/fit_joint.r")
sim <- simulate_joint(seed = 1)
res <- fit_joint(sim, "joint")        # also "edna" or "visual"; 4 chains x 800/800
res$field                             # R², 95% CI width and coverage per species
res$recovery                          # GP hyperparameters, detection scale, noise theta
fit_joint(sim, "joint", conv_free = c("humpback", "pwsd"))   # estimate conversion factors
```

`stan/hsgp_joint.stan` gives each species one field: a 2-D HSGP plus a
bottom-depth spline, in animals/km². Two observation models share it:

- **eDNA:** HSGP4eDNA's eDNA observation model (the submodule is at v1.12):
  qPCR for hake, zero-inflated beta-binomial metabarcoding for all species,
  and the between-bottle and junk noise terms. The metabarcoding assumes
  equal amplification efficiency across species, which holds for MarVer1 and
  these three species (ROADMAP_2026-10-21, section 5). Ported from
  HSGP4eDNA v1.8–v1.12, each off by default and exact when off: dilution
  offsets for both assays and filtered volume per bottle (used whenever the
  eDNA object carries them), per-plate qPCR standard curves estimated from
  the standards (`qpcr_curve = "estimated"`), a per-species sampling-depth
  smooth (`zdepth = "smooth"`) and the bounded Ct-noise SD
  (`ct_sd = "conditional"`), all `format_stan_data_joint()` arguments. With
  the visual likelihood off and these on, the joint model equals
  HSGP4eDNA's own model (`tests/joint_exactness.R`). Not ported: the
  experimental mechanistic metabarcoding overdispersion.
- **Visual:** the line-transect model for humpback and PWSD. Its half-normal
  or hazard-rate key and per-species truncation carry through, and per
  species the Beaufort and visibility covariates (`bft_species`,
  `vis_species`) and negative-binomial counts (`nb_species`) can be switched
  on. With each off the log density is identical to the model without it, and
  with eDNA off it equals the single-species visual models
  (`tests/joint_exactness.R`).

Data switches turn the eDNA or the visual likelihood off, and free the
cetaceans' conversion factors. So the eDNA-only, visual-only and joint fits use
the same priors, basis and starting values, and differ only in the data they
use.

Exactness checks (`tests/joint_exactness.R`):
- **Every option off:** the log density equals the joint model before the
  NB, visibility and eDNA options, exactly.
- **eDNA off:** it equals the sum of the two `hsgp_visual.stan` fits up to a
  constant (PWSD with Beaufort, visibility and NB counts).
- **Visual off, eDNA options on:** it equals HSGP4eDNA's own
  `hsgp_2d_bathysp.stan` up to a constant.

The spread is ≤ 3e-11 over random draws.

`Rscript R/compare_joint.r` runs all four fits on one simulated dataset and
writes the tables and the figure below to `outputs/joint/compare_s<SEED>/`
(`compare_sea_s<SEED>/` with `SEA_STATE=1`). Its defaults are NB counts for
both cetaceans, non-centred junk noise and `adapt_delta = 0.9`.

## Results

**Visual model** (`R/validate_visual.r`, 5 replicates on the 35°N domain, 35
transects / 2,415 segments, 4 chains × 1000/1000, design-based GP priors,
800-function budget, spline df 4; run 2026-10-12):

| Species | Detections | Runtime | Divergences | max Rhat | Field R² | Coverage D / σ_det / E[s] | Coverage gp_sigma / lx / ly |
|---|---|---|---|---|---|---|---|
| Humpback | 808–1529 | 10–18 min | 0 | ≤ 1.010 | 0.84–0.90 | 5/5, 5/5, 5/5 | 5/5, 4/5, 5/5 |
| PWSD | 94–110 | 20–37 min | 0 | ≤ 1.005 | 0.50–0.69 | 5/5, 5/5, 5/5 | 5/5, 5/5, 4/5 |

Median relative bias:
- mean density: +2% for both species;
- `gp_sigma`: −5% (humpback), −1% (PWSD);
- `lx` / `ly`: −1% / +3% (humpback), +15% / +48% (PWSD, which has only
  94–110 detections).

The north–south length-scale is the one to watch: the floor check flags it in
4 of 5 humpback fits and 2 of 5 PWSD fits (posterior near the prior's
117 km floor against a truth of 150 km). That is the narrow 1.28× visual
north–south margin recorded in ROADMAP_2026-10-06, 6.0; the truth stays covered.

**Sea-state covariate** (`SEA_STATE=1 SEEDS=1:3`; truth β_bft = −0.08
humpback, −0.15 PWSD). Each species fitted with and without the covariate.
With the switch off the model's log density matches the single-ESW model
exactly.

| Species | Covariate | Divergences | max Rhat | Field R² | β_bft covered | Field-error slope per Beaufort unit | Runtime |
|---|---|---|---|---|---|---|---|
| Humpback | with | 0 | ≤ 1.002 | 0.84–0.91 | 3/3 | +0.005 to +0.019 | 9–13 min |
| Humpback | without | 0 | ≤ 1.004 | 0.84–0.91 | — | −0.041 to −0.024 | 8–13 min |
| PWSD | with | 0 | ≤ 1.005 | 0.52–0.66 | 3/3 | +0.075 to +0.100 | 26–35 min |
| PWSD | without | 0 | ≤ 1.004 | 0.52–0.66 | — | −0.034 to +0.064 | 25–29 min |

Leaving the covariate out reads poor sighting conditions as low density: the
fitted field shifts by −0.03 to −0.06 (humpback) and −0.04 to −0.11 (PWSD)
log units per Beaufort unit, up to 0.61 log units across the survey. Density
is covered in every fit either way. β_bft is covered but shrunk toward zero
(median relative bias −29% humpback, −37% PWSD).

**Joint model** (`R/compare_joint.r`, `simulate_joint()` seed 1 on the 35°N
domain, 4 chains × 800/800, data-based starting field; run 2026-10-13; the
replicates are being run with sea state, see the explainer's Section 5). The
same model and priors are used in each column; only the
data differ. Each cell is the field R² over all 200 stations and 2,415
segments, with the mean 95% CI width of log density in brackets. Coverage of
the true field is 0.92–0.96 throughout.

| Species | eDNA only | Visual only | Joint | Joint, conversion free |
|---|---|---|---|---|
| Hake | 0.939 (2.60) | — | 0.939 (2.65) | 0.939 (2.64) |
| Humpback | 0.494 (4.17) | 0.925 (1.32) | **0.928 (1.28)** | 0.929 (1.28) |
| PWSD | 0.589 (4.12) | 0.589 (3.17) | **0.745 (2.78)** | 0.744 (2.81) |

- **Hake** is unchanged from eDNA-only, as it should be: it has no visual data.
- **Humpback** (1,690 sightings) is carried by the visual data. Joint is
  slightly better than visual-only in R² and CI width, and much better than
  eDNA-only.
- **PWSD** gains most: each source alone gives R² 0.59 (153 sightings;
  metabarcoding reads in about one replicate in five), and together they give
  0.75 with narrower CIs than either.
- **Conversion factors.** Freed with a deliberately off-centre prior
  (normal(log 100, 2)), the joint model recovers them: humpback
  5.30 [5.05, 5.53] (truth 5.30), PWSD 4.77 [4.45, 5.08] (truth 4.70). eDNA
  alone cannot identify them, because only density × conversion enters its
  likelihood. They assume equal amplification efficiency across species in
  the metabarcoding, which holds for MarVer1 and these three species
  (ROADMAP_2026-10-21, section 5).
- **Convergence.** 0 divergences in every fit; max Rhat 1.007 (visual), 1.019
  (joint), 1.06 (conversion free) and 1.11 for eDNA-only, on θ_bottle, the
  slowly mixing noise term.
- **Hyperparameters** are mostly covered. The misses: humpback's north–south
  length-scale reads 185–187 km [156, 219] against 150 in the visual and joint
  fits, and hake's and humpback's GP SDs sit slightly high.
- **Starting values.** From a flat starting field this dataset's eDNA side got
  stuck at its extreme hake hot spots (hake R² 0.19, θ_bottle 1.6 against 10).
  Joint fits now start the field from the data (`edna_field_start()`), which
  recovers it (ROADMAP_2026-10-06, 6.0). Joint fits also run 4 chains, so a stuck chain
  shows in Rhat.
- **Runtime** on a laptop: eDNA-only, joint and conversion-free about 2 h each;
  visual-only about 45 min.
- **Earlier set-up.** On the old 500 × 1270 km box (Phase 5, seeds 1 and 2,
  3 chains) the pattern was the same: joint best for both cetaceans, hake
  unchanged, conversion factors recovered. One of those runs did not converge
  (one chain with very large junk noise).

![Joint vs single-source field recovery, seed 1 on the 35°N domain: posterior mean vs true log density per species (rows) and fit (columns)](docs/figures/joint_compare_s1.png)

## Real data

The MURI data: the 2019 hake survey's eDNA (qPCR and MarVer1 metabarcoding)
and the 2018 CCES line-transect survey, north of 35°N.

- **Loading.** `R/01_load_muri.r` builds the eDNA and visual objects in the
  simulators' shape (`data/processed/muri_{edna,visual}.rds`, tracked);
  `R/02_check_muri.r` checks them. `R/fit_qpcr_standards.r` fits the hake
  standard curves per plate (`data/processed/qpcr_standard_curves.rds`). The
  choices (target depths, dilutions, read cutoff, contamination, Beaufort ≤ 5,
  visibility, NB counts) are in `notes/muri_data_formatting.md`.
- **Visual fits** (`R/fit_visual_real.r`): humpback 0.023 and PWSD 0.108
  animals per km² on the Beaufort ≤ 5 effort; PWSD detection improves with
  visibility (0.17 per nmi on log σ). The length-scales sit at the basis floor
  (about 40 × 65 km, the Riutort-Mayol limit for the 800-function basis), and
  followed it down when the basis was doubled. Keeping budget 800 is a
  computational compromise (density does not change with the budget).
  Clustering independent between segments is ruled out as the cause;
  clustering correlated over several segments, or real fine-scale density,
  is not (`notes/real_visual_fits_2026-10-22.md`,
  `notes/count_dispersion_check.md`; tests planned in ROADMAP_2026-10-21,
  section 6).
- **eDNA and joint fits** (`R/fit_joint_real.r`, running on stats3): Ct with
  per-plate curves estimated from the 659 hake standards, dilution offsets,
  per-bottle volumes and a sampling-depth smooth; the joint fit estimates the
  cetaceans' conversion factors.

### GP priors

The field priors come from the **survey design only**: sample locations,
domain extent and the basis budget. Observed outcomes never enter
(`gp_prior_from_design()` in HSGP4eDNA's `R/functions.R`):

- **Length-scales:** inverse-gamma on each axis, with 1% of mass below ℓ_min
  and 1% above the domain extent. ℓ_min is the larger of two floors:
  - the median spacing of sample coordinates along the axis;
  - the shortest length-scale a basis of `GP_BASIS_BUDGET` = 800 functions
    can represent (400 before HSGP4eDNA v1.4).

  The basis is sized to ℓ_min, so it can represent everything the prior
  allows.
- **Marginal SD:** `gp_sigma ~ gamma(2.94, 1.66)`, with 1% of mass below 0.25
  and 1% above 5 (raised from 3 in HSGP4eDNA v1.11, as the real visual fits
  sat near the old bound).
- **Floor check:** `floor_check()` flags a fit whose length-scale posterior
  piles up against ℓ_min. This uses no truth, so it works on real data. A flag
  means the budget may be limiting the fit: refit with a larger `M_max` and
  see whether the length-scale moves.

## Layout

```
.
├── ROADMAP_2026-10-21.md   planned work (current)
├── ROADMAP_2026-10-06.md   restructure plan + per-phase results (archived)
├── external/HSGP4eDNA/     submodule: eDNA model, HSGP helpers, simulators
├── R/                      drivers: fit_edna.r, functions_visual.R, functions_joint.R,
│                           fit_visual.r, validate_visual.r, validate_joint_sim.r,
│                           fit_joint.r, compare_joint.r, check_count_dispersion.r;
│                           real data: 01_load_muri.r, 02_check_muri.r, fit_qpcr_standards.r,
│                           fit_visual_real.r, fit_joint_real.r
├── stan/                   hsgp_visual.stan, hsgp_joint.stan, qpcr_standards.stan;
│                           stan/include/visual_functions.stan
├── tests/                  smoke.R (every model), joint_exactness.R, visual_nb.R,
│                           visual_detection.R
├── tools/stats3_reruns.sh  runs the batches on the stats3 server
├── data/                   real survey data (+ data/grpsz/ group-size pools)
├── scripts/data_figures/   real-data maps and exploration (no model fitting)
├── figures/                outputs of scripts/data_figures/
├── presentations/          collaborator slide decks + presenter notes
├── docs/history/           write-ups of earlier model versions
├── notes/                  data-formatting notes, analysis write-ups, the stats3 guide
└── outputs/                generated artifacts (mostly git-ignored)
```

The raw MURI files are private and git-ignored: `qPCR_2019A.rds`,
`M3_compiled_taxon_table_wide.csv`, `Hake_2019_metadata.csv` and
`SITEINFO_CCE_2018_Pp_Lo_Mn.csv` (read by `R/01_load_muri.r`), and the older
`MV1_MURI_df.csv` (used only by the eDNA data-figure scripts). Get them from the
project share.

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
- `presentations/smm2026.qmd` is the 15-minute SMM 2026 talk (goals, real
  data, model structure, simulation results). It renders
  to HTML and PDF, with speaker notes inline. Its figures are built by
  `Rscript presentations/smm2026_figs/make_figures.R` from the processed
  real-data objects and `docs/figures/model_*.svg`.
- `docs/explainer.md` (rendered `docs/explainer.html`) walks through every model
  with simulated and real data, prior and posterior predictive checks, and the
  open issues. Its first sections predate the 35°N study domain (see its status
  notes); later sections are updated to 2026-10-22. It will be brought up to
  date after the current stats3 runs (ROADMAP_2026-10-21, section 3). Figure scripts are in `docs/explainer/`; re-render with
  `quarto render docs/explainer.md --to html --embed-resources`.
- `docs/history/` has two self-contained write-ups. `v3.2_notebook.html` is
  the eDNA sampler-pathology debugging case study, and
  `distance_v4.1_notebook.html` covers the distance-sampling development
  (PRs #34–#37).

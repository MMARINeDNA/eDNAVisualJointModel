# Whale eDNA HSGP simulation study

Simulation-and-recovery scripts for the joint qPCR + metabarcoding eDNA model
with a Hilbert-space approximate Gaussian process (**HSGP**) latent field. Each
scenario simulates eDNA data from a known truth, fits the Stan model, and
checks how well the field and parameters are recovered.

The model links two observation processes to one shared latent animal-density
field per species:

```
log(lambda_si)      = mu_s + f_s(coords_i)          # HSGP latent field
log(lambda_edna_si) = log(lambda_si) + log_zsample_effect + log_conv_factor[s] + log_vol_filtered
  -> qPCR hurdle (Pacific hake only) + zero-inflated Beta-Binomial metabarcoding (all species + a "junk" background)
```

`f_s` is a per-species anisotropic HSGP. The Stan model
(`stan/whale_edna_hsgp_vS1.stan`) is **dimension-general**: the number of GP
axes is set by the data (`D1`, `INDICES`, `m_hsgp`), so the same model serves
1-D, 2-D, or 3-D fields.

---

## Quick start (vS1 reference pipeline)

Run from the repository root, in order:

```r
source("scripts/simulation/01_simulate_whale_edna_vS1.r")   # simulate -> outputs/whale_edna_output_vS1/
source("scripts/simulation/02_plot_simulated_data_vS1.r")   # simulated-data summary PDF
source("scripts/simulation/03_format_stan_data_vS1.r")      # build Stan data (+ pred_points)
source("scripts/simulation/04_run_whale_edna_model_vS1.r")  # compile + sample (cmdstanr)
source("scripts/simulation/05_check_whale_edna_model_vS1.r") # diagnostics + recovery PDF
```

Each script is **self-contained** (loads what it needs from disk) and runs in a
fresh session. **cmdstanr is the canonical sampler**; post-processing uses the
cmdstanr / `posterior` API. Requires `cmdstanr` + a working CmdStan, plus
`tidyverse`, `posterior`, `bayesplot`, `splines`.

---

## Scenarios

All scenarios share the observation model, the "junk" background category, and
a zero-mean GP truth. They differ in the **latent-field dimensions** and the
**water-column sampling design**:

| Scenario | GP field | Sampling | Stations | N samples | Pipeline status |
|----------|----------|----------|----------|-----------|-----------------|
| **vS1**  | 2-D (X, Y)              | surface only        | 250 | 250 | **Full** (01–05) + sweep + validated |
| **vS2**  | 3-D (X, Y, Z_bathy)     | surface only        | 200 | 200 | sim only (01); 3-D demo below |
| **vS3**  | 3-D (X, Y, Z_bathy)     | 6 depths, no pref   | 200 | 831 | sim only (01) |
| **vS4**  | 3-D (X, Y, Z_bathy)     | 6 depths + depth pref | 200 | 831 | sim only (01) |

- **Z_bathy** = bottom depth at a station (a continuous GP covariate, ~66–3100 m).
- **Z_sample** = the water-column depth a sample was filtered at; it drives the
  `zsample_effect` offset (vS4 adds a depth preference), **not** the GP.
- vS3/vS4 drop samples where `Z_sample > Z_bathy`, so N < stations × depths.

Only **vS1** has the full estimation pipeline (02–05). vS2–4 currently ship the
simulation stage only; a 3-D / multi-depth estimation pipeline is demonstrated
(not productionised) by the two demos below.

---

## Files

| File | Role |
|------|------|
| `01_simulate_whale_edna_vS{1,2,3,4}.r` | Simulate a scenario → `outputs/whale_edna_output_vS{n}/` |
| `02_plot_simulated_data_vS1.r` | Simulated-data summary plots (vS1) |
| `03_format_stan_data_vS1.r` | Assemble the Stan data list + prediction grid (vS1) |
| `04_run_whale_edna_model_vS1.r` | Compile + sample with cmdstanr; save fit bundle |
| `05_check_whale_edna_model_vS1.r` | Convergence diagnostics + parameter/field recovery |
| `06_basis_sweep_vS1.r` | Systematic HSGP basis-count sensitivity sweep |
| `07_plot_basis_sweep_vS1.r` | Summarise/plot the sweep + recommendation CSV |
| `demo_A_hsgp3D_vS2.r` | Demo A: 3-D HSGP (X, Y, Z_bathy) on vS2 |
| `demo_B_hsgp2D_bathycov_vS2.r` | Demo B: 2-D HSGP + parametric bottom-depth spline |
| `vS1_functions.R` | Shared callable functions (see below) |

**Stan models** (in `stan/`):
- `whale_edna_hsgp_vS1.stan` — dimension-general HSGP + qPCR/MB observation model.
- `whale_edna_hsgp_2D_bathycov.stan` — 2-D HSGP + a parametric spline of bottom depth (demo B).

**`vS1_functions.R`** provides: `simulate_whale_edna_vS1()` /
`simulate_whale_edna_vS2B()` (seed-parameterised simulators),
`format_stan_data_vS1()` / `format_stan_data_hsgp()` (n-D Stan-data formatters),
`hsgp_basis_rule()` (a-priori basis rule), `bathy_spline_basis()` (centred
natural spline), `format_stan_data_bathycov()`. Used by the sweep and demos;
the numbered vS1 scripts are standalone and don't depend on it.

---

## Key findings

### 1. Junk read-proportion fix (all scenarios)

The metabarcoding "junk" background must be **constant within a sample** and the
target read proportions must reflect the **aliquot-level** draws (`mb_copies`),
not the bottle-level expected copies — otherwise the metabarcoding reads and the
qPCR aliquot draws imply different realisations (a qPCR/MB disagreement). vS1
additionally had a recycling misalignment (a per-sample junk vector added to an
aliquot-long-form vector); fixed by expanding via `mb_sample_idx`. All four
sims now drive the multinomial proportions from `mb_copies` / `mb_junk_copies`.

### 2. HSGP basis-count sweep (vS1)

`06`/`07` sweep the basis count over under- → adequate → over-resolved, across
simulation seeds. Latent-field recovery (posterior-mean `log_lambda` vs truth):

| Config | M | field R² (hake) | `lx` (true 50) | reading |
|--------|---|-----------------|----------------|---------|
| 4×4    | 16  | 0.63 | 112 | oversmoothed |
| 8×4    | 32  | 0.89 | 58  | good; `lx` CI covers |
| **13×6** | **78** | **0.98** | 44 | the rule — plateau |
| **16×8** | **128** | **0.98** | 42 | **default** — on plateau |
| 24×12  | 288 | 0.96 | 25 | overfit (`lx` collapsing) |
| 32×16  | 512 | 0.95 | 20 | overfit |

- **Sweet spot ≈ M 32–128.** R² plateaus at ~0.98 then *declines* as M grows.
- **Overfit signature:** past M≈128, `lx` collapses toward 0.
- **Anisotropy binds on the short axis:** at equal M=64, spending basis on the
  short cross-shore dimension (16×4) gives R²=0.96 vs 0.71 for 4×16.
- The default `HSGP_M = c(16, 8)` (`03`) sits safely on the plateau.

### 3. Adaptive basis rule

`hsgp_basis_rule(lengthscale, half_range, c=1.5)` returns the Riutort-Mayol
faithful-representation minimum `m ≥ 1.75·c/ρ` per dimension (ρ = length scale /
domain half-range). It reproduces the sweep's (14, 6) for vS1 and, sized to the
shortest length scale across species, gives more basis to the short axis
automatically. `format_stan_data_vS1(HSGP_M = "auto")` derives it from the truth.

**3-D caveat:** in vS2–4 the rule wants `(17, 6, 40) → M ≈ 4080`. The Z_bathy
axis alone needs ~40 basis (short bottom-depth length scale over a large range),
which is computationally prohibitive — motivating the two depth demos.

### 4. Handling bottom depth in 3-D: two demos (vS2)

| | **Demo A — 3-D HSGP** | **Demo B — 2-D HSGP + spline covariate** |
|--|----------------------|------------------------------------------|
| Bottom depth | 3rd GP axis (Z basis capped at 20) | centred natural spline in `log_lambda` |
| Basis M | 1680 | **84** (~20× cheaper) |
| Field R² (hake) | 0.83 | **0.97** |
| Depth recovered | `lz` = 147 (true 150) ✓ | spline curve ✓ |
| GP length scales | `lx` ✓, `ly` under, **Rhat ≤ 1.6** | `lx`/`ly`/`gp_sigma` ✓, clean |
| `gp_sigma` | inflated (Z cap) | ✓ |
| Divergences | 0 | 0 |

Both incorporate bottom depth. **Treating it as a parametric covariate (B) is
the better default** — far cheaper, mixes cleanly, and recovers the field, GP
length scales, and the depth curve. The **full 3-D GP (A)** is feasible and
recovers the field but is expensive and harder to sample. (A and B use different
matched truths, so their R² are each vs their own truth — not a head-to-head.)

---

## Outputs

Each scenario writes to `outputs/whale_edna_output_vS{n}/` (git-ignored; the
`simulated_edna_fields_v*.pdf` plots are the exception and are tracked):

- `whale_edna_sim_vS{n}.rds` — the simulation (truth + observed data)
- `stan_data.rds`, `pred_points.rds` — Stan inputs (vS1)
- `whale_edna_vS1_fit.rds` — fit bundle (draws, summaries, prediction data frames)
- `fit_summary_vS1.pdf`, `basis_sweep/`, `hsgp3D_*`, `hsgp_bathycov_*` — diagnostics/recovery

The sweep and demos accept env-var overrides for quick runs, e.g.
`SWEEP_CONFIGS`, `SWEEP_SEEDS`, `DEMO_MX/MY/MZ`, `DEMO_CHAINS/WARMUP/SAMPLE`.

---

## Status & caveats

- **vS1** is validated end to end: healthy sampling (near-zero divergences,
  Rhat ≈ 1.0), latent-field recovery R² ≈ 0.98/0.86/0.68 (hake/dolphin/humpback),
  and GP hyperparameters recovered at an adequate basis.
- **Rare species** (humpback) recover more weakly everywhere — sparse detections.
- **GP intercepts** (`mu_sp`) are partly confounded with the realised zero-mean
  field's sample mean; predictions stay accurate regardless.
- **vS2–4** ship the simulation stage only; the 3-D / multi-depth estimation
  pipeline is demonstrated (demos A/B) but not productionised.

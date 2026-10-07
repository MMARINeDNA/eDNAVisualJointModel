# Restructure roadmap: one eDNA, one visual, one joint HSGP model

Status: proposal, 2026-10-06. Nothing has been deleted yet.

## 1. Where things stand

### HSGP4eDNA (the debugged eDNA machinery)

- `stan/hsgp_nd.stan` is the qPCR hurdle + ZI-Beta-Binomial metabarcoding model
  with an n-D anisotropic HSGP (`phi_nD`, `spd_nD`, `INDICES`). The same file
  serves `gp2d` and `gp3d`.
- `stan/hsgp_2d_bathysp.stan` is the same model with a 2-D HSGP plus a
  bottom-depth spline (`gp2d_bathysp`).
- `R/functions.R` has the simulators (`simulate_nobathy`, `simulate_bathysp`),
  the Stan-data formatters, `hsgp_basis_rule()` and `bathy_spline_basis()`.
- Scripts are organized as scenario × model (`01_sim_<scenario>`,
  `0x_fit_<model>`). Every matched fit is validated except the accepted L1 gap.
- `README.md` has uncommitted edits (the "Known limitations" section). Commit
  them before pinning a version (Phase 0).

### eDNAVisualJointModel (this repo)

| Area | Contents | Status |
|---|---|---|
| `stan/whale_edna_hsgp_v1…v4.1.stan` | Older eDNA-only models (3-D loop-based HSGP) | **Obsolete.** Replaced by `hsgp_nd.stan`. |
| `scripts/older_simulations/` (35 files) | v1–v4.1 pipeline steps | **Obsolete.** |
| `00_pipeline_v*.r` (7 files) | Orchestrators | **Broken**: they source `scripts/01_…`, but those files now live in `scripts/older_simulations/`. |
| `outputs/whale_edna_output_v*/` | Tracked sim PDFs and notebook `.qmd` files | Obsolete. |
| `notebooks/v3*.html` | Write-ups of v3 and the v3.2 debugging work | Historical. |
| `distance/distance_hn_dens_v4.1.stan` + `00_distance_v4.1.R` | Non-spatial line-transect (LT) model | Superseded by v4.1a. Its detection logic is still the reference. |
| `distance/distance_v4.1a.stan` + `00_distance_v4.1a.R` | **Spatial visual model**: LT half-normal detection, group-size covariate, size-bias correction, 3-D HSGP | **Current visual model**, but it uses the old HSGP code (`hsgp_phi` triple loop, M = 3584, x/y/Z_bathy), not the debugged machinery. |
| `distance/distance_hn_dens.stan`, `distance_stan_density.Rmd` | Early trials | Byte-identical duplicates of the files in `Distance sampling trials/`. |
| `Distance sampling trials/` | Early half-normal trials (with PDFs) | Obsolete. |
| `distance/*_grpsz.RData` | Empirical group-size pools | **Gitignored** (the `distance/*` rule), so a fresh clone cannot run the distance pipeline. They can be rebuilt with `scripts/getGrpSz.R`. |
| `scripts/simulation/README.md` | Pointer to HSGP4eDNA | Stale: it still uses `vS1` / `whale_edna_hsgp_vS1.stan` names. |
| `scripts/plot*Data.R`, `getGrpSz.R`, `FinWhales.R`, `DetectionsBySpecies.R` | Real-data figures and exploration | Keep, but give them their own folder. |
| `presentations/`, `figures/`, `data/` | Collaborator decks, real-data maps, real data | Keep. |
| `README.md` | ~500 lines describing v1–v4.1 | Mostly stale. |

**There is no combined eDNA + visual model in either repo yet.** Despite the
repo name, the `whale_edna_hsgp_*` models are eDNA-only (qPCR + MB). The
"joint" in their history meant qPCR + metabarcoding. The combined model is new
work, not a cleanup.

## 2. Target end state

The guiding rule is **one copy of each piece of machinery**. HSGP4eDNA remains
the canonical home of the eDNA model and the HSGP helpers. This repo consumes
them; it does not copy them.

```
eDNAVisualJointModel/
├── README.md                 short: the three models, how to run them, results table
├── external/HSGP4eDNA/       git submodule, pinned commit (eDNA model + HSGP helpers)
├── stan/
│   ├── include/              visual-only Stan functions (detection fn, ESW, size bias)
│   ├── hsgp_visual.stan      HSGP + line-transect observation model
│   └── hsgp_joint.stan       shared latent field → eDNA obs + LT obs
├── R/
│   ├── functions_visual.R    LT simulator (transects, segments, sightings) + formatter
│   ├── functions_joint.R     joint formatter, species map eDNA↔LT
│   ├── 01_sim_<scenario>.r   one field → eDNA obs + LT obs (reuses HSGP4eDNA simulators)
│   ├── 0x_fit_edna.r         thin driver around HSGP4eDNA's model (comparison baseline)
│   ├── 0x_fit_visual.r
│   ├── 0x_fit_joint.r
│   └── 0x_compare.r          eDNA-only vs visual-only vs joint recovery on the same sim
├── tests/smoke.R             fast compile+sample+recovery check of all three models
├── data/                     real data (+ data/grpsz/*.rds tracked)
├── scripts/data_figures/     plotLTData.R, plot*Data.R, getGrpSz.R, exploration
├── presentations/
└── figures/
```

The **eDNA model** lives only in HSGP4eDNA. This repo's `fit_edna` driver
calls `external/HSGP4eDNA/stan/hsgp_nd.stan` directly. Shared Stan functions
(`phi_nD`, `spd_nD`, `lambda_nD`, `zi_beta_binomial_lpmf`, and ideally the eDNA
likelihood as an `_lp` function) are moved in HSGP4eDNA into
`stan/include/*.stan`. Then `hsgp_visual.stan` and `hsgp_joint.stan` pull them
in with `#include`, using
`cmdstan_model(..., include_paths = c("stan/include", "external/HSGP4eDNA/stan/include"))`.

Shared latent field in the joint model (animals/km²):

```
log λ_s(x) = μ_s + f_s(X, Y) [+ h_s(Z_bathy)]                 # from HSGP4eDNA
eDNA:   log λ_edna = log λ_s + log zsample_eff + log conv_s + log vol   → qPCR / MB
visual: λ_groups_s = λ_s / E[group size_s];  seg_count ~ Poisson(2 L esw_pop λ_groups)
```

Hake has only eDNA data. Humpback and PWSD have both eDNA and LT data, which
needs an explicit species-index map. The **visual and joint models both use
`gp2d_bathysp`**: a 2-D HSGP over (X, Y) plus a spline on bottom depth, the
structure of HSGP4eDNA's `hsgp_2d_bathysp.stan`. Plain `gp2d` is kept as the
special case with no spline basis (`K_bathy = 0`), which is useful for smoke
tests and the `nobathy_surface` scenario. Decided 2026-10-06.

## 3. Guardrail: "works at every step"

Two levels of checks are defined in Phase 1 and run at the end of every phase.

1. **Smoke test** (`tests/smoke.R`, target < 10 min). For each available model
   it simulates a small scenario (about 60 stations, about 10 transects),
   uses a small basis (6×4), compiles, and samples 2 chains × 200/200. It
   asserts that compilation succeeds, there are 0 divergences, max Rhat < 1.1,
   and field R² for hake/humpback exceeds a loose floor. The test runs before
   every commit that touches `R/` or `stan/`.
2. **Validation run** at each milestone. This is the full basis
   (`hsgp_basis_rule()`) with production chains. The results row
   (N, runtime, divergences, Rhat, field R², hyperparameter coverage) is
   recorded in the README results table, as HSGP4eDNA does.

A phase is done when its smoke test passes in a fresh clone (`git clone
--recurse-submodules` + `Rscript tests/smoke.R`) and its validation row is
recorded. Every deletion happens on a branch whose smoke test passes, behind
the tag `pre-restructure` so nothing is lost.

## 4. Phases

### Phase 0: Freeze and connect (no behavior change)

1. In HSGP4eDNA, commit the pending README edits and tag `v1.0`.
2. Here, tag the current `main` as `pre-restructure`. All old versions stay
   recoverable from the tag.
3. Add HSGP4eDNA as a submodule at `external/HSGP4eDNA`, pinned to `v1.0`.
4. Fix the gitignore hole: move the group-size pools to `data/grpsz/*.rds`
   (tracked; they are tiny), update `getGrpSz.R` and both distance scripts.

✅ Check: `00_distance_v4.1a.R` with `ITER_WARMUP=200 ITER_SAMPLING=200` runs
in a fresh clone. HSGP4eDNA's `04_fit_gp2d.r` runs from inside the submodule.

**Done 2026-10-06, with these deviations:**
- Step 1 waited for another HSGP4eDNA session to land `7592a94` (a README-only
  change: it closes the `bathysp_surface` fit gap and documents the
  `bathygp_depthspref` limitation). `v1.0` is tagged on `7592a94` and the
  submodule is pinned there. The tag exists locally only; push it to
  HSGP4eDNA's origin.
- The `pre-restructure` tag exists locally on `0aaa860`. It is not pushed yet.
- Also removed a stray committed gitlink (`.claude/worktrees/jovial-thompson-a04c95`,
  with no `.gitmodules` entry) that would break `--recurse-submodules` clones,
  and gitignored `.claude/worktrees/`.
- The group-size pools are now at `data/grpsz/{humpback,pwsd}.rds`. They are
  identical to the old `.RData` files (683 humpback, 110 PWSD sightings).

**Phase 0 check results** (fresh `--recurse-submodules` clone):

| Check | Result |
|---|---|
| HSGP4eDNA `nobathy_surface` × `gp2d` (01 → 05), run in the submodule | ✅ 24 min, 0 divergences, 0 treedepth hits, Rhat ≤ 1.01, E-BFMI ≈ 0.85 |
| `distance/00_distance_v4.1.R` (non-spatial) | ✅ 0 divergences; every truth inside its 95% CI, both species. Sporadic `ibeta_derivative` overflow messages only. |
| `distance/00_distance_v4.1a.R` (spatial), 200/200 | ⚠️ Runs end to end, but humpback has 50% divergences and PWSD detection σ collapses to its 0.001 floor |
| `distance/00_distance_v4.1a.R` (spatial), 1000/1000 | ❌ Humpback chain 3 stuck: step size ≈ 1e-5 (vs ≈ 6e-3 on the other chains), every draw at treedepth 12, lp ≈ 6×10⁶ away from the other chains. Stopped after 3.5 h at sampling draw 794 / 1000. |

The simulated truths in the v4.1 run differ from the tracked
`outputs/distance_v4.1/distance_v4.1_recovery.csv` (April 2026) even though the
seed is the same. The cause is not the group-size change, because the
`as.integer()` pools are identical. It is most likely RNG changes across R or
package versions.

**Consequence:** `distance_v4.1a` is **not a working baseline**. The spatial
visual model is unvalidated. The detection-function machinery in
`distance_hn_dens_v4.1.stan` is validated. Phases 1 and 3 below are adjusted
accordingly.

### Phase 1: Delete obsolete eDNA code and add the smoke harness

1. Delete `stan/whale_edna_hsgp_v*.stan`, `scripts/older_simulations/`,
   `00_pipeline_v*.r`, `outputs/whale_edna_output_v*/`,
   `Distance sampling trials/`, and the duplicate
   `distance/distance_hn_dens.stan` / `distance_stan_density.Rmd`.
2. Move the real-data scripts to `scripts/data_figures/`.
3. Write `tests/smoke.R` with the eDNA case: a thin `R/fit_edna.r` that
   calls HSGP4eDNA's simulator, formatter, and `hsgp_nd.stan` via the
   submodule. Also include the existing non-spatial `distance_v4.1` case
   unchanged, as the visual baseline that later phases must preserve. Do not
   include `distance_v4.1a`, which does not sample reliably (see the Phase 0
   results).
4. Rewrite `README.md` down to the target state, marking the parts that are
   not built yet.

✅ Check: the smoke test passes for eDNA and the non-spatial v4.1 visual model. Deleted
files are recoverable from `pre-restructure`.

**Done 2026-10-06.**

- Deleted 72 files: the v1–v4.1 eDNA Stan models, `scripts/older_simulations/`,
  `00_pipeline_v*.r`, `outputs/whale_edna_output_v*/`,
  `Distance sampling trials/`, the duplicate distance trial files,
  `notebooks/v3_notebook.html` + `notebooks/README.md`, and
  `scripts/simulation/README.md`.
- Moved the 9 real-data scripts (including the harbor porpoise scripts from
  #61) to `scripts/data_figures/`. Moved the v3.2 and distance v4.1 notebooks
  to `docs/history/`.
- The eDNA driver is `R/fit_edna.r`. Per the gp2d_bathysp decision, it drives
  `hsgp_2d_bathysp.stan` via HSGP4eDNA's own `09_fit_gp2d_bathysp.r`, not
  `hsgp_nd.stan`.
- `distance/00_distance_v4.1.R` gained the env overrides `OUTPUT_DIR`,
  `CHAINS`, `ITER_WARMUP` and `ITER_SAMPLING`. Its defaults are unchanged. It
  now also saves its divergence diagnostics and max Rhat.
- `.gitignore` drops the obsolete rules. Its `stan/` rule now whitelists
  `.stan` files at any depth, so `stan/include/` will be tracked.
- The presentation presenter notes still mention old paths. That is left for
  Phase 6, since those notes need re-rendering anyway.

Smoke test (`Rscript tests/smoke.R`): **13/13 checks pass, 5.6 min**.

| Case | Setup | Result |
|---|---|---|
| edna | `bathysp_smoke`: 100 stations, basis 8×4, spline df 4, 2 chains × 200/200 | 0 divergences, max Rhat 1.02, field R² hake 0.96 / humpback 0.59 / PWSD 0.66 (hake checked > 0.5) |
| visual | v4.1, 2 chains × 300/300 | 0 divergences, max Rhat ≤ 1.016; σ and D truths inside their 95% CIs for both species |

### Phase 2: Refactor shared Stan functions (in HSGP4eDNA)

1. **(Required.)** Move the shared functions out of `hsgp_2d_bathysp.stan`, the
   reference eDNA model for this repo, and of `hsgp_nd.stan`. They go into
   `stan/include/hsgp_functions.stan` (`phi_nD`, `spd_nD`, `lambda_nD`) and
   `stan/include/edna_functions.stan` (`zi_beta_binomial_lpmf`). Also move the
   qPCR + MB likelihood into an `edna_obs_lp(...)` function, so the joint
   model calls exactly the same code instead of copying about 300 lines. The
   spline design matrix is already built in R (`bathy_spline_basis()` in
   `R/functions.R`), so the visual and joint models call that directly.
2. *(Optional tidy-up, not needed by this repo.)* `hsgp_2d_bathysp.stan` is
   `hsgp_nd.stan` plus about 10 lines of spline term. It is already
   dimension-general, so relaxing `K_bathy` to `<lower=0>` would let it
   replace `hsgp_nd.stan` (`K_bathy = 0` gives `gp2d` / `gp3d`). That is
   HSGP4eDNA's call, since it keeps 3-D and spline fits side by side.
3. Turn `log_lik_qpcr` / `log_lik_mb` back on in generated quantities (L4 needs
   them anyway, for model comparison).

✅ Check: HSGP4eDNA's matched fits reproduce with identical posteriors for a
fixed seed, or within MC error. Re-run at least the `nobathy_surface × gp2d`
and `bathysp_surface × gp2d_bathysp` cells. Tag `v1.1` and bump the submodule.

**Done 2026-10-07** ([HSGP4eDNA#1](https://github.com/MMARINeDNA/HSGP4eDNA/pull/1),
tag `v1.1` = `a6aa2b5`; submodule bumped).

- `stan/include/hsgp_functions.stan` provides `lambda_nD`, `spd_nD`,
  `phi_nD`, `hsgp_phi()` and `hsgp_sqrt_spd()`.
- `stan/include/edna_functions.stan` provides `zi_beta_binomial_lpmf`, the
  per-observation `edna_qpcr_loglik()` / `edna_mb_loglik()`, and the MB
  shape helpers.
- Step 1's `edna_obs_lp` became two per-observation log-likelihood
  functions. The model block uses `target += sum(...)` and generated
  quantities uses the same calls for `log_lik_qpcr` / `log_lik_mb`, which
  also covers step 3.
- Step 2 (merging the two models) was not done; it is optional.
- Consumers compile with
  `include_paths = c(..., "external/HSGP4eDNA/stan/include")` and use
  `#include hsgp_functions.stan` / `#include edna_functions.stan` inside
  `functions { }`.

| Check | Result |
|---|---|
| `log_prob` + `grad_log_prob`, old vs new, 20 random unconstrained points | ≤ 1.3e-14 relative for `gp2d`, `gp3d` (D1 = 3), `gp2d_bathysp` |
| `nobathy_surface` × `gp2d` re-run | 0 divergences, 0 treedepth hits; posterior means match v1.0 within MC error. Four params at Rhat 1.01–1.02 (v1.0: none above 1.01), which is chain variation since the log density is identical |
| `bathysp_surface` × `gp2d_bathysp` re-run | 0 divergences, max Rhat 1.01; field R² 0.974 / 0.563 / 0.739 (v1.0: 0.98 / 0.56 / 0.74) |
| `log_lik` | All finite; `loo` runs (99% Pareto k good on a 60-station fit) |

### Phase 3: One visual model on the debugged machinery

1. Move the LT simulator out of `00_distance_v4.1a.R` into
   `R/functions_visual.R`. It should take its latent field from HSGP4eDNA's
   simulators (same domain, same `default_gp_params()`, same coordinate
   normalization), so the eDNA and LT data can share one truth.
2. Write `stan/hsgp_visual.stan`. Copy the detection function,
   population-average ESW, group-size model, and size-bias correction
   verbatim from the **validated** `distance_hn_dens_v4.1.stan` into
   `stan/include/visual_functions.stan`. Use the spatial wiring in
   `distance_v4.1a.stan` (per-segment `lambda_groups`, segment PPCs,
   `log_lik`) as a template only, because v4.1a never sampled reliably.
   Replace its 3-D triple-loop HSGP with the `#include`d `gp2d_bathysp`
   latent field: a 2-D HSGP over (X, Y) from `phi_nD` / `spd_nD`, plus the
   bottom-depth spline. Take the basis size from `hsgp_basis_rule()`.

   **Why v4.1a failed (diagnosed 2026-10-06).** `hsgp_weights` evaluates the
   SE spectral density with length-scales and frequencies in km and m. The
   eigenfunctions (`phi1d`) are orthonormal on the *normalized* domain
   [−L, L]. The per-dimension Jacobian 1/`coord_scale` is missing, so every
   basis weight is √(250·635·1750) ≈ 16,700× too large. In the 200/200 fits,
   essentially 100% of segments in every chain sat at the log λ clamps (+15
   or −10). The clamps have zero gradient, so chains stalled (step size down
   to ~5e-7), and detection σ collapsed to absorb the impossible encounter
   rates. The GP length-scales only *looked* recovered: their priors are
   tight and centred on the truth, and the saturated field passed no
   information back to them. HSGP4eDNA avoids this by sampling `gp_l_raw` in
   normalized units and converting to km only in generated quantities.
   Porting onto that code fixes the bug by construction.

   Two lessons for the new model:
   - Drop the `fmin`/`fmax` clamps on log λ and log σ. They hide this class
     of bug and create zero-gradient regions.
   - Use the same weakly-informative `gp_l` priors as HSGP4eDNA instead of
     tight priors centred on the truth, so recovery actually tests something.
3. Validate in two steps:
   - (a) Fit the **non-spatial special case** (`M` tiny / `gp_sigma` → 0) to
     the same data as `distance_hn_dens_v4.1.stan`. It should match the old
     ESW, σ₀, β_size, and density estimates.
   - (b) Run full spatial recovery for humpback and PWSD.
4. Once (a) and (b) pass, delete `distance/` entirely (both old models,
   scripts, and `outputs/distance_v4.1*`). Keep the distance notebook (see
   the decision on historical write-ups).

✅ Check: the smoke test has the new visual case and the old case is removed.
The validation row is recorded.

**Done 2026-10-07.**

New files:

| File | Contents |
|---|---|
| `stan/include/visual_functions.stan` | `hn_sigma`, `hn_esw`, `group_size_pmf`, `esw_population`, and the per-observation `lt_detection_loglik` / `lt_segment_loglik` |
| `stan/hsgp_visual.stan` | Single species. Field `log λ = μ + f(X, Y) + B(Z_bathy)·β`, in animals/km² like the eDNA model. Expected detected groups = λ / E[s] · 2L · ESW_pop, where E[s] is the mean of the modelled group-size pmf. Data switches `use_gp` and `K_bathy`. Emits `log_lik_det`, `log_lik_seg` and `pp_seg_count`. |
| `R/functions_visual.R` | `simulate_field_bathysp()`, the field at any locations. It uses HSGP4eDNA's `default_gp_params()`, `bathy_spline_basis()` and the same true spline coefficients. Also `lt_design()`, `lt_species_params()` (v4.1 detection and priors), `simulate_lt_sightings()` (v4.1 logic), `simulate_visual()` and `format_stan_data_visual()`. Coordinates are normalized by the domain extents, so the joint model can share them. |
| `R/fit_visual.r` | `fit_visual()`, recovery summaries, and script mode |
| `R/validate_visual.r` | Multi-replicate recovery study |

Removed: `distance/` and `outputs/distance_v4.1/`. The distance notebook HTML
stays in `docs/history/`.

Deviations from the plan:

- **σ clamp kept.** The `log σ ≤ 30` clamp is kept in `hn_sigma`. It
  engages only for σ > 1e13 km, where ESW is flat at w, and it prevents
  `exp()` overflow in the S_max integration loop. The harmful v4.1a clamps
  were on log λ; there are none here (`poisson_log_lpmf` is used).
- **`mu_sp` prior.** It is `normal(-5, 2)` on log animals/km², weakly
  informative and not centred on any truth. The other field priors are
  HSGP4eDNA's (see the open issue below).
- **Bottom depth in the simulation.** Z_bathy at segments is drawn
  independently of (X, Y), as in HSGP4eDNA's `simulate_bathysp()`, so the
  spline is identifiable. Real bathymetry is spatially smooth, which will
  partly confound the spline with the GP. That is a Phase 6 concern.
- **Duplicated truth coefficients.** `BATHY_BETA_TRUE` duplicates the true
  spline coefficients hard-coded inside HSGP4eDNA's `simulate_bathysp()`.
  Phase 4 needs one field shared with the eDNA observations, so it should
  factor the field generation out of `simulate_bathysp()` upstream.

**Check (a): non-spatial reduction vs v4.1, on the same data** (3 seeds × 2
species, 4 × 1000/1000; the new model's `mu_sp` prior is matched to v4.1's
group-density prior):

| Quantity | Max \|relative difference\| in posterior median | CI width ratio new / v4.1 |
|---|---|---|
| σ at s_centre | 0.19% | 0.96–1.04 |
| ESW_pop | 0.9% (PWSD consistently ~0.8% lower) | 0.95–1.04 |
| group density | 1.1% | 0.94–1.05 |
| β_size (PWSD) | 3.0% | 0.98–1.00 |

There were 2 divergences for v4.1 and 1 for the new model, all on PWSD seed
2. The small PWSD ESW shift comes from putting the density prior on animals
rather than groups, which couples it weakly to the group-size parameters.

**Check (b): spatial recovery** (`SEEDS=1:5 Rscript R/validate_visual.r`; basis
14×6 = 84 from `hsgp_basis_rule()`, spline df 4; 4 × 1000/1000; 28 min total):

| Species | Detections | Divergences | max Rhat | min E-BFMI | Field R² | Coverage D / σ / E[s] / μ | Median rel. bias D |
|---|---|---|---|---|---|---|---|
| Humpback | 549–917 | 0,1,0,0,0 | ≤ 1.007 | 0.70 | 0.84–0.93 | 5/5, 5/5, 5/5, 5/5 | −3% |
| PWSD | 41–100 | 1,0,0,0,0 | ≤ 1.003 | 0.79 | 0.59–0.75 | 5/5, 4/5, 5/5, 4/5 | −4% |

Smoke test: the visual case is now `hsgp_visual` (both species, basis 8×4, 2 ×
300/300), replacing distance v4.1. **12/12 checks pass, 6.9 min.** It also
caught a name-capture bug: `summarise_visual()` used a bare `rhat`, which a
caller's variable could shadow. The fix namespaces it as `posterior::rhat`.

**Open issue: GP hyperparameter priors.** HSGP4eDNA's field priors pull the
GP hyperparameters for both whales:

| Parameter | Prior | Truth | Result |
|---|---|---|---|
| `gp_sigma` | `gamma(8, 4)`, mean 2 | 1.0 / 1.3 | Median ~40% high; coverage 3/5 (humpback), 4/5 (PWSD) |
| `lx` | `gp_l_raw ~ gamma(10, 16)`, mean ≈ 156 km | 50 km | +17% (humpback, 3/5), +99% (PWSD, 2/5) |
| `ly` | same | 300 km | +21–24%; coverage 4/5, 5/5 |

These priors were tuned on hake's dense eDNA data. In v3.2 the tight
`gp_sigma` prior was introduced to escape a low-σ trap. Density and the field
are still recovered, but the joint model should settle the priors first,
since eDNA and visual data will share them. Options:
- keep them for consistency;
- loosen them for the whale species;
- make them species-specific data.

### Phase 4: Joint simulator

`R/01_sim_<scenario>.r` draws **one** latent field per species and generates:

- eDNA samples (HSGP4eDNA's design: stations, qPCR, MB with junk), and
- LT effort and sightings (transects, half-normal detection, group sizes) for
  the cetaceans.

Start with `bathysp_surface`, the matched truth for `gp2d_bathysp`. Use
`nobathy_surface` for the smoke test. This simulator is the only one in the
repo; the eDNA-only and visual-only fits use subsets of its output.

✅ Check: the eDNA subset fit with HSGP4eDNA's model and the LT subset fit with
`hsgp_visual.stan` each still recover the field. This proves the shared sim is
consistent with both single-source models before any joint fitting.

### Phase 5: Joint model

1. Write `stan/hsgp_joint.stan`: one `gp2d_bathysp` field per species, then
   `edna_obs_lp(...)` (from the include) for all S species plus the visual
   likelihood (from the include) for the LT species. It needs a species map
   (e.g. `lt_species_idx = {2, 3}`).
2. Hold `conv_factor` fixed as data first, matching HSGP4eDNA.
3. Run `0x_compare.r` on one simulated dataset:
   - eDNA-only vs visual-only vs joint.
   - Field R² and CI width should satisfy joint ≥ max(single-source) for
     humpback and PWSD.
   - Hake should be unchanged from eDNA-only.
4. *(Science step.)* Free `log_conv_factor` for the cetaceans. The LT data
   pin down absolute density, so the joint model should identify the eDNA
   shedding/conversion rate that eDNA alone cannot. Report its recovery.

✅ Check: the smoke test now covers all three models, the validation table has
eDNA / visual / joint rows, and the comparison figure is produced.

### Phase 6: Extend and tidy

- Add the remaining scenarios to the joint model: `bathysp_depthspref` (water
  column + depth preference). Also add a `bathygp_*` truth as a deliberate
  mis-specification check.
- Add the real-data entry point (`01_load_muri.r`) that produces the same
  object shape as the simulator.
- Write a final README results table. Re-render or retire presentations that
  reference v-numbers.

## 5. Decisions for you

1. **How this repo depends on HSGP4eDNA.** I recommend a git submodule (pinned,
   reproducible, and Stan `#include` works across it). The alternatives are a
   sibling-checkout path set by an env var (simpler, but unpinned, so it breaks
   when HSGP4eDNA changes) or vendoring copies (simplest, but it recreates the
   duplication we are removing).
2. **Historical write-ups** (`notebooks/v3*.html`, `distance_v4.1_notebook`,
   the tracked sim PDFs). I recommend keeping the v3.2 debugging case study
   and the distance notebook under `docs/history/` because they document
   decisions, and deleting the rest (the tag preserves them).
3. **Latent-field structure.** ~~`gp2d` first~~. Decided 2026-10-06: the
   visual and joint models use `gp2d_bathysp` (2-D HSGP + bottom-depth
   spline), with `gp2d` as the `K_bathy = 0` special case. Order is unchanged:
   Phase 1 → Phase 2 (makes the HSGP and eDNA functions includable) → Phase 3 (visual).
4. **Phase 2 changes HSGP4eDNA itself** (includes, the optional model merge,
   `log_lik`). That repo has its own validation table and collaborators. If
   you'd rather leave it untouched, the only fallback is to vendor its Stan
   functions into this repo's `stan/include/`, plus a smoke-test check that
   they still match upstream. A monolithic `.stan` file can't be `#include`d
   piecemeal.

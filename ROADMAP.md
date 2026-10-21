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

**The demonstration goal: three tiers** (added 2026-10-20). The model is to be
shown to work on, in order:

1. **An idealised simulation**, *in progress* (Phases 4–6). The truth is
   drawn from the model's own structure (a GP plus a bottom-depth spline), on
   uniform eDNA stations and evenly spaced east–west transects across the
   study domain. This is the joint comparison in Phase 5 and 6.0, now on
   the 35°N domain.
2. **A realistic simulation**, *planned* (6.4). The truth is previously
   modelled density surfaces for PWSD, humpback and hake. Visual and eDNA
   observations are simulated from them at the **actual** ship transects and
   eDNA station locations. The test is how well the model recovers the "true"
   surfaces under the real sampling design, when the truth is not a draw from
   the model.
3. **Real data**, *paused* (6.1, 6.3) while additional data are sourced.

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

**Resolved 2026-10-07: design-based GP priors** (`R/gp_priors.R`).

*Principle (decided 2026-10-07):* priors depend only on the survey **design**
(sample locations, domain extent) and the compute budget, never on observed
outcomes. A version that also used detection counts, through an "information
floor", was dropped for this reason. Species with the same design therefore
get identical priors.

*Length-scales.* Inverse-gamma on each axis, with 1% of mass below ℓ_min and
1% above ℓ_max (Betancourt's boundary-avoiding construction). ℓ_max is the
domain extent along the axis. ℓ_min is the larger of two floors:
- the design floor, the median gap between distinct sample coordinates;
- the compute floor, the shortest scale a basis of `GP_BASIS_BUDGET` = 400
  functions can represent.

When the design floor would need a bigger basis, both axes are scaled by a
common factor. The basis is sized to ℓ_min, so it covers the whole prior.

*Marginal SD.* `gp_sigma ~ gamma(4.05, 3.37)`, with 1% of mass below 0.25 and
1% above 3.

*Why a budget of 400.* It is a compute budget, not a statistical quantity.
Uncapped, the design floor asks for 484–2178 basis functions. With 400, the
floor stays below the 50 km simulation length-scale for every simulated data
type (lx floor 23–33 km). A cap of 128, or HSGP4eDNA's "M well under the
number of observations" rule (M = 100 for 200 eDNA stations), would put the
floor at 41–66 km and exclude the truth. That rule came from overfitting
under a gamma prior that allowed ℓ → 0, which the inverse-gamma floor
prevents.

For the transect design the rule gives ℓ_min = (23.5, 120) km and a basis of
28×14 = 392 for both species.

*Results* (`SEEDS=1:5 Rscript R/validate_visual.r`; the old-prior run is
archived as `outputs/visual/validation_legacy_priors/`):

| | Old (HSGP4eDNA) priors, M = 84 | Design priors, M = 392 |
|---|---|---|
| Humpback gp_sigma / lx / ly coverage | 3/5, 3/5, 4/5 | **5/5, 5/5, 5/5** |
| Humpback gp_sigma / lx median bias | +40% / +17% | **+3% / −6%** |
| PWSD gp_sigma / lx / ly coverage | 4/5, 2/5, 5/5 | **5/5, 4/5, 5/5** |
| PWSD gp_sigma / lx / ly median bias | +38% / +99% / +24% | **+5% / +48% / −3%** |
| Field R² humpback / PWSD | 0.84–0.93 / 0.59–0.75 | 0.87–0.93 / 0.60–0.79 |
| Density coverage, bias | 10/10, −3 / −4% | 10/10, −3 / −3% |
| Divergences | 1 each in 2 fits | 1 each in 2 fits |
| Runtime per fit (humpback / PWSD) | 1.8–2.6 / 3.0–3.6 min | 2.1–4.8 / 5.2–7.9 min |

*Floor check* (`floor_check()`, a truth-free basis-adequacy test). It flags a
fit when P(ℓ < 1.25 ℓ_min) in the posterior exceeds max(3 × the prior's
probability, 5%). It flagged 1 of 10 fits: PWSD seed 3, lx, with 11%
posterior vs 3.6% prior. Refitting at M_max = 800 moved lx from 39 [25, 67]
to 36 [22, 61] km and left field R² at 0.60. The flag cleared and runtime
went from 6.7 to 11.6 min. So it was a false alarm for the field, but it
shows that lx is mildly prior-sensitive with sparse data.

Procedure for real data: on a flag, refit at 2× the budget, and keep the
larger basis only if the field or density changes.

*For the joint model.* Build the prior from the union of eDNA station and LT
segment locations. Port `inv_gamma` length-scale priors into
`hsgp_joint.stan`.

*Upstream ([HSGP4eDNA#3](https://github.com/MMARINeDNA/HSGP4eDNA/pull/3),
tag `v1.3`; submodule bumped).* The design-based prior is now HSGP4eDNA's
default, and `R/gp_priors.R` is gone (both models use the upstream
`gp_prior_from_design()` / `floor_check()`). v1.3 also adds reproducible fit
seeds and two noise terms:
- **Between-bottle noise:** exact negative binomial, one θ shared by all
  species.
- **Junk noise:** a per-sample junk term with its own θ_junk.

With both, θ covers its truth of 10 in every matched cell. `hsgp_joint.stan`
must carry both terms.

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

**Done 2026-10-07.**

- *Upstream ([HSGP4eDNA#2](https://github.com/MMARINeDNA/HSGP4eDNA/pull/2), tag `v1.2`; submodule bumped).*
  `simulate_bathysp()` is split into
  `simulate_field_bathysp(locs, gp_params, df)`, which draws the field at any
  locations, plus `simulate_edna_obs()` and `zsample_effect_matrix()`.
  `BATHY_BETA_TRUE` and `aniso_cov_2d()` are now top level. Output is
  **bit-identical** (`identical()`) for `bathysp_surface` (seeds 202, and 7
  with 120 stations) and `bathysp_depthspref` (seed 203).
- *Visual simulator.* It now calls the upstream field function through a thin
  `species_field()` wrapper, so the duplicated `BATHY_BETA_TRUE` is gone.
  `simulate_visual()` truth changed by ≤ 8.9e-16 (matrix × matrix vs per-species
  matrix × vector products). Sightings are identical, so the Phase 3 and
  design-prior validations stand.
- *`R/functions_joint.R`.* `simulate_joint()` builds the eDNA stations (as in
  `simulate_bathysp()`) and the LT segments. It draws **one** field over all
  their locations, with the spline basis from all locations, then simulates
  both observation types. Its `$edna` half has the shape of
  `simulate_bathysp()` output and its `$visual` half the shape of
  `simulate_visual()` output, so every single-source tool works unchanged.
  `joint_field_agreement()` checks that the two halves share the field.
- *`R/fit_edna.r`* accepts a pre-built `sim`.
- *`R/validate_joint_sim.r`* is the check below.
- *`tests/smoke.R`* gains a `joint_sim` case (no fits): field sharing and
  both formatters. 16/16 checks pass.

Check (`Rscript R/validate_joint_sim.r`, 2 seeds; 55 min):

| | Seed 1 | Seed 2 | Single-source reference |
|---|---|---|---|
| GP field cor, station/segment pairs < 10 km (all species) | 0.995–0.996 | 0.995–0.997 | — |
| eDNA (HSGP4eDNA 09, 17×6, df 4, 3 × 400/400) field R² hake / humpback / PWSD | 0.98 / 0.73 / 0.68 | 0.98 / 0.58 / 0.65 | 0.97 / 0.56 / 0.74 |
| eDNA divergences, max Rhat | 0, 1.009 | 0, 1.027 | 0, 1.01 |
| Visual humpback R² (detections) | 0.90 (688) | 0.87 (422) | 0.87–0.93 (549–917) |
| Visual PWSD R² (detections) | 0.54 (37) | 0.58 (43) | 0.60–0.79 (41–100) |
| Visual divergences, floor flags, D covered, GP hyperparameters covered | 0, 0, 2/2, 6/6 | 0, 0, 2/2, 6/6 (σ_det missed for PWSD) | — |

The low PWSD visual R² tracks its few detections in these two draws. **What
this sets up for Phase 5:** PWSD is weak in both sources (eDNA R² ≈ 0.65–0.68,
visual ≈ 0.54–0.58), so it is where joint fitting should help most. Humpback
is strong visually but only moderate in eDNA.

Notes for Phase 5:
- The eDNA formatter normalizes coordinates by the *data* range; the visual
  formatter uses the *domain* extents. The joint formatter must use one
  normalization (the domain) for stations and segments.
- The design-based prior should be built from the union of station and
  segment locations.

### Phase 4b: Sea-state detection covariate (visual model)

The visual model's detection function has one half-normal scale per species,
plus the PWSD group-size term, so ESW is the same on every segment. On real
surveys Beaufort sea state is usually the strongest detection covariate,
especially for small cetaceans like PWSD. If sea state is uneven across the
area, one ESW reads poor conditions as low density, and that bias lands in the
spatial field. This has to be in place before Phase 5 runs on real data.
(Starts after the HSGP4eDNA junk-noise validation finishes.)

1. **Model.** `log σ = log σ₀ + β_size·(s − s̄) + β_bft·(bft − b̄)`, with the
   group-size term only for PWSD. Use segment-level Beaufort (`avgbft`) for
   both the segment's ESW and its sightings' distances. Sighting-level Beaufort
   is the alternative, and it could replace the segment value for the distance
   term later. ESW becomes per segment: for PWSD, the size-averaged ESW is
   computed once per binned Beaufort level (e.g. 0.5 steps), not per segment.
   That keeps the cost at a handful of size-pmf loops.
2. **Switch and prior.** Add a `use_bft_covar` switch (0 = current model,
   exactly). β_bft gets a weakly informative prior, negative-leaning: detection
   falls as sea state rises.
3. **Simulator.** `simulate_visual()` draws segment Beaufort with a spatial
   trend (e.g. rougher offshore and to the north), so leaving the covariate out
   visibly biases the field. Detection uses the true β_bft.
4. **Checks.**
   - With the switch off, the log-density is identical to the current model.
   - On the new simulator, the fit with the covariate recovers β_bft, σ₀, D
     and the field.
   - The fit without the covariate shows the expected bias, recorded as a
     validation row.
   - The smoke test gains a sea-state case.
5. **Joint model.** Phase 5 carries the same detection function. The real-data
   loader keeps `avgbft`, `avgvis` and swell on segments, and sighting-level
   conditions on the observations.

**Done.**
- **Model and formatter.** `use_bft_covar` and the binned levels are in
  `stan/hsgp_visual.stan`. The likelihood functions in `visual_functions.stan`
  now take a per-detection log σ and per-detection / per-segment ESW.
  `format_stan_data_visual()` switches the covariate on when the segments have
  a `bft` column (0.5-step bins, centre Beaufort 3). β_bft ~ normal(−0.1, 0.2).
- **Exactness.** With the switch off, log density and gradient match the
  single-ESW model to 0 (both species, 3 random points). With it on and
  β_bft = 0, the difference is exactly β_bft's prior term (−0.125).
- **Simulator.** `simulate_visual(sea_state = TRUE)` draws Beaufort after the
  field, so the field matches the default draw for the same seed. The default
  draws are unchanged. Truths: β_bft = −0.08 (humpback), −0.15 (PWSD).
- **Validation** (`SEA_STATE=1 SEEDS=1:3`, README results).
  - β_bft is covered 3/3 for both species. It is weakly identified for PWSD.
  - Without the covariate, the field shifts −0.03 to −0.12 log units per
    Beaufort unit, in the expected direction.
  - Everything else keeps its coverage.
  - Cost: PWSD fits are about 2.5× slower.
- **Smoke test.** New `visual_bft` case: humpback with the covariate, checking
  that β_bft is covered.
- **Side finding.** Loading Stan model methods into an R session
  (`compile_model_methods = TRUE`) changes which field a seed simulates. The
  draw is still valid, but not the same one. Simulate before loading model
  methods in test scripts.

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

**Done.**

- **Model.** `stan/hsgp_joint.stan` gives each species one field (2-D HSGP +
  bottom-depth spline, animals/km²), shared by two observation models:
  - **eDNA:** HSGP4eDNA v1.3, unchanged: qPCR + ZI-beta-binomial MB, one shared
    θ for bottle noise, own θ_junk for junk noise.
  - **Visual:** for `lt_species_idx = {2, 3}`. Single ESW per species, with the
    detection key (half-normal / hazard-rate) and per-species truncation from
    main's visual model.

  Switches `use_edna` / `use_visual` give eDNA-only and visual-only fits from
  the same model, priors, basis and starting values. `n_conv_free` frees
  conversion factors (step 4).
- **Data.** `format_stan_data_joint()` normalizes stations and segments by the
  domain. It builds the design prior from the union of their locations
  (21×19 basis), and reuses both single-source formatters for the observation
  data.
- **Exactness.**
  - Joint with `use_visual = 0` equals `hsgp_2d_bathysp.stan` up to a constant.
  - Joint with `use_edna = 0` equals the sum of the two `hsgp_visual.stan` fits
    up to a constant, with half-normal and hazard-rate humpback.
  - Spreads ≤ 3e-11; re-checked after merging main's detection-interface
    changes.
- **Comparison** (`R/compare_joint.r`, seeds 1 and 2; README results table):
  - **Hake:** unchanged.
  - **Humpback:** joint beats the best single source (visual) in R²
    (0.92 vs 0.90; 0.90 vs 0.87) and CI width in both seeds.
  - **PWSD:** joint has the best R² in both seeds (0.66 vs 0.66; 0.81 vs
    0.75). CIs are narrower in seed 2 but 4% wider than eDNA-only in seed 1,
    where the visual data shorten PWSD's along-coast length-scale (435 → 199
    km; truth 300, covered).
  - So the "joint ≥ max(single source)" check holds for R² everywhere, and for
    CI width except PWSD in seed 1.
- **Step 4: free conversion factors.** Prior normal(log 100, 2). The joint
  model recovers humpback 5.33 [5.09, 5.57] (truth 5.30) and PWSD 4.94
  [4.46, 5.42] (truth 4.70). The field is unchanged.
- **Mixing caveat.**
  - The first seed-2 joint run did not converge: max Rhat 1.66, one chain in a
    large-junk-noise mode (θ_junk ≈ 2 vs ≈ 17), hake R² 0.84.
  - A rerun on the same data converged and is the reported run.
  - Rhat catches it. Candidate fixes if it recurs: a non-centred (lognormal)
    junk term, or more chains.
- **Cost.** Joint and eDNA-only fits take 2–3 h (3 × 800/800); visual-only
  12–15 min.
- **Smoke test.** New `joint` case; 26 checks pass on the merged code.
- **Not yet in the joint model.**
  - The sea-state covariate: `format_stan_data_joint()` warns and ignores
    `use_bft` (MURI's PWSD requests it).
  - A real-data entry point (Phase 6).
- **ESW note.** Results were fitted before main's binned group-size ESW was
  merged. The binned sum matches the exact one to ~1e-5, far below anything
  compared here.

### Phase 6: Real data and extensions

The real-data loader (`R/01_load_muri.r`, `R/02_check_muri.r`), the visual
model's hazard-rate key, detection covariates and per-species truncation, and
the explainer (`docs/explainer.md`) came in from a parallel session (#78, #80,
#81). The explainer's real-data fits set most of the agenda below.

#### 6.0 Study domain: north of 35°N (code done; re-runs pending)

Every analysis is restricted to the MURI study area north of 35°N.

- **One domain everywhere.** HSGP4eDNA's `STUDY_DOMAIN` (UTM 10N, south-west
  corner (−10, 3870) km, 690 × 1770 km, southern edge on 35°N) is used by:
  - every HSGP4eDNA simulator (`simulate_bathysp()`,
    `simulate_nobathy()`, the `bathygp_*` / `nobathy` scenario
    scripts);
  - this repo's visual and joint simulators (`LT_DOMAIN` points at it);
  - the real-data loader, which rebuilds the box from the data and stops if
    it differs.
- **Real data.**
  - All 551 eDNA bottles (37.9–48.4°N) are kept.
  - Effort is kept where the segment midpoint is at or north of 35°N: 1,641
    of 2,779 segments, 7,561 of 12,967 km. Sightings on dropped segments are
    removed and counted in `sight_flow`.
  - The domain shrinks from 1360 × 2600 km to 690 × 1770 km.
- **Simulations.** The box grows from 500 × 1270 km (about 37.7–49.2°N) to
  690 × 1770 km. There are still 200 uniform stations; the visual design is
  now 35 transects of 69 segments (2,415 segments, was 25 × 50 = 1,250).
- **Length-scales and basis budget (decided 2026-10-12).** The domain is
  1.9× larger in area. At the old 400-function budget, the design floor for
  200 stations would rise to about 43 × 130 km.
  - **Truths:** now 50 km east–west × 150 km north–south for every species.
    Was 50 × 300 km, with PWSD 40 km east–west in the `bathygp` scenarios.
    150 km is closer to the real visual fits (85–130 km north–south).
  - **Budget:** `GP_BASIS_BUDGET` 400 → 800.

  Floors at 800, with the margin below the truths:

  | Design | Floor (east–west × north–south) | Margin |
  |---|---|---|
  | eDNA | 30 × 93 km | about 1.6× |
  | joint | 30 × 94 km | about 1.6× |
  | visual | 23 × 117 km | 2.2×, but only **1.28×** north–south |

  The visual north–south floor is limited by the transect spacing, so the
  simulated visual design now has 35 east–west transects (about 51 km apart,
  2,415 segments; was 25). Runtime grows roughly in proportion to the budget:
  expect eDNA and joint fits of about 4–6 h. The 800 budget also applies to
  the real-data fits.
  - **Why the budget matters.** The budget is a compute choice, but it sets
    the prior's lower length-scale bound wherever the design is denser than
    the basis, which is every set-up so far. So it acts as the resolution
    limit. A squared-exponential length-scale ℓ means correlation 0.61 at
    distance ℓ and about 0.05 at 2.45 ℓ, so hot spots about 2ℓ across.
    Resolving 50 × 50 km everywhere would take about 900 functions;
    25 × 25 km about 3,400.

**Re-run list, in order.** These run on stats3 with
`tools/stats3_reruns.sh` (stages 1–5); `notes/stats3_guide.md` (and its PDF)
covers the set-up, from cloning and copying the private data to pulling the
results back.

**Priority for the October 2026 presentation** (decided 2026-10-13): working
eDNA-only, visual-only and joint fits on simulations of the real study
domain. So the joint comparison (item 6) runs first, on several seeds:
seeds 2 and 3, then 4 and 5 (`CMP_SEEDS`), then seed 1 with the data-based
starting field (item 5). The full validation (items 3, 5) and the explainer (item 7) follow
after the presentation. *(Done: seed 1 ran in time for the talk; the replicates
moved to the sea-state comparison, see item 6.)*

1. **Real data objects.** `Rscript R/01_load_muri.r` (needs
   `data/MV1_MURI_df.csv` and `data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv`, which
   stay private and git-ignored), then `R/02_check_muri.r`. *Done for the
   loader in #85*: the committed `data/processed/muri_{edna,visual}.rds` now
   hold the 35°N data (with humpback's group-size effect switched on). The
   stats3 run repeats it and checks that it reproduces them.
   - Per-species truncation (95th percentile) is recomputed on the kept
     sightings.
   - The detection-function choices (hazard-rate humpback; half-normal +
     Beaufort PWSD) came from AIC on all effort; re-check them on the kept
     set.
2. ✅ **Smoke test** (`tests/smoke.R`) and `tests/visual_detection.R`: 26/26
   checks on main (Mac, 2026-10-12) and on the v1.5 code. stats3's stage 1
   repeats them.
3. **HSGP4eDNA 2-D validation** (on v1.5):
   - `nobathy_surface` × gp2d and `bathysp_surface` / `bathysp_depthspref` /
     `bathygp_surface` × gp2d_bathysp, giving a new README §4 table;
   - the §5 noise-term and §6 junk tables if they are to be kept current;
   - release a v1.5.x with the updated tables.

   The §2 basis sweep and the 3-D cells are v1.x history; mark them as on
   the old box rather than re-run.
4. ✅ **Visual validation** (Mac, 2026-10-12). `SEEDS=1:5` and `SEA_STATE=1
   SEEDS=1:3`: every fit converged (0 divergences, max Rhat ≤ 1.010), with
   density covered in all of them; the README Results tables are updated. The
   floor check flags the north–south length-scale in 4/5 humpback and 2/5 PWSD
   fits (the 1.28× margin above), with the truth still covered.
5. **Joint simulator check.** `Rscript R/validate_joint_sim.r`. Run on the Mac
   (v1.5 code): the visual halves are fine in both seeds (humpback R² 0.92,
   PWSD 0.55–0.60), and seed 2's eDNA half too (R² 0.96 / 0.63 / 0.60, max
   Rhat 1.01). **Seed 1's eDNA half fails** (hake R² 0.16–0.19, max Rhat
   1.7–1.8, θ_bottle ≈ 1.6 against 10), on main's model as well, and with 4
   chains × 800/800, so it is neither the v1.5 changes nor the short warmup.
   - Seed 1 has extreme hake hot spots (log density up to 9.5, about
     13,000/km²). At its 17 hottest stations every chain sits at hake *and*
     junk copies near zero, even though Ct ≈ 26 and hake has ≈ 99% of the
     reads there.
   - Two features of the likelihood make that state possible: the Ct-noise
     SD grows without limit as copies fall, and at tiny copies the read
     model conditions on presence and its overdispersion collapses.
   - But it is not a better fit: the qPCR alone fits each hot spot about 40
     log-likelihood units better at the true copies. So it looks like a
     sampler trap from the flat starting field.
   - ✅ **Fixed by a data-based starting field** (2026-10-13). Seed 1 refit
     with 4 × 800/800 from a field built from the data (Ct through the
     standard curve, reads relative to junk, smoothed onto the HSGP basis),
     against the same fit started at the true field:

     | Seed 1, eDNA | hake R² | θ_bottle (truth 10) | θ_junk | lp__ | max Rhat |
     |---|---|---|---|---|---|
     | flat start | 0.19 | 1.6 | 1.0 | −15,190 to −15,450 | 1.72 |
     | data-based start | **0.978** | **8.9** | 13.8 | −13,190 to −13,210 | 1.09 (θ_bottle) |
     | start at the truth | 0.978 | 9.0 | 13.6 | −13,180 to −13,230 | 1.06 |

     `fit_joint()` now starts from the data (`edna_field_start()`, eDNA in the
     likelihood only; visual-only fits start as before). Still to do:
     the same start in HSGP4eDNA's `09` driver (used by `fit_edna()` and
     `validate_joint_sim.r`).
6. **Joint comparison.** `R/compare_joint.r` with `joint_conv`, one job per
   seed and mode on stats3 (stage 3; seeds 2 and 3, then 4 and 5, then 1;
   seeds 2 and 3 started before the data-based start). **Seed 1 done** (Mac,
   2026-10-13, data-based start, 4 × 800/800): every fit converged (0
   divergences; max Rhat 1.007 visual, 1.019 joint, 1.06 joint_conv, 1.11
   eDNA-only on θ_bottle). Field R² hake 0.939 in every mode; humpback 0.494
   eDNA / 0.925 visual / **0.928 joint**; PWSD 0.589 / 0.589 / **0.745**;
   conversion factors recovered (humpback 5.30 [5.05, 5.53], truth 5.30; PWSD
   4.77 [4.45, 5.08], truth 4.70). In the explainer (Section 5) and the SMM
   talk. This
   rebuilds the README Phase 5 table and `docs/figures/joint_compare_s<seed>.png`.
   **Replicates (2026-10-21):** the calm-sea seeds 2 and 3 were stopped
   unfinished; the replicate comparison is now the sea-state one (6.2, item
   5), run on seeds 1–5 after the negative-binomial counts (6.2, item 4).
7. **Explainer.** `docs/explainer/01`–`09`: simulated and real fits (eDNA,
   visual at 800 and 1600 functions) and all figures. The simulated and real
   visual fits are done (Mac, 2026-10-12; stage 4's `expl_visual_sim`,
   `expl_real_*`); the eDNA ones and the figures remain. The real eDNA fit
   waits for the 6.1 data. Then regenerate
   `docs/explainer.{md,html}`; several quoted numbers change, e.g. the eDNA
   region's share of the domain, now about 39% of the width and 66% of the
   height.

#### 6.1 Data needed (from collaborators; blocks a real eDNA analysis)

- qPCR Ct values and the standard curve (slope, intercept, detection
  efficiency κ, Ct-noise parameters). The real eDNA fit currently drops all
  hake qPCR.
- Volumes filtered, aliquot and elution volumes, and whether filtered volume
  varies by bottle (the model takes one value).
- Template dilution factors. About 1,160 qPCR replicates were run at 0.1–0.5×
  and metabarcoding at 1–100×, and they are not yet in the likelihood.
- Dates and cruises for the eDNA samples. They and the visual survey are
  treated as one snapshot.
- **MarVer1 amplification efficiencies** from mock communities (per species,
  including hake, humpback and PWSD, ideally with the non-target background),
  and the **PCR cycle number** for each metabarcoding run. Questions for the
  lab:
  - Are there MarVer1 mock-community runs with known template proportions for
    these species, and what were the reads per species?
  - Was the cycle number the same for every library, and is it known per
    sample?
  - How were the 1–100× metabarcoding dilutions chosen, and do diluted and
    undiluted templates go through the same cycle number?

  Why it matters: see 6.2 item 6.
- *(Decision.)* How to treat the 19 sampling depths: an offset, a covariate,
  or separate surfaces.

#### 6.2 Model work (in order)

1. ✅ **HSGP4eDNA posterior predictive reads** (v1.5). `pp_mb_reads` drew
   beta-binomial reads without the zero-inflation (presence) step, so every
   replicate had reads. It now uses `edna_mb_aliquot_reads_rng()`, which draws
   presence first. On a 200-station simulation the share of replicates with
   reads from the new draws matches the likelihood's analytic P(reads > 0) and
   the data (hake 0.803 / 0.803 / 0.807 observed; humpback 0.140 / 0.140 /
   0.142; PWSD 0.185 / 0.184 / 0.188); without the presence step it was 1 for
   every species.
2. **Priors.**
   - **Intercepts: tried, not adopted** (HSGP4eDNA v1.5, opt-in
     `PRIOR_MU=species`). `prior_mu_sp_default()` gives normal(−5, 2) for the
     cetaceans, as the visual and joint models use, and normal(2, 1.5) for the
     rest. On one `bathysp_surface` simulation (budget 400, 4 chains) it made
     the cetacean intercepts *worse*: the data put them slightly low, the
     shared prior's upward pull partly offset that, and normal(−5, 2) adds to
     it. The eDNA default stays the shared normal(2, 1.5); the joint
     formatter keeps its normal(−5, 2) for the cetaceans (now through the same
     helper, values unchanged).

     | Prior | humpback μ (truth −4.61) | PWSD μ (truth −3.91) | PWSD field bias |
     |---|---|---|---|
     | normal(−5, 2) | −4.76 [−5.41, −4.23] | −4.77 [−5.75, −3.96] (missed) | −0.41 |
     | shared normal(2, 1.5) | −4.52 [−5.02, −4.00] | −4.26 [−5.05, −3.46] | −0.25 |

   - **GP marginal SD:** the prior's 99% bound (3) is below the real visual
     fits (2.8–3.6). *(Decision pending: new upper bound, e.g. 5.)*
3. **Junk-noise mixing.** One Phase 5 chain stuck at θ_junk ≈ 2 with the
   centred junk term.
   - ✅ Joint fits now default to **4 chains** (`fit_joint()`,
     `compare_joint.r`), so a stuck chain shows in Rhat.
   - **Non-centred junk term: tried, not adopted** (HSGP4eDNA v1.5, opt-in
     `JUNK_DIST=lognormal`; the species' bottle noise stays exact NB). On one
     `bathysp_surface` simulation (eDNA only, budget 400, 4 chains) recovery
     was the same to two decimals, with no clear gain:

     | Junk term | max Rhat | E-BFMI | runtime | θ_bottle (truth 10) | θ_junk |
     |---|---|---|---|---|---|
     | gamma (default) | 1.08 (θ_bottle) | 0.21–0.27 | 55 min | 6.5 [3.7, 12.1] | 19.2 |
     | lognormal | 1.02 | 0.14–0.18 | 72 min | 6.6 [3.6, 13.5] | 17.6 |

     The stuck chain was intermittent, so one fit cannot test it; the stats3
     joint comparison (4 chains) will show whether it recurs.

4. **Overdispersion for visual segment counts** (*next*, decided 2026-10-21:
   the replicate runs wait for it, so that they use the model that goes
   forward).
   - First simulate clustered sightings (negative-binomial counts) and check
     whether the Poisson model reproduces the real fits' floor-following short
     length-scales, and whether a negative-binomial count term removes them.
   - Then carry it into the visual and joint models: segment counts
     negative binomial around the current Poisson mean, with an
     overdispersion per species and a switch whose off position is exactly
     the Poisson model (checked as in `tests/joint_exactness.R`).
   - Then rerun the sea-state joint comparison on seeds 1–5 (item 5): seeds 1
     and 2 again, on the new model, and seeds 3–5 for the first time.
5. ✅ **Beaufort covariate in the joint model** (merged in #95, 2026-10-21). PWSD
   needs it on real data: AIC on the real distances chose half-normal + group
   size + Beaufort, and the sea-state validation (README) shows that leaving it
   out reads rough sea as low density, shifting PWSD's field by −0.04 to −0.11
   log units per Beaufort unit (up to about 0.6 across the survey). Until now
   the joint model had one ESW per species, and `format_stan_data_joint()`
   warned and dropped `use_bft`. It was left out of Phase 5 to keep the
   simulated case simple.

   **Structure** (ported from `stan/hsgp_visual.stan`, which is validated):
   - **Stan** (`stan/hsgp_joint.stan`): per line-transect species, a switch
     `use_bft_covar[l]` and a slope `beta_bft` on log σ with its own prior
     (the visual model's normal(−0.1, 0.2)):
     log σ = log σ_l + β_size,l (s − s_centre) + β_bft,l (b − 3), where b is
     the segment's Beaufort. Segments are binned to Beaufort levels (0.5
     wide) shared by all species, because they share the segments;
     `seg_bft` / `det_bft` give each segment's and detection's level. The
     size-averaged ESW is computed per species and level (`esw_bft[l]`)
     instead of once per species; segments and detections use their level's
     ESW and σ. A species without the covariate computes one ESW and uses it
     at every level, which is exactly the old model. `p_det` becomes the mean
     over segments.
   - **R**: `format_stan_data_joint()` builds the shared levels from the
     segments' `bft` column and sets `use_bft_covar` per species from its
     `use_bft` setting (instead of warning); `simulate_joint()` gains
     `sea_state = TRUE`, drawing segment Beaufort with the visual simulator's
     `draw_sea_state()` after the field, so the field and eDNA data match the
     `sea_state = FALSE` draw for the same seed; `fit_joint()` starts each
     Beaufort slope near 0; `compare_joint.r` gains `SEA_STATE=1` and a
     choice of which species use the covariate.
   - **Checks** (new `tests/joint_exactness.R`): with the covariate off for
     every species the log density must equal the current joint model's
     exactly; with eDNA off and the covariate on for PWSD it must equal the
     sum of the single-species visual fits (`hsgp_visual.stan`, PWSD with
     `use_bft`) up to a constant. A sea-state case in `tests/smoke.R`.

   **Re-runs:**
   - **Nothing already run changes.** With the covariate off the model is
     identical (the exactness check above), so the calm-sea seed-1 joint
     comparison stands.
   - **New: joint comparison with sea state.** `simulate_joint(sea_state =
     TRUE)`, fitted visual-only, joint and joint with free conversion
     factors, each with the covariate on for both cetaceans and, for
     contrast, off (the eDNA-only fit does not use the visual data and does
     not change). This shows whether the joint model recovers β_bft and
     what leaving it out costs when eDNA is also in the fit.
     - **With the covariate: seeds 1 and 2 done** (Mac, 2026-10-20/21, all
       four modes). Seeds 3–5 were stopped part-way and wait for the
       negative-binomial counts (item 4), after which all five seeds run on
       the new model. Seeds 1–2: 0 divergences, max Rhat
       ≤ 1.06; joint best for both cetaceans (PWSD R² 0.70–0.74 against
       0.59–0.62 for the better single source); humpback β_bft recovered
       closely (−0.068 to −0.079, truth −0.08), PWSD's covered but noisy
       (−0.089 to −0.227, truth −0.15); conversion factors recovered. Tables
       in the explainer, Section 5.
     - **Without the covariate (contrast):** not run yet.
   - ✅ **Docs:** the explainer's Section 5 has the seed 1–2 sea-state
     comparison, and its [REVIEW] note no longer lists the covariate as
     missing; the README points to it.
   - **Real data (6.3):** the real joint fit uses it for PWSD (humpback had no
     Beaufort effect by AIC).
6. **Metabarcoding amplification efficiency** (needs the 6.1 mock-community
   data). Neither the model nor the simulator has species-specific PCR
   efficiency: expected read shares follow aliquot copies relative to junk, so
   MarVer1 is assumed to amplify every species and the junk background
   equally. A species amplifying with efficiency α_s over N cycles gets reads
   inflated by about (1 + α_s)^N (Shelton et al. 2023).
   - **Where it does not matter.** With a fixed cycle number the bias is a
     constant multiplier per species: it shifts a species' whole surface, not
     its shape (field R², length-scales, hot spots unchanged). For hake, qPCR
     fixes the absolute copies and a hake–junk read bias is absorbed by the
     free junk intercept.
   - **Where it does.** The cetaceans' copies are inferred through hake qPCR →
     hake/junk reads → cetacean/junk reads, so a hake–cetacean efficiency
     difference scales their copies by the efficiency ratio. With fixed
     conversion factors (eDNA-only) that biases their absolute density. With
     free conversion factors (joint) it is absorbed but confounded: c_s
     becomes copies per animal × relative amplification bias, so it is not a
     shedding rate and does not transfer to another primer or protocol.
   - **Where it leaks into the pattern.** If the cycle number or dilution
     differs between samples, the bias varies between samples and enters the
     spatial field; check this first (6.1).
   - **How to add it.** A per-species offset log a_s = N log(1 + α_s) on the
     metabarcoding side only (not qPCR), with α_s given a prior from the
     mock-community estimates. Without mock data α_s is not identifiable
     separately from c_s. Simulate it first: draw α_s, fit with and without
     the term, and check what c_s and the cetacean densities recover.
   - Related, from the CalCOFI comparison (`docs/calcofi_hake_comparison.md`):
     test whether copies scale one-for-one with density by freeing an
     elasticity (copies ∝ λ^b); CalCOFI found eDNA amount carried little
     information beyond presence for whales.

#### 6.3 Real joint fit and extensions

- **Joint model on the real data**, with the explainer's predictive checks.
  Needs 6.0's re-runs, 6.2, and ideally 6.1. Joint fits only borrow strength
  where the eDNA bottles are.
- **Remaining scenarios in the joint model.** `bathysp_depthspref` (water
  column + depth preference), and a `bathygp_*` truth as a deliberate
  mis-specification check.
- **Wrap-up.** Write a final README results table. Re-render or retire
  presentations that reference v-numbers.

#### 6.4 Realistic simulation (validation tier 2; planned, details to follow)

Tier 2 of the demonstration (section 2). It sits between the idealised
simulation and the real data: the sampling is real, and the truth is known but
was not generated by the model.

- **Truth: previously modelled density surfaces** for PWSD, humpback and hake,
  put on the 35°N study domain (UTM 10N, km from the domain origin) in the
  model's units (animals per km²; hake's unit to decide).
- **Visual observations** simulated from those surfaces on the **actual 2018
  effort**: the 1,641 segments north of 35°N, with their real lengths and
  Beaufort. Detection, group sizes and truncation come from the real-data
  detection-function fits and the empirical group-size pools.
- **eDNA observations** simulated on the **actual 2019 stations**: the 551
  bottles, with their real locations, sampling depths and replicate structure
  (qPCR and metabarcoding). The conversion factors, volumes, bottle and junk
  noise, and water-column treatment are either the idealised simulator's or
  values informed by the 6.1 data, whichever exist at the time.
- **Fits:** eDNA-only, visual-only, joint and joint with free conversion
  factors, the same four as the idealised comparison. The model is now
  mis-specified (the truth is not a GP plus a spline), so this tests
  robustness as well as sampling.
- **What to measure:**
  - recovery of the true surface on a **prediction grid** over the whole
    domain, not only at sampled locations: R², error and 95% interval
    coverage, and how they degrade away from effort and stations;
  - abundance in regions or strata;
  - recovery of the conversion factors;
  - several replicate observation draws from the same truth.
- **Code needed:**
  - a simulator that takes truth surfaces and the real designs (segments,
    stations, sampling depths from `data/processed/`) and returns the same
    object shape as `simulate_joint()`, so the formatters,
    `fit_joint()` and `compare_joint.r` apply unchanged;
  - prediction on a grid in the joint model (it has none yet; the visual
    model's `pred_coords` is the template);
  - a comparison script that scores the grid predictions.
- **Depends on:** the Beaufort covariate in the joint model (6.2.5; the real
  effort carries Beaufort). Ideally also negative-binomial segment counts
  (6.2.4), if the real sightings are clustered beyond Poisson.
- **To decide (details to follow):**
  - which surfaces, and their source, years, season, resolution and units
    (hake: biomass or numbers, and a conversion to animals per km²);
  - whether to rescale them to plausible eDNA detection rates;
  - how to treat sampling depth (19 real depths) and the one-year gap between
    the 2018 visual and 2019 eDNA surveys (one static truth, or different
    surfaces per survey year).

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

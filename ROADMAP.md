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
3. Write `tests/smoke.R` with the eDNA case: a thin `R/0x_fit_edna.r` that
   calls HSGP4eDNA's simulator, formatter, and `hsgp_nd.stan` via the
   submodule. Also include the existing non-spatial `distance_v4.1` case
   unchanged, as the visual baseline that later phases must preserve. Do not
   include `distance_v4.1a`, which does not sample reliably (see the Phase 0
   results).
4. Rewrite `README.md` down to the target state, marking the parts that are
   not built yet.

✅ Check: the smoke test passes for eDNA and the non-spatial v4.1 visual model. Deleted
files are recoverable from `pre-restructure`.

### Phase 2: Refactor shared Stan functions (in HSGP4eDNA)

1. Move the functions block of `hsgp_nd.stan` and `hsgp_2d_bathysp.stan` into
   `stan/include/hsgp_functions.stan` and `stan/include/edna_functions.stan`.
   Optionally move the eDNA likelihood into an `edna_obs_lp(...)` function, so
   the joint model calls exactly the same code.
2. **(Required.)** Fold `hsgp_2d_bathysp.stan` into `hsgp_nd.stan` behind a
   `K_bathy` spline-basis count (0 = off). Make the spline's design-matrix
   construction an include as well, so the visual and joint models share it.
   HSGP4eDNA then also has a single eDNA model.
3. Turn `log_lik_qpcr` / `log_lik_mb` back on in generated quantities (L4 needs
   them anyway, for model comparison).

✅ Check: HSGP4eDNA's matched fits reproduce with identical posteriors for a
fixed seed, or within MC error. Re-run at least the `nobathy_surface × gp2d`
and `bathysp_surface × gp2d_bathysp` cells. Tag `v1.1` and bump the submodule.

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
   Phase 1 → Phase 2 (makes the spline model includable) → Phase 3 (visual).
4. **Phase 2 changes HSGP4eDNA itself** (includes, the optional model merge,
   `log_lik`). That repo has its own validation table and collaborators. If
   you'd rather leave it untouched, the only fallback is to vendor its Stan
   functions into this repo's `stan/include/`, plus a smoke-test check that
   they still match upstream. A monolithic `.stan` file can't be `#include`d
   piecemeal.

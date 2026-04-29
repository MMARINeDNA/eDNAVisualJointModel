# eDNAVisualJointModel

Joint modelling of visual survey and environmental DNA (eDNA) data for
marine megafauna off the US West Coast. The workflow simulates realistic
eDNA observations (qPCR + metabarcoding) for three species on a spatial
grid, and fits an anisotropic Hilbert-Space GP (HSGP) model to the
simulated data in Stan.

The repo also contains a parallel **distance-sampling** pipeline (the
visual-survey side of the joint model, fit standalone for development
and PPC), self-contained **notebooks** that walk through each model
version, a **GP presentation** suitable for sharing with collaborators,
and a small set of **data figures** summarising the real survey data.

## Directory layout

```
.
├── README.md                           This file
├── LICENSE
├── eDNAVisualJointModel.Rproj
│
├── 00_pipeline_v1.r                    End-to-end orchestrators (run from root)
├── 00_pipeline_v2.r                    `Rscript 00_pipeline_v{N}.r` sources
├── 00_pipeline_v3.r                     each step in order.
├── 00_pipeline_v3.1.r                  v3 + reparameterised model fit
├── 00_pipeline_v3.2.r                  v3 hake-only / surface / qPCR-only debug
├── 00_pipeline_v4.r
├── 00_pipeline_v4.1.r                  v4 + v3.2 lessons + zero-mean GP
│
├── stan/                               Stan model source
│   ├── whale_edna_hsgp_v1.stan
│   ├── whale_edna_hsgp_v2.stan
│   ├── whale_edna_hsgp_v3.stan
│   ├── whale_edna_hsgp_v3.1.stan       reparameterised v3
│   ├── whale_edna_hsgp_v3.2.stan       v3 hake-only / qPCR-only
│   ├── whale_edna_hsgp_v4.stan
│   └── whale_edna_hsgp_v4.1.stan       v4 + v3.2 fixes
│
├── scripts/                            R pipeline steps (sim → plot → format → run → check)
│   ├── 01_simulate_whale_edna_v1.r         simulate eDNA data
│   ├── 01_simulate_whale_edna_v2.r
│   ├── 01_simulate_whale_edna_v3.r
│   ├── 01_simulate_whale_edna_v3.1.r
│   ├── 01_simulate_whale_edna_v3.2.r
│   ├── 01_simulate_whale_edna_v4.r
│   ├── 01_simulate_whale_edna_v4.1.r
│   ├── 02_plot_simulated_data_v2.r         plot the simulated truth
│   ├── 02_plot_simulated_data_v3.r
│   ├── 02_plot_simulated_data_v3.1.r
│   ├── 02_plot_simulated_data_v3.2.r
│   ├── 02_plot_simulated_data_v4.r
│   ├── 02_plot_simulated_data_v4.1.r
│   ├── 03_format_stan_data_v1.r            assemble stan_data list
│   ├── 03_format_stan_data_v2.r
│   ├── 03_format_stan_data_v3.r
│   ├── 03_format_stan_data_v3.1.r
│   ├── 03_format_stan_data_v3.2.r
│   ├── 03_format_stan_data_v4.r
│   ├── 03_format_stan_data_v4.1.r
│   ├── 04_run_whale_edna_model_v1.r        compile + fit Stan model
│   ├── 04_run_whale_edna_model_v2.r
│   ├── 04_run_whale_edna_model_v3.r
│   ├── 04_run_whale_edna_model_v3.1.r
│   ├── 04_run_whale_edna_model_v3.2.r
│   ├── 04_run_whale_edna_model_v4.r
│   ├── 04_run_whale_edna_model_v4.1.r
│   ├── 05_check_whale_edna_model_v1.r      diagnostics, PPC, recovery plots
│   ├── 05_check_whale_edna_model_v2.r
│   ├── 05_check_whale_edna_model_v3.r
│   ├── 05_check_whale_edna_model_v3.1.r
│   ├── 05_check_whale_edna_model_v3.2.r
│   ├── 05_check_whale_edna_model_v4.r
│   ├── 05_check_whale_edna_model_v4.1.r
│   ├── DetectionsBySpecies.R               Ad-hoc data exploration
│   └── FinWhales.R
│
├── outputs/                            Generated artifacts (per-version)
│   ├── whale_edna_output_v{1,2,3,3.1,3.2,4,4.1}/   Per-version eDNA output:
│   │       whale_edna_sim_v{N}.rds           Sim output (gitignored)
│   │       simulated_edna_fields_v{N}.pdf    Tracked multi-page sim plot
│   │       stan_data.rds                     stan_data list from step 03 (gitignored)
│   │       whale_edna_fit.rds                CmdStanR fit object from step 04 (gitignored)
│   │       *.png, *.csv, session_info.txt   diagnostics from step 05
│   │       v{N}_notebook.qmd                 source for the per-version write-up
│   └── distance_v{N}/                  Distance-sampling per-version output (see below)
│
├── distance/                           Line-transect distance-sampling pipeline
│   ├── 00_distance_v{N}.R                  Self-contained simulate → fit → diagnose
│   ├── 00_distance_v4.1a.R                 Spatial HSGP variant (cetaceans only)
│   ├── distance_hn_dens_v{N}.stan          Non-spatial HN detection + density model
│   ├── distance_v4.1a.stan                 Spatial HSGP version of the LT model
│   ├── humpback_grpsz.RData                Empirical group-size data used by sim
│   └── pwsd_grpsz.RData
│
├── notebooks/                          Self-contained HTML model write-ups
│   ├── README.md                           Index + regeneration instructions
│   ├── v3_notebook.html                    v3 HSGP joint qPCR / MB fit
│   ├── v3.2_notebook.html                  v3.2 debugging case study
│   └── distance_v4.1_notebook.html         Distance-sampling pipeline (PRs #34–#37)
│
├── presentations/                      Quarto slide decks (RevealJS) + presenter notes (HTML)
│   ├── gp_models.{qmd,html}                              GP intro deck
│   ├── gp_models_presenter_notes.{qmd,html}              ...with detailed companion notes
│   ├── project_goals.{qmd,html}                          Project goals + state-space overview
│   └── project_goals_presenter_notes.{qmd,html}          ...with detailed companion notes
│
├── figures/                            Real-data summary maps (PNG)
│
├── Data/                               Real survey data
│       effort.csv, sightings.csv            Line-transect (CCES 2018)
│       hake_qPCR_MURI_df.csv                qPCR (tracked, ~200 KB)
│       MV1_MURI_df.csv                      MARVER1 metabarcoding (gitignored, ~30 MB)
│
└── scripts/                            (continued) data-summary plot scripts
        plotLTData.R, ploteDNAData.R         multi-species + multi-method overview maps
        plotPWSDData.R, plotHumpbackData.R   per-species LT + eDNA two-panel figures
```

## Pipeline

Each version has a five-step pipeline. Run end-to-end with:

```
Rscript 00_pipeline_v4.r
```

…or run individual steps in order:

```
Rscript scripts/01_simulate_whale_edna_v4.r     # -> outputs/whale_edna_output_v4/whale_edna_sim_v4.rds
Rscript scripts/02_plot_simulated_data_v4.r     # -> outputs/whale_edna_output_v4/simulated_edna_fields_v4.pdf (3 pages)
Rscript scripts/03_format_stan_data_v4.r        # -> outputs/whale_edna_output_v4/stan_data.rds
Rscript scripts/04_run_whale_edna_model_v4.r    # -> outputs/whale_edna_output_v4/whale_edna_fit.rds
Rscript scripts/05_check_whale_edna_model_v4.r  # -> diagnostics in outputs/whale_edna_output_v4/
```

Each script expects the **project root as the working directory**. Artifacts
(`.rds`, `whale_edna_output_v*/`) are `.gitignored` — they are regenerable.

## Versions

### V1 — initial simulation + model

- **Domain**: Oregon / Washington coast only (UTM 10N, ~300 × 400 km).
- **Bathymetry**: logistic shelf-slope profile as a pure function of
  Easting `X` (isobaths run due N-S).
- **Observation design**: 150 stations × 5 sample depths (0, 50, 150, 300,
  500 m); marker MARVER3 at 5000 reads/replicate; `vol_filtered = 2.0 L`.
- **Species structure**: bathymetric habitat preference only (no latitude
  effect); shelf-slope-break hake, shallow-shelf humpback, deep-slope PWSD.

### V2 — code review refinements (with Ole)

- Domain and bathymetry: unchanged from V1.
- Marker switched to MARVER1; read depth 50 000; `vol_filtered = 2.5 L`;
  explicit aliquoting step (`vol_aliquot = 2 µL` of a 100 µL elution).
- Explicit animals → copies conversion (`conv_factor = 10`).
- `zsample_pref` values shifted so baseline eDNA multiplier is ≈ 1× rather
  than suppressing most depths.
- Vectors carry an `_si` suffix (sample × species index order).
- Stan model (V2): qPCR Ct hurdle on hake, zero-inflated Beta-Binomial
  metabarcoding for all three species, `phi` parameterised in log total
  copies.

### V3 — extended domain, rotated bathymetry, realistic species distributions

(Files: `*_v3.*`. Original v3 model fit, kept alongside v3.1 for
side-by-side comparison.)

### V3.1 — v3 + reparameterised model fit

(Files: `*_v3.1.*`. Same simulation, plotting, and downstream
structure as v3, but the model fit chain is reparameterised to fix
the v3 sampler pathology — `HSGP_M = c(5, 5, 5)` instead of
`c(10, 8, 8)`; the BB-phi `fmax(..., 0)` hinge replaced with
`log1p_exp`; `kappa` fixed as data instead of sampled; `gp_l` and
`gamma_phi` priors substantially tightened. Outputs land in
`outputs/whale_edna_output_v3.1/`
so the two pipelines don't collide.)

### V3.2 — hake-only, surface-only, qPCR-only debug case

(Files: `*_v3.2.*`. Stripped-down debug version of v3 to isolate
the HSGP latent-field recovery problem from the metabarcoding
likelihood. Same domain, bathymetry, station design, and hake
habitat preference as v3.)

- **Species**: Pacific hake only (S=1). Humpback / PWSD removed.
- **Sample depth**: surface only (`Z_sample = 0`). The water-column
  multiplier collapses to 1 everywhere, so `log_zsample_effect` is all
  zeros and contributes nothing to the likelihood.
- **Sim matches model exactly** (Option A). The deterministic habitat
  preference function (`zbathy_pref` + `y_pref`) that v3 added to the
  GP mean has been removed. The simulated `f_s` is now a pure
  zero-mean GP draw with covariance `K(lx, ly, lz)`, so the named
  length-scales are unambiguously the truth and length-scale recovery
  is a clean test.
- **Observations**: qPCR only. The metabarcoding likelihood and all
  MB-only parameters (`beta0_phi` / `gamma0_phi` / `gamma1_phi`) are
  removed from the Stan model.
- **Fixed qPCR calibration**: `kappa` and `sigma_ct` both promoted
  from sampled parameters to data (joining `alpha_ct` / `beta_ct`).
  The standard-curve fit is treated as known. `sigma_ct` was
  consistently inflated when sampled (~0.79 vs truth 0.4) because of
  a `Ct ~ log(integer count)` vs `Ct ~ log(expected count)`
  discreteness mismatch between sim and model; pinning it removes
  that nuisance.
- **Simulation `sigma_ct`** = 0.4 (the explicit PCR-noise SD added
  in the sim).
- **Model `sigma_ct`** = **0.6** (data-fixed, deliberately *higher*
  than the sim's PCR-noise truth). The sim's Ct values also carry
  discrete-count noise from `Ct ~ log(integer aliquot count)`, which
  contributes ~0.46 of additional SD at the simulated mean aliquot
  copies. Pinning `sigma_ct = 0.4` (PCR-only) made the model fit the
  field too tightly to absorb that extra noise — 50% max-treedepth
  saturation, length-scales biased low. 0.6 ≈ the honest
  `√(0.4² + 0.46²)` combined PCR + discrete-count noise.
- **`max_treedepth`** bumped from 12 to 14 — the geometry is still
  somewhat stiff with `M = 3584` even after the `sigma_ct` bump, and
  the extra trajectory headroom prevents truncated transitions.
- **Sampled parameters**: only `mu_sp`, `gp_sigma`, `gp_l`, `z_beta`.
  With kappa, sigma_ct, and the MB block all out of the parameter
  block, the only things the posterior has to identify are `mu_sp`
  and the latent GP.
- **`log_vol_filtered` in the eDNA log-mean**: the previous v3.2 model
  computed `log_lambda_edna = log_lambda + log_zsample_effect +
  log_conv_factor`, missing the `log(vol_filtered) ≈ 0.92` factor that
  the simulation uses when generating bottle copies
  (`E[bottle] = conv_factor · λ · zsample_effect · vol_filtered`). This
  forced the posterior `mu_sp` upward by ~0.7-0.9 to compensate; now
  passed in as data and included in `log_lambda_edna`.
- **`gp_l` priors centred on the v3 truth**: the previous priors `gp_lx
  ~ N(50, 40)`, `gp_ly ~ N(150, 80)`, `gp_lz ~ N(300, 150)` had `ly` and
  `lz` mu's at half / double the simulated truth (left-overs from an
  earlier sim configuration). Now `gp_lx ~ N(50, 30)`, `gp_ly ~ N(300,
  100)`, `gp_lz ~ N(150, 50)` — centred on truth and ~30% tighter to
  fight the multi-modality the previous mismatch was producing
  (Rhat 1.4+ on `gp_l[1,1]` and `gp_l[1,3]`).
- **HSGP basis count**: `M = (14, 8, 32)`. Bumped from `(10, 8, 8)`
  in v3 → `(10, 8, 20)` after the normalisation fix → `(14, 8, 32)`
  paired with Option A. Sized to Riutort-Mayol's faithful-
  representation rule `m >= 1.75 c / ρ_ℓ`, which with `c = 1.5` and
  the corrected half-ranges (250 km, 635 km, 1750 m) requires
  `m >= (13, 6, 31)` for the true `(lx, ly, lz) = (50, 300, 150)`.
- **`gp_sigma` prior**: `gamma(8, 4)` (mode 1.75, mean 2, sd 0.71).
  Tightened from `gamma(4, 2)` after the latter was found to be bimodal
  across chains (some chains got trapped in a low-`gp_sigma` mode where
  `f ≈ 0` and `sigma_ct` absorbed the variance). With `M = 3584` the
  posterior has a near-degenerate `gp_sigma ↔ z_beta` ridge; the
  tighter prior makes the low-`gp_sigma` mode essentially unreachable.
  Original switch from `half_normal(0, 1.5)` to gamma was needed
  because half-normal's mode at zero was visibly trapping `gp_sigma`
  at the boundary; gamma puts zero density at 0.
- **Coordinate normalisation**: `coord_centre` / `coord_scale` derived
  from the actual v3 domain extents (`X_km_max / 2`, `Y_km_max / 2`,
  `3500 / 2`) so all normalised coords land in `[-1, 1]`.
  Pre-existing bug carried over from v1/v2 fixed.
- Outputs land in `outputs/whale_edna_output_v3.2/`.

- **Domain**: extended to cover the full US West Coast — San Francisco
  (37.77°N) to the US/Canada border (~49°N), 500 km × 1270 km in UTM 10N.
- **Bathymetry**: rotated cross-shore axis `X_prime` (25° in normalised
  space) so isobaths run NW-SE. A pure 45° rotation confounds latitude
  with bathymetry in this elongated domain.
- **Observation design**: **300 stations** stratified 50/30/20 shelf /
  slope / offshore on the rotated coordinate; **3 sample depths** (0, 150,
  500 m); samples with `Z_sample > Z_bathy` dropped.
- **Species spatial structure**:
  - Pacific hake — shelf-slope break, broad in latitude.
  - Humpback whale — southern + inshore shelf.
  - Pacific white-sided dolphin — northern + offshore slope.
- **Model form**: bathymetric + latitude habitat structure is absorbed
  into the GP via a non-zero GP mean, so the exposed model is
  `log(λ) = μ + f`.
- **Plotting**: the 1 km grid in `02_plot_simulated_data_v3.1.r` evaluates
  the closed-form GP mean at every cell (501 × 1271 = 637k cells) — no
  kriging.

### V4 — variable per-sample replication + junk reads + retuned densities

- **Sampling design**: most samples still have the full 3 qPCR + 3
  metabarcoding replicates, but ~20% of samples (chosen independently for
  each process) are reduced to 1 or 2 reps, reflecting realistic
  per-sample loss / QC failure. Controlled by `prop_reduced_qpcr_rep` and
  `prop_reduced_mb_rep` in `scripts/01_simulate_whale_edna_v4.r`.
- **Data layout**: the v4 sim now emits the qPCR and metabarcoding
  observed data already in **long form** — one row per aliquot, with
  explicit `qpcr_sample_idx` / `mb_sample_idx` vectors and per-sample
  replicate counts `n_qpcr_rep_i` / `n_mb_rep_i`. The v4 format script
  is a thin pass-through to the Stan data list.
- **"Junk" MB background**: real MARVER1 data is dominated by
  non-target species (plankton, bacteria, other vertebrates, primer
  artifacts). v4 models this with a lumped junk category that claims a
  random 5–95% of each MB aliquot's read depth, leaving only the
  remaining reads to be distributed among the three target species.
  Junk does not appear in qPCR (species-specific primers). Stored as
  `mb_junk_reads` and `mb_pi_junk` in the sim output for transparency;
  the Stan model sees only target reads.
- **Realistic MB read-depth distribution**: per-aliquot total read
  depth (junk + target) is drawn from a two-component log-normal
  mixture — 80% of aliquots from a tight "typical run" around 75k
  reads (sdlog 0.25), 20% from a wider "problem run" component
  (meanlog log(40000), sdlog 1.0). Hard-clipped to `[1000, 250000]`.
  In practice ~70% of aliquots fall in 50k–100k reads, with some
  dropping below 10k or exceeding 200k — matching what we see across
  real MARVER1 sequencing runs. Parameters live at the top of
  `scripts/01_simulate_whale_edna_v4.r` as `mb_reads_*`.
- **Species-specific conversion factor**. Whales shed far more eDNA
  per animal than hake, so `conv_factor` is a length-S vector rather
  than a scalar. Current values (copies/L per animal/km² per litre
  filtered): hake = 10, humpback = 125, PWSD = 40. Passed to Stan as
  `log_conv_factor[S]` (vector) and indexed by species in the
  `log_lambda_edna[i, s] = log_lambda[i, s] + log_zsample_effect[i, s]
  + log_conv_factor[s]` computation.
- **Realistic animal densities** (animals/km²). Mean λ across the
  study area: hake ≈ 12, humpback ≈ 0.005, PWSD ≈ 0.05 — whales are
  orders of magnitude rarer than hake at real field densities.
  Controlled by `mu` values: hake = log(6), humpback = log(0.004),
  PWSD = log(0.032).
- **Detection rates preserved**. Despite the density change, sample-
  level detection still matches the real-data targets: qPCR hake
  ≈ 95% / MB hake ≈ 100% / MB humpback ≈ 20% / MB PWSD ≈ 20% —
  because the species-specific conv_factor rescales so that expected
  bottle copies (cf × λ × vol_filtered) stays in the same ballpark.
  Humpback and PWSD distributions remain independent. Sample-level
  detection rates are printed at the end of the sim run.
- **Stan model**: only change from v3 is `log_conv_factor` going from
  scalar to `vector[S]`. The long-form `*_sample_idx` structure from
  v3 already accommodates variable replication per sample.
- **Domain, bathymetry, species habitat structure, priors**: unchanged
  from v3.

### V4.1 — v4 + v3.2 lessons + zero-mean GP

(Files: `*_v4.1.*`. v4 simulation upgraded with the lessons learned
in v3.2 / v3.2-bump-sigma-ct-and-treedepth. Output dir:
`outputs/whale_edna_output_v4.1/`.)

- **Sampling design**: 500 stations × 6 sample depths
  (0, 50, 150, 250, 350, 500 m), up from 300 × 3. Samples with
  `Z_sample > Z_bathy` dropped, leaving N ≈ 2000–2400 in practice.
- **Sim matches model exactly (Option A)**. Habitat-preference fields
  (`zbathy_pref_*`, `y_pref_*`) deleted; the simulated `f_s` is a
  pure zero-mean GP draw with covariance `K(lx, ly, lz)`. Length-scale
  recovery against the named scales is a clean test.
- **MB junk reads still simulated** for read-depth realism, but —
  same as v4 — the Stan model only sees the 3 target species
  (`mb_total` is target-reads-only).
- **Stan model**: kappa, sigma_ct moved from `parameters` to `data`;
  `log_vol_filtered` added to `log_lambda_edna`; `gp_sigma` prior
  switched from `half_normal` to `gamma(8, 4)`; `prior_kappa_*` /
  `prior_sigma_ct_*` removed from data block.
- **Format script**: `coord_centre`/`coord_scale` derived from sim
  domain extents (so all normalised coords sit in `[-1, 1]`);
  `HSGP_M = (14, 8, 32)` per Riutort-Mayol's threshold for the
  named length-scales; `gp_l` priors centred on truth
  (`lx ~ N(50, 30)`, `ly ~ N(300, 100)`, `lz ~ N(150, 50)`);
  `sigma_ct = 0.7` (inflated above the sim's 0.4 to absorb
  discrete-count noise from `Ct ~ log(integer count)`).
- **Run script**: `iter_warmup = 1000`, `iter_sampling = 1000`,
  `max_treedepth = 14`. kappa / sigma_ct dropped from `init_fn`.
- **Check script**: scalar trace plots / parameter-recovery scatter
  drop kappa/sigma_ct (no longer sampled). New per-species
  posterior-decomposition section produces 1D and 2D marginals of the
  latent field. New per-species prior-vs-posterior density overlays
  for `mu_sp`, `gp_sigma`, `gp_l[1..3]` (15 panels: 3 species × 5
  parameters).

## Distance sampling

`distance/` holds a self-contained line-transect distance-sampling
pipeline. It is the **visual-survey half** of the joint model, factored
out of the joint Stan model so the detection / encounter-rate logic
can be developed, tested, and PPC'd standalone before being merged
back. The same v4.1 zero-mean GP is used to draw the underlying
density surface, so anything learned here transfers directly to the
joint pipeline.

Per species (humpback whale, Pacific white-sided dolphin) one run does:

1. Simulate a zero-mean GP **animal-density** surface (animals/km²)
   with the v4.1 domain, bathymetry, and species-specific GP
   hyperparameters.
2. Convert to a **group-density** surface using a per-species mean
   group size (humpback = 2, PWSD = 50).
3. Lay down systematic E-W parallel transects, split into 10 km
   segments, and simulate sightings from a **species-specific
   half-normal detection function**. PWSD's detection function has a
   group-size covariate, `log σ_g = log σ_0 + β_size · (s_g − s̄)`,
   so larger groups are more detectable.
4. Draw PWSD group sizes from a log-normal fit to the empirical
   `pwsd_grpsz.RData`; humpback group sizes are constant at 2.
5. Fit `distance/distance_hn_dens_v4.1.stan` (one fit per species).
   The model includes a Jensen-corrected ESW, the size-bias correction
   `+ Σ log(esw_i) − n · log(esw_pop)` in the group-size likelihood,
   and accepts species-specific priors on `log_lambda_s` as data.
6. Write diagnostic plots, a parameter-recovery table, and
   prior-vs-posterior density overlays for animal density `D = λ_s·µ_s`
   to `outputs/distance_v4.1/`.

Run from the project root:

```
Rscript distance/00_distance_v4.1.R
```

The full development arc (PRs #34–#37) is documented in
`notebooks/distance_v4.1_notebook.html`.

### v4.1a — spatial HSGP variant

`distance/distance_v4.1a.stan` + `distance/00_distance_v4.1a.R` extend
v4.1 by replacing the single scalar `log_lambda_s` parameter with a 3-D
anisotropic HSGP latent field over segment locations, mirroring the
spatial structure of `stan/whale_edna_hsgp_v3.2.stan`:

- `log λ_groups(x, y, z) = µ_sp + f(x, y, z)`, `f ∼ GP(0, K(ℓx, ℓy, ℓz))`
- HSGP basis `M = (14, 8, 32) = 3584` (same as v4.1), boundary `c = 1.5`,
  per-species `gp_l` priors centred on the named length-scales
  ((50, 300, 100) for humpback; (40, 300, 300) for PWSD).
- Detection function, group-size sub-model, Jensen-corrected ESW, and
  size-bias correction are identical to v4.1.
- Generated quantities adds spatial PPCs: per-segment `lambda_groups`,
  `f_seg` (spatial component), posterior-predictive segment counts, and
  per-segment / per-detection log-likelihoods for `loo`.

Outputs land in `outputs/distance_v4.1a/`. Each species-fit produces a
six-panel PDF with the truth + posterior `λ` surfaces, per-segment
truth-vs-posterior scatters for `λ` and `f`, the count PPC, and a GP
hyperparameter-recovery facet.

Run from the project root (full fit takes ~1 hr per species):

```
Rscript distance/00_distance_v4.1a.R
# or for a quick wiring smoke-test:
ITER_WARMUP=200 ITER_SAMPLING=200 Rscript distance/00_distance_v4.1a.R
```

## Notebooks

`notebooks/` collects self-contained HTML write-ups, one per model
version, that explain the simulation, the statistical model (with
equations), and the diagnostic + posterior-predictive output. The
HTML files are inlined (Quarto's `embed-resources: true`) so they can
be opened in a browser, dropped into Slack, or emailed without
bringing additional assets along.

| File | Covers |
|---|---|
| `v3_notebook.html` | v3 simulation + HSGP joint qPCR / metabarcoding fit. |
| `v3.2_notebook.html` | v3.2 debugging case study — eight PRs of iterative diagnosis on the v3 sampler pathology, walked through chronologically. |
| `distance_v4.1_notebook.html` | Line-transect distance sampling pipeline (PRs #34–#37). Spatial GP density → group-density surface → simulated sightings via half-normal detection function → non-spatial Stan distance fit, with per-species PPC plots and parameter recovery. |

Each notebook's source `.qmd` lives next to its model artefacts in the
matching `outputs/whale_edna_output_v{N}/` (or `outputs/distance_v{N}/`)
folder. See `notebooks/README.md` for regeneration instructions and
the `pdftoppm` recipe used to extract the per-version "simulated
truth" page-by-page figures from the multi-page sim PDFs.

## Presentations

`presentations/` holds Quarto / RevealJS slide decks aimed at
collaborators / audience-members rather than at people running the
code.

| File | Covers |
|---|---|
| `gp_models.html` | Practical introduction to Gaussian processes for spatial modelling: what a GP is, GPs as priors over spatial fields, how to fit a GP, and the HSGP approximation (Riutort-Mayol et al. 2023) that this project actually uses. |
| `project_goals.html` | Project goals + conceptual diagram of the joint state-space model (single latent density surface → eDNA + line-transect observation submodels). Brief slides on the GP, qPCR/metabarcoding, and LT components, plus the per-species LT + eDNA data figures for humpback and PWSD. |

Each deck has a **companion presenter-notes HTML document** at
`<deck>_presenter_notes.html`. Same section headers as the slides, but
in plain HTML (not RevealJS) with the longer-form prose, equations, and
references that don't fit on the slides. Use the deck in a room and
the notes for self-study or while preparing the talk.

Render either with `quarto render presentations/<file>.qmd`; the
rendered HTML for both formats is tracked in the repo so they're easy
to share without a Quarto install.

## Data figures

A handful of standalone scripts in `scripts/` produce summary maps of
the **real** survey data (no model fitting). They share a consistent
visual style (light-grey land, grey effort lines or dark crosses for
sampling locations, red sized circles for sightings/detections) and
all write to `figures/`.

| Script | Output | Contents |
|---|---|---|
| `scripts/plotLTData.R` | `figures/lt_data_pwsd_humpback.png` | Two-panel LT effort + sightings: PWSD (sqrt-scaled group sizes 1 – 452), humpback (linear 1 – 8). |
| `scripts/ploteDNAData.R` | `figures/edna_data_4panel.png` | Four-panel eDNA: hake qPCR, hake / PWSD / humpback MARVER1. Crosses = sampled locations; circles sized by qPCR positives or total target reads. |
| `scripts/plotPWSDData.R` | `figures/lt_edna_pwsd.png` | PWSD: LT (left) + eDNA (right) side by side, independent legends. |
| `scripts/plotHumpbackData.R` | `figures/lt_edna_humpback.png` | Humpback: LT (left) + eDNA (right) side by side, independent legends. |

All four scripts run from the project root with `Rscript scripts/<name>.R`.

## Versioning convention

- Every file in `scripts/`, `stan/`, and the root `00_pipeline_v{N}.r`
  has an explicit `_v{N}` suffix (lowercase `v`).
- V{N} simulation writes `outputs/whale_edna_output_v{N}/whale_edna_sim_v{N}.rds`.
- V{N} plotting reads `outputs/whale_edna_output_v{N}/whale_edna_sim_v{N}.rds`
  and writes
  `outputs/whale_edna_output_v{N}/simulated_edna_fields_v{N}.pdf`
  (multi-page PDF, one panel
  per page).
- V{N} Stan-data formatting writes
  `outputs/whale_edna_output_v{N}/stan_data.rds`.
- V{N} model runner writes
  `outputs/whale_edna_output_v{N}/whale_edna_fit.rds` and CmdStan files
  to the same directory.
- V{N} checking script writes all diagnostic plots and CSVs into
  `outputs/whale_edna_output_v{N}/`.

When starting a new iteration (`v5`), copy the highest-version scripts
(including `00_pipeline_v4.r`) to their `_v5` counterparts and update
this README with what changed and why.

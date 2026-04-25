# eDNAVisualJointModel

Joint modelling of visual survey and environmental DNA (eDNA) data for
marine megafauna off the US West Coast. The workflow simulates realistic
eDNA observations (qPCR + metabarcoding) for three species on a spatial
grid, and fits an anisotropic Hilbert-Space GP (HSGP) model to the
simulated data in Stan.

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
├── 00_pipeline_v4.r
│
├── stan/                               Stan model source
│   ├── whale_edna_hsgp_v1.stan
│   ├── whale_edna_hsgp_v2.stan
│   ├── whale_edna_hsgp_v3.stan
│   ├── whale_edna_hsgp_v3.1.stan       reparameterised v3
│   └── whale_edna_hsgp_v4.stan
│
├── scripts/                            R pipeline steps (sim → plot → format → run → check)
│   ├── 01_simulate_whale_edna_v1.r         simulate eDNA data
│   ├── 01_simulate_whale_edna_v2.r
│   ├── 01_simulate_whale_edna_v3.r
│   ├── 01_simulate_whale_edna_v3.1.r
│   ├── 01_simulate_whale_edna_v4.r
│   ├── 02_plot_simulated_data_v2.r         plot the simulated truth
│   ├── 02_plot_simulated_data_v3.r
│   ├── 02_plot_simulated_data_v3.1.r
│   ├── 02_plot_simulated_data_v4.r
│   ├── 03_format_stan_data_v1.r            assemble stan_data list
│   ├── 03_format_stan_data_v2.r
│   ├── 03_format_stan_data_v3.r
│   ├── 03_format_stan_data_v3.1.r
│   ├── 03_format_stan_data_v4.r
│   ├── 04_run_whale_edna_model_v1.r        compile + fit Stan model
│   ├── 04_run_whale_edna_model_v2.r
│   ├── 04_run_whale_edna_model_v3.r
│   ├── 04_run_whale_edna_model_v3.1.r
│   ├── 04_run_whale_edna_model_v4.r
│   ├── 05_check_whale_edna_model_v1.r      diagnostics, PPC, recovery plots
│   ├── 05_check_whale_edna_model_v2.r
│   ├── 05_check_whale_edna_model_v3.r
│   ├── 05_check_whale_edna_model_v3.1.r
│   ├── 05_check_whale_edna_model_v4.r
│   ├── DetectionsBySpecies.R               Ad-hoc data exploration
│   └── FinWhales.R
│
├── outputs/                            Generated artifacts
│   ├── simulated_edna_fields_v2.pdf        Tracked plots of simulated data
│   ├── simulated_edna_fields_v3.pdf        (multi-page PDF, one panel per page)
│   ├── simulated_edna_fields_v3.1.pdf
│   ├── simulated_edna_fields_v4.pdf
│   ├── whale_edna_sim_v{1,2,3,3.1,4}.rds   Sim outputs (gitignored)
│   └── whale_edna_output_v{1,2,3,3.1,4}/   Stan fit artifacts (gitignored):
│           stan_data.rds                     stan_data list from step 03
│           whale_edna_fit.rds                CmdStanR fit object from step 04
│           *.png, *.csv, session_info.txt   diagnostics from step 05
│
├── Data/                               Real survey data (effort, sightings, detections)
└── Distance sampling trials/           Self-contained DS exploration (separate Stan model)
```

## Pipeline

Each version has a five-step pipeline. Run end-to-end with:

```
Rscript 00_pipeline_v4.r
```

…or run individual steps in order:

```
Rscript scripts/01_simulate_whale_edna_v4.r     # -> outputs/whale_edna_sim_v4.rds
Rscript scripts/02_plot_simulated_data_v4.r     # -> outputs/simulated_edna_fields_v4.pdf (3 pages)
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
`outputs/whale_edna_output_v3.1/` and `outputs/whale_edna_sim_v3.1.rds`
so the two pipelines don't collide.)

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

## Versioning convention

- Every file in `scripts/`, `stan/`, and the root `00_pipeline_v{N}.r`
  has an explicit `_v{N}` suffix (lowercase `v`).
- V{N} simulation writes `outputs/whale_edna_sim_v{N}.rds`.
- V{N} plotting reads `outputs/whale_edna_sim_v{N}.rds` and writes
  `outputs/simulated_edna_fields_v{N}.pdf` (multi-page PDF, one panel
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

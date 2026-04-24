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
│
├── stan/                               Stan model source
│   ├── whale_edna_hsgp_v1.stan
│   ├── whale_edna_hsgp_v2.stan
│   └── whale_edna_hsgp_v3.stan
│
├── scripts/                            R pipeline steps (sim → plot → format → run → check)
│   ├── 01_simulate_whale_edna_v1.r         simulate eDNA data
│   ├── 01_simulate_whale_edna_v2.r
│   ├── 01_simulate_whale_edna_v3.r
│   ├── 02_plot_simulated_data_v2.r         plot the simulated truth
│   ├── 02_plot_simulated_data_v3.r
│   ├── 03_format_stan_data_v1.r            assemble stan_data list
│   ├── 03_format_stan_data_v2.r
│   ├── 03_format_stan_data_v3.r
│   ├── 04_run_whale_edna_model_v1.r        compile + fit Stan model
│   ├── 04_run_whale_edna_model_v2.r
│   ├── 04_run_whale_edna_model_v3.r
│   ├── 05_check_whale_edna_model_v1.r      diagnostics, PPC, recovery plots
│   ├── 05_check_whale_edna_model_v2.r
│   ├── 05_check_whale_edna_model_v3.r
│   ├── DetectionsBySpecies.R               Ad-hoc data exploration
│   └── FinWhales.R
│
├── outputs/                            Generated artifacts
│   ├── simulated_edna_fields_v2.png        Tracked plots of simulated data
│   ├── simulated_edna_fields_v3.png
│   ├── whale_edna_sim_v{1,2,3}.rds         Sim outputs (gitignored)
│   └── whale_edna_output_v{1,2,3}/         Stan fit artifacts (gitignored):
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
Rscript 00_pipeline_v3.r
```

…or run individual steps in order:

```
Rscript scripts/01_simulate_whale_edna_v3.r     # -> outputs/whale_edna_sim_v3.rds
Rscript scripts/02_plot_simulated_data_v3.r     # -> outputs/simulated_edna_fields_v3.png
Rscript scripts/03_format_stan_data_v3.r        # -> outputs/whale_edna_output_v3/stan_data.rds
Rscript scripts/04_run_whale_edna_model_v3.r    # -> outputs/whale_edna_output_v3/whale_edna_fit.rds
Rscript scripts/05_check_whale_edna_model_v3.r  # -> diagnostics in outputs/whale_edna_output_v3/
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
- **Plotting**: the 1 km grid in `02_plot_simulated_data_v3.r` evaluates
  the closed-form GP mean at every cell (501 × 1271 = 637k cells) — no
  kriging.

## Versioning convention

- Every file in `scripts/`, `stan/`, and the root `00_pipeline_v{N}.r`
  has an explicit `_v{N}` suffix (lowercase `v`).
- V{N} simulation writes `outputs/whale_edna_sim_v{N}.rds`.
- V{N} plotting reads `outputs/whale_edna_sim_v{N}.rds` and writes
  `outputs/simulated_edna_fields_v{N}.png`.
- V{N} Stan-data formatting writes
  `outputs/whale_edna_output_v{N}/stan_data.rds`.
- V{N} model runner writes
  `outputs/whale_edna_output_v{N}/whale_edna_fit.rds` and CmdStan files
  to the same directory.
- V{N} checking script writes all diagnostic plots and CSVs into
  `outputs/whale_edna_output_v{N}/`.

When starting a new iteration (`v4`), copy the highest-version scripts
(including `00_pipeline_v3.r`) to their `_v4` counterparts and update
this README with what changed and why.

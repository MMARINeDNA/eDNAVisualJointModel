# eDNAVisualJointModel

Joint modelling of visual survey and environmental DNA (eDNA) data for
marine megafauna off the US West Coast. The workflow simulates realistic
eDNA observations (qPCR + metabarcoding) for three species on a spatial
grid, and fits an anisotropic Hilbert-Space GP (HSGP) model to the
simulated data in Stan.

## Directory layout

```
.
├── README.md                       This file
├── LICENSE
├── eDNAVisualJointModel.Rproj
│
├── stan/                           Stan model source
│   ├── whale_edna_hsgp_V1.stan
│   └── whale_edna_hsgp_V2.stan
│
├── scripts/                        R pipeline scripts
│   ├── 01_simulate_whale_edna_V1.r     Simulate eDNA data (versions V1–V3)
│   ├── 01_simulate_whale_edna_V2.r
│   ├── 01_simulate_whale_edna_V3.r
│   ├── 02_plot_simulated_data_V2.r     Diagnostic plots of simulated data
│   ├── 02_plot_simulated_data_V3.r
│   ├── 03_run_whale_edna_model_V1.r    Compile + fit Stan + diagnostics
│   ├── 03_run_whale_edna_model_V2.r
│   ├── DetectionsBySpecies.R           Ad-hoc data exploration
│   └── FinWhales.R
│
├── outputs/                        Generated artifacts
│   ├── simulated_edna_fields_V2.png    Tracked plots of simulated data
│   ├── simulated_edna_fields_V3.png
│   ├── whale_edna_sim_V{1,2,3}.rds     Sim outputs (gitignored — regenerable)
│   └── whale_edna_output_V{1,2}/       Stan fit artifacts (gitignored)
│
├── Data/                           Real survey data (effort, sightings, detections)
└── Distance sampling trials/       Self-contained DS exploration (separate Stan model)
```

## Running the pipeline

All scripts expect the project root as the working directory. From the
repository root:

```r
# Simulate (writes outputs/whale_edna_sim_V3.rds)
Rscript scripts/01_simulate_whale_edna_V3.r

# Plot the simulated truth (writes outputs/simulated_edna_fields_V3.png)
Rscript scripts/02_plot_simulated_data_V3.r

# Compile + fit Stan model, write diagnostics to outputs/whale_edna_output_V2/
Rscript scripts/03_run_whale_edna_model_V2.r
```

The `.rds` and `whale_edna_output_*/` artifacts are `.gitignored` — they are
regenerable from the scripts.

## Versions

The simulation, plotting, runner, and Stan files are all versioned so that
each iteration of the model is reproducible. Older versions are kept for
comparison and because the Stan models target different simulation outputs.

### V1 — initial simulation + model

- **Domain**: Oregon / Washington coast only (UTM 10N, ~300 × 400 km).
- **Bathymetry**: logistic shelf-slope profile as a pure function of
  Easting `X` (isobaths run due N-S).
- **Observation design**: 150 stations × 5 sample depths (0, 50, 150, 300,
  500 m); metabarcoding marker MARVER3 at 5000 reads/replicate;
  `vol_filtered = 2.0 L`.
- **Species structure**: bathymetric habitat preference only (no latitude
  effect); shelf-slope-break specialist hake, shallow-shelf humpback,
  deep-slope PWSD.
- **Files**: `scripts/01_simulate_whale_edna_V1.r`,
  `scripts/03_run_whale_edna_model_V1.r`, `stan/whale_edna_hsgp_V1.stan`.

### V2 — code review refinements (with Ole)

- **Domain and bathymetry**: unchanged from V1.
- **Observation design refinements**:
  - Marker switched to MARVER1 (larger read depth — 50000 reads/replicate).
  - `vol_filtered` increased to 2.5 L; explicit aliquoting step
    (`vol_aliquot = 2 µL` drawn from a 100 µL elution).
  - Explicit animals → copies conversion (`conv_factor = 10`).
  - `zsample_pref` values shifted so baseline eDNA ≈ 1× rather than
    suppressing most depths.
- **Naming**: vectors carry an `_si` suffix (sample × species index order)
  for clarity downstream.
- **Stan model (V2)**: expanded observation model — qPCR Ct hurdle on
  hake, zero-inflated Beta-Binomial metabarcoding for all three species,
  `phi` parameterised in log total copies.
- **Plotting**: 3-row figure (surface density, Z_bathy vs λ, water-column
  effect).
- **Files**: `scripts/01_simulate_whale_edna_V2.r`,
  `scripts/02_plot_simulated_data_V2.r`,
  `scripts/03_run_whale_edna_model_V2.r`, `stan/whale_edna_hsgp_V2.stan`.

### V3 — extended domain, rotated bathymetry, realistic species distributions

- **Domain**: extended to cover the full US West Coast — San Francisco
  (37.77°N, ~4,180,000 m Northing) to the US/Canada border (~49°N,
  ~5,450,000 m Northing), 500 km × 1270 km in UTM Zone 10N.
- **Bathymetry**: no longer a pure function of `X`. A rotated cross-shore
  axis `X_prime` is used, applied in **normalised space** so the rotation
  angle has a consistent geometric meaning regardless of the tall domain
  aspect ratio. The rotation is **25°** — a pure 45° rotation confounds
  latitude with bathymetry in this elongated domain (at 45°, every
  northern station ends up "inshore" and every southern station
  "offshore", which breaks realistic species preferences).
- **Observation design**:
  - **300 stations** stratified 50/30/20 shelf / slope / offshore on the
    rotated coordinate (via rejection sampling in X/Y).
  - **3 sample depths** per station (0, 150, 500 m).
  - Samples with `Z_sample > Z_bathy` are dropped (can't sample below
    the seafloor).
- **Species spatial structure** (now realistic for the expanded range):
  - **Pacific hake** — shelf-slope break, broad in latitude.
  - **Humpback whale** — southern + inshore shelf (central CA to southern OR).
  - **Pacific white-sided dolphin** — northern + offshore slope
    (north of Cape Mendocino).
- **Model form**: bathymetric + latitude habitat structure is **absorbed
  into the GP** via a non-zero GP mean, so the exposed model is simply
  `log(λ) = μ + f`. Mathematically identical to a mean-function form
  (`μ + mean_fn + GP_residual`), but matches what the Stan model fits.
- **Plotting**: Row 1 evaluates the closed-form GP mean at every cell of
  a 1 km grid (501 × 1271 = 637k cells) — no kriging, just deterministic
  evaluation on a regular grid. Rows 2 and 3 unchanged in spirit.
- **Bug fix carried into V3**: `rmultinom` gets a `NaN` probability vector
  when all species have zero copies in an aliquot. This was latent in V2
  (rare) but common in V3 because strong habitat mismatch produces true
  zeros; fallback to uniform `pi_edna` for those rows.
- **Files**: `scripts/01_simulate_whale_edna_V3.r`,
  `scripts/02_plot_simulated_data_V3.r`. No V3 Stan runner yet — the V2
  Stan model is the intended target once ported to read V3 outputs.

## Versioning convention

- Every file in `scripts/` and `stan/` has an explicit `_V{N}` suffix.
- V{N} simulation writes `outputs/whale_edna_sim_V{N}.rds`.
- V{N} plotting reads `outputs/whale_edna_sim_V{N}.rds` and writes
  `outputs/simulated_edna_fields_V{N}.png`.
- V{N} model runner writes diagnostics to
  `outputs/whale_edna_output_V{N}/`.

When starting a new iteration, copy the highest-version files (e.g.
`*_V3.r` → `*_V4.r`) and update this README with what changed and why.

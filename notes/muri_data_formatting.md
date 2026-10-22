# MURI real-data formatting: assumptions, gaps, open questions

Code: `R/01_load_muri.r` (builds objects), `R/02_check_muri.r` (checks, summary table,
formatter smoke test). Outputs in `data/processed/`:

| file | contents |
|---|---|
| `muri_edna.rds` | all 551 bottles, `simulate_bathysp()` shape |
| `muri_visual.rds` | all 2779 effort segments, `simulate_visual()` shape |
| `etopo2022_60s_cce.rds` | cached bathymetry grid (not committed) |

No real-data analysis was run (only sanity fits, below). All three formatter calls run on the
objects (eDNA N=551; humpback/pwsd n_seg=2779). Nothing is excluded for lying outside a box:
the domain is built from the data (see Shared).

## Summary table

| | |
|---|---|
| eDNA bottles / stations | 551 / 182 |
| qPCR replicates (hake) | 2663 (detection rate 0.776, 2066 detected) |
| MB replicates | 1709 (reads in rep: hake 0.644, humpback 0.034, pwsd 0.048) |
| Visual segments / sections / effort | 2779 / 1015 / 12,967 km |
| Truncation `w` (95th percentile, per species) | humpback 6.45 km, pwsd 4.64 km |
| Sightings used within `w` (humpback / pwsd) | 548 / 89 (of 577 / 94 on effort with a distance) |
| Detection function | humpback: hazard-rate + group size; pwsd: half-normal + group size + Beaufort |

Sighting flow (all `sightings.csv` rows of the species -> on effort -> has pdist -> within `w`):
humpback 683 -> 646 -> 577 -> 548; pwsd 112 -> 98 -> 94 -> 89.

## Where each field comes from

**Visual**
- Segments: every row of `effort.csv` (all `ds_type = lneff`, cruise 1651 = 2018 CCES), ordered by
  start time. `X, Y` = projected `mlon/mlat` (segment midpoint). `seg_l` = `length_km`.
  `transect_id` = `section_id` (the file's continuous-effort section, 1015 groups, many with one
  segment). Alternatives (per day, per leg) are possible; nothing in the model uses it yet.
- Extra segment columns: `avgbft, avgvis, avgswellhght, mode, efftype, eswsides, speed_kph`,
  begin/end/mid times, `section_id, segnum, row_id, lat, lon, X_utm, Y_utm`.
- Truncation (decided): 5% of detections dropped per species, i.e. `w` = 95th percentile of the
  on-effort perpendicular distances: **humpback 6.454 km, pwsd 4.644 km** (`meta$w`, a named
  vector). `observed[[sp]]$obs` holds ALL on-effort sightings with a distance and group size, and
  `seg_count` is `tabulate(obs$seg_id)` of that full set. `format_stan_data_visual()` drops
  detections beyond `meta$w[[sp]]` and rebuilds `seg_count`. (The 5% rule applies to the two
  visual species; hake has no visual data.)
- Sightings: `sightings.csv` rows with `species_name` 76 (humpback) or 22 (pwsd). "On effort" =
  sighting time lies inside an effort segment's [begin, end] (`sightings.csv` has no on/off-effort
  flag; one sighting matched two segments, first taken). `seg_id` = that segment.
- **Perpendicular distance** = `pdist` in `data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv`, joined on
  sighting number x species code. I did not derive distance from lat/lon. `sightings.csv` and
  `effort.csv` contain no radial distance or bearing, so I could not re-check `pdist`;
  it is the survey's own perpendicular distance, **in km** (confirmed by the data owner; correct as is).
  SITEINFO lat/lon equal the sightings' to 1e-6 deg.
- `size` = `gsspbest` rounded to integer (identical to `group_size`; 2 NA rows are not among the
  used sightings). SITEINFO `Sppgs` has 50 fractional values (mean of observers); not used.
- Sighting covariates kept in `obs`: `bft, vis, swellhght, mixed, nsp, sightno, datetime, lat, lon, row_id`.

**eDNA**
- A bottle = `location_id` (identical across the qPCR and MB files; one lat/lon/depth each).
  `station` = distinct lat/lon (182). `sample_id` = row in `samples`.
- 11 bottles have two qPCR `sample` IDs (e.g. 225/226 at location 29). The MB file has no sample
  ID, so both were pooled into one bottle (IDs kept in `qpcr_sample_ids`). The two IDs can differ a
  lot (location 29: ~11 vs ~0.4 copies/uL), so they may be real duplicate bottles; needs a decision.
- qPCR: one entry per row of `hake_qPCR_MURI_df.csv`; `qpcr_detect = detected` (= copies non-NA).
- MB: `hake` = *Merluccius productus*, `humpback` = *Megaptera novaeangliae*, `pwsd` =
  *Lagenorhynchus obliquidens*; `junk` = replicate total (sum over all 295 taxa in
  `MV1_MURI_df.csv`) minus the three. Genus-level `Lagenorhynchus` (1302 reads in total, 1 replicate in
  ~1700 positive) is left in junk; so is *M. angustimanus* (46,845 reads). Total min 1, median
  34,184; 140 replicates have < 1000 reads; none have 0.
- Extras in `observed` (ignored by the formatter): `qpcr_copies_ul, qpcr_dilution,
  qpcr_inhibition, qpcr_plate, qpcr_sample_id, mb_pcr_replicate_id, mb_dilution, mb_rep`.

**Shared**
- UTM 10N (EPSG:32610). **Study area: north of 35°N** (`LAT_MIN`). Effort segments whose
  midpoint is south of 35°N, and their sightings, are dropped; every eDNA bottle is north of it.
  **Domain = bounding box of the kept effort and all eDNA bottles** (segment begin, middle and end
  points, and every bottle) plus a 10 km margin to the west, east and north, rounded to 10 km; the
  southern edge is 35°N at the central meridian (rounded down: N0 = 3870 km, ~34.96°N).
  `X = (E - E0)/1000`, `Y = (N - N0)/1000` with `(E0, N0) = (-10000, 3870000)` m (stored in
  `meta$utm_origin`), giving **X in [0, 690] km, Y in [0, 1770] km** (`meta$domain`, and
  `X_km_max` / `Y_km_max` in the eDNA `meta`). This box is HSGP4eDNA's `STUDY_DOMAIN`, which the
  simulators also use; the loader stops if the data give a different box. (Before the 35°N
  restriction the box was 1360 x 2600 km with N0 = 3040000, covering all effort down to 27.4°N.)
  `X_utm`, `Y_utm` (metres) are unchanged. The kept effort reaches about -121°W, so the UTM
  distortion far from zone 10's central meridian that the full effort had (to -114.5°W) is
  mostly gone.
- `Z_bathy` (m, positive down) for bottles and segments = **NOAA ETOPO 2022 60-arc-second bedrock**
  (downloaded with `marmap::getNOAA.bathy`), bilinear. Check against `sightings.csv` BATH field at
  2326 sightings: r = 0.998, median |diff| 12 m vs the file's ETOPO1 value, 19 m vs S2004.
  No NA; segment range 23-4909 m, bottle range 44-3308 m (the simulator draws >= 50 m).
- `zsample_effect` = 1; `bathy_spline = NULL`.

## Values I could not fill (all are placeholders or NA)

| field | status |
|---|---|
| `truth$qpcr_params` (kappa, alpha_ct, beta_ct, gamma_0, gamma_1, sigma_0) | **NA**; no standard curve in the data |
| `observed$qpcr_ct` | **all NA**; Ct not in the files. The formatter turns NA into 0, so detected replicates would be fit as Ct = 0: do not fit the qPCR part until Ct is supplied (`meta$ct_available = FALSE`) |
| `meta$vol_filtered` | 2.5 L, simulator placeholder; not in the data |
| `meta$vol_aliquot` | 2 uL, simulator placeholder; not in the data |
| `meta$vol_elution` | NA (model assumes 100 uL) |
| `meta$conv_factor` | simulator placeholder (as instructed) |
| eDNA `datetime`, `cruise` | NA; neither eDNA file has date, time or cruise |
| blanks / negatives | none present in either eDNA file, so no separate blanks object |

## Answers to the open questions

1. **qPCR.** `hake_qPCR_MURI_df.csv` has `hake_copies_ul` (= `conc`), `detected`, `dilution`,
   `inhibition_rate`, plate (`qPCR`, 30 plates H7-H39) and sample; no Ct anywhere in `data/`.
   Raw Ct and the standard curve (slope, intercept, efficiency, Ct noise) are needed from whoever ran
   the plates. Not back-calculated.
   - `dilution` is 1 for 1501 replicates, and 0.1/0.2/0.5 for 1162 (134 bottles). Within those
     bottles the same bottle is run at several dilutions (e.g. 0.1/0.2/0.5/1), so each bottle has
     up to 14 replicates. Whether copies/uL are already corrected for dilution is unknown; the
     model treats all replicates as undiluted aliquots.
   - `inhibition_rate`: median 0.10, range -1.58 to 1.50; 9.7% of replicates are > 0.5. I applied
     no correction or filter. Whether the diluted replicates, or only the least inhibited dilution, should
     be used is a decision for you.
   - 597 replicates are non-detects (copies NA).
2. **Filtered volume.** Not in any file, so I can't say if it varies. Placeholder only.
3. **Sampling depths** (`Z_sample`, m; bottles): 0 (172), 50 (152), 150 (74), 300 (73), 500 (65),
   and 14 singletons/pairs at 3, 40, 47, 48 (2), 280, 446, 464, 467, 475, 485-495. Two bottles have `Z_sample > Z_bathy`
   (ETOPO resolution; probably shelf stations); not dropped. `zsample_effect` = 1.
4. **Metabarcoding.** The files contain no marker/primer name, no ASV information and no
   assignment method; only taxon-level read counts after assignment. Taxa = 295 (mostly fish; also
   *Homo*). Replicate structure: 1709 PCR replicates over 551 bottles, 1-6 per bottle (3 in 417
   bottles). Each carries a template `Dilution` (1, 5, 10, 19, 20, 50, 100; fold convention, unlike
   the qPCR fraction) and `Rep` 1-6. 936 replicates are undiluted. Whether the `Nreads` are raw or
   filtered/rarefied is unknown; "junk" is relative to the 295 retained taxa only.
5. **Domain and dates.**
   - Domain: see Shared. All 551 bottles (X 256-522, Y 324-1496 km) are kept. Of the 2779 effort
     segments (12,967 km), 1641 (7,561 km) have midpoints north of 35°N and are kept (X 13-680,
     Y 2-1758 km); sightings on the dropped segments are counted in `sight_flow`
     (`n_south_of_lat_min`). The bottles span about 39% of the domain's width and 66% of its
     height. **The processed objects in data/processed/ still hold the pre-restriction data until
     the loader is re-run** (see ROADMAP_2026-10-06, Phase 6).
   - Dates: effort 2018-06-28 to 2018-12-03 (2018 CCES, cruise 1651); sightings 2018-07-01 to
     2018-10-16. eDNA: no dates. Not matched, as agreed (different cruises).

## Other flags

- **Truncation.** Perpendicular distance, km (SITEINFO, all rows): humpback n=577, quantiles
  50/80/90/95/99/max = 2.15/4.5/5.31/6.45/8.4/9.77; pwsd n=94 = 1.28/2.86/3.88/4.64/5.41/6.41.
  Decision: drop the largest 5% per species (w = 6.454 / 4.644 km; 29 and 5 detections dropped).
  This follows Len Thomas's CalCOFI analysis (https://rpubs.com/lenthomas/1268708), which truncates
  per species at roughly the upper 5-10% of distances.
- **Sightings without distance.** 69 humpback and 4 pwsd on-effort sightings have no `pdist`
  (not in SITEINFO; they span all Beaufort states and modes, so SITEINFO's filter is unknown).
  They are excluded, so counts may be biased low. 37 humpback and 14 pwsd sightings fall
  outside every effort segment (off effort) and are excluded.
- **Mixed-species sightings** are one row per species, sharing `sightno`, with `mixed = TRUE`, `nsp`,
  and a species-specific `gsspbest`. Used rows: 6 humpback, 16 pwsd.
- **One-sided effort.** 122 segments have `eswsides = 1` (one side observed); with a two-sided strip
  `seg_l` over-states their effort. Not adjusted; the column is kept.
- **Effort mode.** `mode` C (2436) / P (343) and `efftype` F (1114) / N (1665) are kept, not
  filtered. `avgbft` goes up to 6 (Beaufort not filtered).
- **Group-size pools** (`data/grpsz/*.rds`): humpback pool n=683 (all humpback sightings, including
  those not used), mean 1.72, median 1, max 8; real used sightings mean 1.69, median 1, max 8: consistent.
  pwsd pool n=110, mean 42.8, median 21, q90 81, max 452; real used mean 46.2, median 24,
  q90 107, max 452: same shape, real is slightly larger.
- **Submodule.** `external/HSGP4eDNA` is not checked out in this worktree; I tested against a
  clone at the pinned commit e8d0432 outside the repo and did not touch the submodule.
- `data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv` and `data/MV1_MURI_df.csv` are untracked in the main
  checkout; I copied them into this worktree's `data/` for the run and did not commit them.

## Detection function

Selection followed Len Thomas's CalCOFI recipe (forward selection by AIC with the Distance
package, `nadj = 0`, candidates group size, Beaufort, visibility, swell, effort mode), run on the
on-effort sightings truncated at the 95th percentile, using **segment-level** `avgbft`, `avgvis`,
`avgswellhght`, `mode`:

| species | key | steps (AIC change) | AIC selection | used |
|---|---|---|---|---|
| humpback (n=548) | hazard-rate (~20 AIC better than half-normal at w = 5 km) | size -3.2 | group size | group size |
| pwsd (n=89) | half-normal | visibility -6.8, Beaufort -1.9, size -1.2 | visibility + Beaufort + group size | **group size + Beaufort** (visibility dropped by choice) |

Beaufort and group size for pwsd are marginal (< 2 AIC each); kept because the rule was "follow
the AIC". In the sanity fit visibility was also weak (+0.10 per SD, 90% interval -0.02 to 0.21).
Swell and effort mode were never selected. There is no ship / observer-height covariate (one DAS
file, `CCES1651.das`). These are exploratory fits (small n for pwsd, modes and conditions pooled).
Per-species settings live in `truth$sp_params`: `det_key` (humpback 1, pwsd 0) and `use_bft`
(humpback FALSE, pwsd TRUE).

Model support (`stan/hsgp_visual.stan`, after merging main's Phase 4b):
- **Beaufort** (main's mechanism): segments carry `bft` (here the segment mean `avgbft`), binned
  to 0.5 steps; `log sigma = log sigma0 + beta_size (s - s_c) + beta_bft (bft - 3)`; the
  size-averaged ESW is computed once per level and used for both detections and segments.
- **Detection key** (`det_key`): 0 half-normal, 1 hazard-rate `g(x) = 1 - exp(-(sigma/x)^b)`,
  ESW by a 400-cell midpoint rule, `log b ~ normal(log 3, 0.7)`.
- **Binned group-size sum:** the ESW sum over group sizes uses sizes 1-30 individually then 10%
  log bins (a 677-detection pwsd test fit went from 15 to 2 min, ESW error 4e-5).
- With `det_key = 0` the model is main's model.

### Validation

- `tests/visual_detection.R`: parameter recovery (99% intervals) and per-level ESW vs R's
  `integrate()` for hazard-rate + size and half-normal + size + Beaufort on simulated data.
- `tests/smoke.R` (includes main's visual_bft case).
- **Real-data sanity fits** (before the merge, non-spatial, 2 chains x 400, no divergences; NOT a
  final analysis): humpback mean detection probability 0.44 (90% 0.36-0.52), consistent with
  `Distance`'s hazard-rate fit (0.43 at w = 6.5 km); hazard-rate shape b = 0.98 (90% 0.78-1.23),
  i.e. a heavy tail despite a prior centred on 3. pwsd Beaufort effect about -0.55 per SD.

### Caveats

- cmdstanr does not recompile when only `stan/include/*.stan` changes; pass `force_recompile = TRUE`.
- Beaufort is the segment mean, not the value at each sighting (the model needs it for every segment).
- g(0) = 1 is still assumed (no availability or perception bias correction).
- The loader asserts that every segment and bottle is inside the domain; if data are added, the
  domain is recomputed from them.

## Correction: group size for humpback (2026-10-12)

The AIC selection (table above) chose group size for humpback, and this note said it was
used. It was not: `R/01_load_muri.r` set `det_key` and `use_bft` for humpback but left
`use_size_covar` at the simulator's value (off), so the humpback fits ran without a size
effect on detection. The loader now sets `use_size_covar = 1L` for humpback and the processed
objects were regenerated (on the 35°N domain). On the current set-up, switching it on changes
humpback's mean group size from 1.70 to 1.69 and mean density by -0.8%; `beta_size` is 0.029
[-0.043, 0.101], so the data say little about it. PWSD was correct throughout (size effect on:
beta_size 0.004, population mean group size 37.0 against a detected mean of 46.6, matching the
37.4 mean of groups seen within 1 km). Cutting PWSD's group-size support from 1000 to 500 moves
its mean group size by -1.3% and density by -0.5%.

## eDNA from the lab files (2026-10-21)

This supersedes the eDNA parts above (the summary table's eDNA rows, "Where each field comes
from: eDNA" and answers 1-4), which describe the earlier `hake_qPCR_MURI_df.csv` /
`MV1_MURI_df.csv` objects. The visual object is unchanged.

**Inputs** (private, git-ignored; `RAW_DIR` points the loader at them):
`qPCR_2019A.rds` (hake qPCR, analysis A: wells with Ct, standards, controls),
`M3_compiled_taxon_table_wide.csv` (reads, taxa x replicate columns, six markers),
`Hake_2019_metadata.csv` (per sample). `R/fit_qpcr_standards.r` must run before
`R/01_load_muri.r`; it writes `data/processed/qpcr_standard_curves.rds`.

**Choices** (Ole's answers are in `docs/dataQs.md`):
- **Bottles** = MarVer1 samples that are not controls or field negatives, north of 35°N; matched
  across files on sample number (qPCR `52193_282` = metadata / M3 `52193-282`). 569 bottles;
  56 flagged `30EtOH` / `30EtOHpaired` are dropped (`ETOH=keep` keeps them, flagged in
  `etoh_wash`); 513 kept, at 167 stations. 16 station x depth locations have two bottles.
- **Depth**: `Z_sample` = `depth_cat` from analysis A (target depth: 0 169, 50 146, 150 65,
  300 67, 500 66 bottles), also as an ordered factor `depth_cat` and `depth_idx`; recorded depth
  kept as `depth_recorded` (recorded minus target: -36 to +3 m). `location` = station x
  `depth_cat`.
- **Bottom depth** (`Z_bathy`): the metadata's `bottom.depth.consensus`, as the lab recommends;
  ETOPO 2022 at the bottle is kept as `Z_bathy_etopo`. Consensus minus ETOPO: median -7.6 m,
  median absolute difference 13 m, 90% within 46 m, extremes -150 and +320 m (steep slope
  stations). Visual segments still use ETOPO (there is no consensus value off the stations).
  5 bottles have `Z_sample > Z_bathy`: all 500 m targets at stations 474-499 m deep, recorded
  at 464-491 m, i.e. near-bottom samples.
- **qPCR**: all analysis-A hake wells of these bottles, as provided (cleaned by the lab: one
  dilution per sample, no inhibited wells, late Ct kept). 1,532 wells on 25 plates; 3 per bottle
  for 458 bottles; detection rate 0.858; Ct 26.7-44.7; dilution 1 (1,139 wells) or 0.2 (393).
  `qpcr_ct` is NA for non-detects.
- **Standard curves**, one per plate: `stan/qpcr_standards.stan` (the vignette's multi-plate model
  and priors: per-plate beta0, beta1 <= -1/ln 2; shared gamma0, gamma1, sigma0, kappa;
  lambda = 2 uL x copies/uL) on the 659 hake standards with inhibition < 0.75 (106 dropped, all
  10^5 standards among them), 33 plates, 26 non-detects. 4 x 1000 draws, no divergences,
  max R-hat 1.003. alpha_ct 36.1-39.9, beta_ct 1.448-1.493 (efficiency 0.95-1.00), gamma0 0.26,
  gamma1 -0.37, sigma0 0.14, kappa 0.90. Checks (`outputs/qpcr/checks.txt`): 97.0% of detected
  standards inside their 95% predictive interval; detection at 1 copy/uL 0.865 observed vs 0.836
  fitted. Each well carries its plate's posterior-mean curve (`qpcr_alpha_ct`, `qpcr_beta_ct`,
  `qpcr_plate_idx`; `truth$qpcr_curves`); `truth$qpcr_params` holds the plate means
  (alpha 37.6, beta 1.46) and the shared terms, for the current single-curve formatter.
- **Metabarcoding**: MarVer1 columns of the M3 table (column names give sample, dilution fold,
  replicate and run). Replicates with no reads after contaminant removal are dropped; where a
  bottle was run at several folds (52 bottles), the fold with the most reads is kept (Ole's
  rule); a duplicated sample x fold x replicate keeps the larger. 1,438 replicates (3 for 402
  bottles); folds kept 1 (787), 5 (195), 10 (419), 19-100 (37). 7 bottles have no replicate with
  reads and keep only their qPCR (`meta$flow$bottles_without_mb`). 136 replicates have < 1000 reads.
- **Dilution offsets**: `qpcr_log_dilution` = log(fraction); `mb_log_dilution` = -log(fold).
  Neither is used by the HSGP4eDNA v1.5 likelihood yet.
- **Contamination**: *Gallus* (127,162 of 98.4 M MarVer1 reads) is removed from both the junk
  column and the replicate total (`mb_contaminant_reads`); the list in the loader also has human,
  cattle and pig, none present. Junk = total - contaminants - the three targets.
- **Volumes, dates**: `vol_filtered` per bottle (2.5 L for 511; one 2.0, one 2.4); `meta$vol_filtered`
  is the median, for the formatter. `date` 2019-07-02 to 2019-08-17. Template volume 2 uL for
  qPCR (the vignette); the MarVer1 template volume is assumed the same (asked). Elution volume
  100 uL, assumed: common to both assays, it only rescales `conv_factor`. `cruise` = "52193"
  (project ID; one survey). `conv_factor` still the simulator placeholder.

**Checks**: `R/02_check_muri.r` passes (formatter N 513, N_qpcr 1532, N_mb 1438; all bottles in
the domain). Bottle-level agreement, hake qPCR detected vs metabarcoding reads present:
both 392, qPCR only 82, metabarcoding only 5, neither 27.

**Lab `copies_ul`** (kept as `qpcr_copies_ul`; not used by the model, which works from Ct). They
are an exact per-plate log-linear function of Ct with shallower slopes than ours (3.13-3.39 Ct per
log10 copy, about 105% efficiency, against our 3.33-3.41 with the doubling bound), as the
vignette's tangent on machine-software curves would predict. So they run 11-40% below copies
from our curves, the gap growing with Ct, and they are per uL of the template as run (not
corrected for dilution: the 1:5 wells are not 5x higher). Not a question for the lab unless we
start using them.

**Negatives**: hake is detected in 73 of 114 evaporated-DI and 25 of 72 lab-DI field-negative
wells (up to 12 and 10 copies/uL), 1 of 174 extraction controls and 0 of 3 purification
negatives. Not used yet; a question for the lab.

**Detection function, 35°N re-check**: the forward AIC selection on the kept effort (548 humpback,
89 PWSD detections within w) reproduces the table above: humpback hazard-rate + group size
(hazard-rate 17 AIC better than half-normal); PWSD half-normal + visibility (-6.8) + Beaufort
(-1.9) + group size (-1.2). The PWSD set in use (group size + Beaufort, no visibility) is 9.7 AIC
worse than that and 0.2 better than no covariates. Decision (2026-10-21): add visibility to the visual
model and use the AIC set for PWSD (ROADMAP_2026-10-21, section 5 item 5).

## Visual effort and covariates (2026-10-22)

- **Beaufort ≤ 5.** Segments with mean Beaufort above 5 are dropped (`BFT_MAX`
  in `R/01_load_muri.r`): 24 of 1,641 segments (120 of 7,561 km) and 1
  humpback sighting (`sight_flow$n_bft_above_max`); 1,617 segments and
  7,441 km remain. The domain is built from all effort north of 35°N
  before this filter, so it and the eDNA object are unchanged. Humpback's
  truncation distance moves from 6.454 to 6.461 km (547 detections within
  it); PWSD is unchanged (4.644 km, 89). Effort by Beaufort (rounded segment
  mean, before the cut): 0: 146 km, 1: 622, 2: 886, 3: 1,749, 4: 2,174,
  5: 1,896, 6: 88.
- **Visibility** (`avgvis`, segment mean; `vis` in the segments): the
  observers' estimate of the distance in **nautical miles** at which a
  dolphin could be seen surfacing against the water (SWFSC DAS format,
  "DasAll Events and Associated Data Fields", updated 2020 by J. Carretta,
  weather event "W" field 5; bundled with the swfscDAS R package as
  `inst/DAS_Format.pdf`). The OBIS-SEAMAP metadata only call it "a measure
  of air/atmospheric clarity", range 0.5–7.5. Recorded in about 0.5 nmi
  steps; segments range 1–7.5 nmi (1.9–13.9 km), median 6 nmi. The model
  treats it as a continuous covariate on log σ (per nmi). Visibility below
  the truncation distance (where it could cut detection off, not only scale
  it) is limited: 660 km of the 7,441 km for humpback (below 6.46 km = 3.5 nmi)
  and 59 km for PWSD (below 4.64 km = 2.5 nmi).


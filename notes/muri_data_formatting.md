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
| Detection function | humpback: hazard-rate + group size; pwsd: half-normal + group size + visibility + Beaufort |

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
- UTM 10N (EPSG:32610). **Domain = bounding box of all effort and eDNA bottles** (segment begin,
  middle and end points, and every bottle) plus a 10 km margin, rounded to 10 km:
  `X = (E - E0)/1000`, `Y = (N - N0)/1000` with `(E0, N0) = (-10000, 3040000)` m (stored in
  `meta$utm_origin`), giving **X in [0, 1360] km, Y in [0, 2600] km** (`meta$domain`, and
  `X_km_max` / `Y_km_max` in the eDNA `meta`). The origin moved from the earlier (100000, 4180000)
  so that the domain still starts at 0, which the downstream code assumes; `X_utm`, `Y_utm`
  (metres) are unchanged. Margins to the edges are 12-16 km. Because the effort runs to
  about -114.5 W, parts of it are 6-8 degrees from zone 10's central meridian, where UTM
  distances are distorted by up to a few percent (not corrected).
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
   - Domain: see Shared. All 551 bottles (X 256-522, Y 1154-2326 km) and all 2779 segments
     (X 16-1344, Y 12-2586 km) are inside, so nothing is dropped. The domain is 1360 x 2600 km, so
     the same HSGP basis budget (17 x 23 = 391 functions) covers the area more coarsely than
     before, and most of the box has no eDNA data (bottles occupy about 20% of its width).
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

## Detection function (added after the first version)

Selection followed Len Thomas's CalCOFI recipe (forward selection by AIC with the Distance
package, `nadj = 0`, candidates group size, Beaufort, visibility, swell, effort mode), run on the
on-effort sightings truncated at the 95th percentile, using **segment-level** `avgbft`, `avgvis`,
`avgswellhght`, `mode` (the values the Stan model has for every segment):

| species | key | steps (AIC change) | selected |
|---|---|---|---|
| humpback (n=548) | hazard-rate (chosen; ~20 AIC better than half-normal at w = 5 km) | size -3.2 | group size |
| pwsd (n=89) | half-normal (chosen) | visibility -6.8, Beaufort -1.9, size -1.2 | visibility + Beaufort + group size |

Beaufort and group size for pwsd are marginal (< 2 AIC each); kept because the rule was "follow
the AIC". Swell and effort mode were never selected. There is no ship / observer-height covariate
(one DAS file, `CCES1651.das`). These are exploratory fits (small n for pwsd, mode and
conditions pooled); per-species `det_key` and `det_covars` live in `truth$sp_params`.

Model support: `stan/hsgp_visual.stan` has `det_key` (0 half-normal, 1 hazard-rate
`g(x) = 1 - exp(-(sigma/x)^b)`, ESW by a 400-cell midpoint rule) and detection covariates
`Z_det` (centred, scaled; segment-level, rounded to Beaufort 1 and visibility 0.5 so the number of
unique covariate rows stays small). The ESW, and so the encounter rate, is computed per covariate
row. With `det_key = 0` and no covariates it is the original model.

### Validation of the new detection options

- **Same model when switched off.** With `det_key = 0` and no covariates the log-posterior mode
  matches the original model (humpback `lp__` identical to 7 digits; pwsd within 0.01, a flat surface).
- **Recovery on simulated data** (`tests/visual_detection.R`, non-spatial): hazard-rate + size
  (log sigma 0.98 vs 0.92 true, b 3.3 vs 3, size effect 0.08 vs 0.10) and half-normal + size +
  two covariates (log sigma 0.402 vs 0.405, covariates 0.31 / -0.34 vs 0.3 / -0.3). Stan's
  covariate-row ESW agrees with R's `integrate()` to 1e-7 (hazard-rate) and 4e-5 (binned size sum).
- **Speed.** The group-size sum for the ESW is binned (sizes 1-30 individually, then 10% log bins),
  which cut a 677-detection pwsd fit from 15 to 2 min with unchanged results.
- **Real-data sanity fits** (non-spatial, 2 chains x 400, no divergences; NOT a final analysis):
  - humpback: mean detection probability 0.44 (90% 0.36-0.52), consistent with `Distance`'s
    hazard-rate fit (0.43 at w = 6.5 km). Hazard-rate shape b = 0.98 (90% 0.78-1.23), i.e. a
    heavy tail, despite a prior centred on 3.
  - pwsd: mean detection probability 0.49. Beaufort effect -0.55 per SD (90% -0.80 to -0.34);
    visibility +0.10 per SD (90% -0.02 to 0.21), much weaker than its AIC ranking suggested.

### Caveats

- cmdstanr does not recompile when only `stan/include/*.stan` changes; pass `force_recompile = TRUE`.
- Covariates are segment-level (not sighting-level) and rounded (Beaufort to 1, visibility to 0.5),
  and centred / scaled by effort-weighted mean and SD over the segments in the object passed to the
  formatter; `log_sigma` is the scale at the mean covariates. Subsetting segments changes that
  reference point.
- g(0) = 1 is still assumed (no availability or perception bias correction).
- The loader asserts that every segment and bottle is inside the domain; if data are added,
  the domain is recomputed from them.
- The covariate centring is computed over all segments in the object passed to the formatter.

# MURI real-data formatting: assumptions, gaps, open questions

Code: `R/01_load_muri.r` (builds objects), `R/02_check_muri.r` (checks, summary table,
formatter smoke test). Outputs in `data/processed/`:

| file | contents |
|---|---|
| `muri_edna.rds` | all 551 bottles, `simulate_bathysp()` shape |
| `muri_visual.rds` | all 2779 effort segments, `simulate_visual()` shape |
| `muri_edna_indomain.rds`, `muri_visual_indomain.rds` | same, restricted to the model box (visual only changes; `seg_id` renumbered) |
| `etopo2022_60s_cce.rds` | cached bathymetry grid (not committed) |

No models were fitted. All three formatter calls run on the full objects
(eDNA N=551; humpback/pwsd n_seg=2779) and on the in-domain subsets.

## Summary table

| | |
|---|---|
| eDNA bottles / stations | 551 / 182 |
| qPCR replicates (hake) | 2663 (detection rate 0.776, 2066 detected) |
| MB replicates | 1709 (reads in rep: hake 0.644, humpback 0.034, pwsd 0.048) |
| Visual segments / sections / effort | 2779 / 1015 / 12,967 km |
| Sightings used (humpback / pwsd) | 500 / 92 (w = 5 km) |
| In-domain subset | 1132 segments, 5,273 km, 269 humpback + 51 pwsd sightings |

Sighting flow (all `sightings.csv` rows of the species -> on effort -> has pdist -> <= 5 km):
humpback 683 -> 646 -> 577 -> 500; pwsd 112 -> 98 -> 94 -> 92.

## Where each field comes from

**Visual**
- Segments: every row of `effort.csv` (all `ds_type = lneff`, cruise 1651 = 2018 CCES), ordered by
  start time. `X, Y` = projected `mlon/mlat` (segment midpoint). `seg_l` = `length_km`.
  `transect_id` = `section_id` (the file's continuous-effort section, 1015 groups, many with one
  segment). Alternatives (per day, per leg) are possible; nothing in the model uses it yet.
- Extra segment columns: `avgbft, avgvis, avgswellhght, mode, efftype, eswsides, speed_kph`,
  begin/end/mid times, `section_id, segnum, row_id, lat, lon, X_utm, Y_utm, in_domain`.
- Sightings: `sightings.csv` rows with `species_name` 76 (humpback) or 22 (pwsd). "On effort" =
  sighting time lies inside an effort segment's [begin, end] (`sightings.csv` has no on/off-effort
  flag; one sighting matched two segments, first taken). `seg_id` = that segment.
- **Perpendicular distance** = `pdist` in `data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv`, joined on
  sighting number x species code. I did not derive distance from lat/lon. `sightings.csv` and
  `effort.csv` contain no radial distance or bearing, so I could not re-check `pdist`;
  I assume it is the survey's own perpendicular distance, **in km**. SITEINFO lat/lon equal the
  sightings' to 1e-6 deg. **Please confirm provenance and units.**
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
- UTM 10N (EPSG:32610), `X = (E-100000)/1000`, `Y = (N-4180000)/1000`.
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
5. **Domain box and dates.**
   - eDNA: all 551 bottles inside (X 146-412, Y 14-1186).
   - Visual: **1647 of 2779 segments (7,694 of 12,967 km) are outside the box**: 148 west
     (X < 0, includes the ~128.8 W effort), 1144 east (X > 500), 325 south (Y < 0, south of
     ~37.7 N), 30 north (Y > 1270). 231 of 500 humpback and 41 of 92 pwsd sightings are on
     out-of-box segments. With the full object, `format_stan_data_visual()` runs but
     `coords` span -2.8..3.9, outside [-1, 1] (HSGP is not valid there), so use the in-domain
     subset or move the domain.
   - Dates: effort 2018-06-28 to 2018-12-03 (2018 CCES, cruise 1651); sightings 2018-07-01 to
     2018-10-16. eDNA: no dates. Not matched, as agreed (different cruises).

## Other flags

- **Truncation.** Perpendicular distance, km (SITEINFO, all rows): humpback n=577, quantiles
  50/80/90/95/99/max = 2.15/4.5/5.31/6.45/8.4/9.77, 13% > 5 km; pwsd n=94 = 1.28/2.86/3.88/4.64/5.41/6.41,
  2% > 5 km. w = 5 km (kept) drops 77 humpback + 2 pwsd; the 95th percentile suggests ~6.5 km for
  humpback and ~4.6 km for pwsd (a single w could be 6.5).
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

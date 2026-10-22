# Real visual fits on the 35°N data: NB counts, visibility, basis budget (2026-10-22)

`R/fit_visual_real.r` on `data/processed/muri_visual.rds` (1,641 segments,
7,561 km), with the stored settings: humpback hazard-rate + group size + NB
counts; PWSD half-normal + group size + Beaufort + visibility + NB counts.
4 chains × 1000/1000 unless noted. Outputs in `outputs/visual/real/`
(git-ignored); the budget-800 runs without visibility and the Poisson runs
are from `R/check_count_dispersion.r` (4 × 500/500,
`outputs/visual/count_dispersion/`).

## Humpback (548 detections)

| Run | Basis (functions) | Floor (km) | Length-scales (km, median) | Floor flags | GP SD | Mean density (/km², 90%) | φ (90%) |
|---|---|---|---|---|---|---|---|
| Poisson, budget 800 (500/500) | 792 | 42 × 65 | 33 × 47 | both | 2.23 | 0.0226 [0.019, 0.027] | |
| NB, budget 800 (500/500) | 792 | 42 × 65 | 35 × 50 | both | 2.12 | 0.0229 [0.019, 0.028] | 2.9 [1.8, 5.6] |
| NB, budget 1600 | 1581 | 30 × 46 | 23 × 32 | both | 1.91 | 0.0230 [0.019, 0.028] | 8.3 [3.2, 136] |

Budget 1600: 68 min, no divergences, max R-hat 1.003, E-BFMI ≥ 0.88.
Detection: hazard-rate shape b = 0.97, p_det = 0.44, mean group size 1.69,
β_size 0.03 [−0.05, 0.11].

## PWSD (89 detections)

| Run | Basis | Floor (km) | Length-scales (km) | Flags | Mean density (/km², 90%) | β_vis (90%) | β_bft (90%) | φ (90%) |
|---|---|---|---|---|---|---|---|---|
| Poisson, no visibility, 800 (500/500) | 792 | 42 × 65 | 42 × 72 | both | 0.094 [0.069, 0.128] | | | |
| NB, no visibility, 800 (500/500) | 792 | 42 × 65 | 46 × 77 | both | 0.098 [0.070, 0.137] | | | 1.1 [0.5, 8.6] |
| NB, no visibility, 800 | 792 | 42 × 65 | 46 × 76 | both | 0.099 [0.071, 0.136] | | −0.10 [−0.26, 0.03] | 1.1 [0.5, 5.4] |
| NB + visibility, 800 | 792 | 42 × 65 | 46 × 78 | both | 0.107 [0.076, 0.153] | 0.17 [0.08, 0.27] | −0.16 [−0.29, −0.04] | 1.0 [0.45, 4.4] |
| NB + visibility, 1600 | 1581 | 30 × 46 | 29 × 69 | both | 0.106 [0.076, 0.149] | 0.17 [0.08, 0.26] | −0.16 [−0.29, −0.04] | 1.7 [0.6, 27] |

Budget 800 with visibility: 138 min (105 detection-condition cells against
13 Beaufort levels), 1 divergence, max R-hat 1.002; without visibility (same
settings otherwise, 4 × 1000/1000): 47 min, 1 divergence, max R-hat 1.002. Budget 1600: 139 min, no
divergences, max R-hat 1.003. Mean group size 37.7, p_det 0.49.

## What this shows

- **Visibility matters for PWSD in the spatial model too.** β_vis is
  clearly positive (0.17 per visibility unit on log σ, 90% interval
  0.08–0.27, against 0.21 in the distance-only AIC fit), and with it in,
  Beaufort's effect is clearer (−0.16 against −0.10 [−0.26, 0.03] without
  it). Mean density rises by about 8% (0.099 → 0.107), and where it rises
  follows visibility: the ratio of the two surfaces (with ÷ without) is 1.09
  in the southern third of the domain (mean segment visibility 5.1) and 1.03
  in the middle and northern thirds (5.9, 6.0); its correlation with segment
  visibility is −0.62. The two surfaces still correlate at 0.999, so the
  effect is a modest south-north tilt, not a different map. The cost is
  run time: 105 detection-condition cells instead of 13 make the fit about
  3 times slower (138 against 47 min).
- **The length-scales follow the basis floor.** Doubling the budget lowers
  the floor from 42 × 65 to 30 × 46 km, and humpback's posterior follows it
  on both axes (35 × 50 → 23 × 32 km), as on the old domain. PWSD's
  east-west length-scale follows it too (46 → 29 km); its north-south one
  settles near 70 km, above the floor.
- **Density does not depend on it** (humpback 0.0229 → 0.0230; PWSD 0.107 →
  0.106), nor do the detection parameters.
- **Overdispersion and fine-scale field trade off.** With the larger basis,
  humpback φ rises from 2.9 to 8.3 (the field takes up part of what looked
  like clustering at budget 800), confirming that segment counts alone
  cannot separate the two.

## Reading

The humpback field (and PWSD's along the east-west axis) varies on scales at
or below about 25–30 km, five or six 5 km segments, finer than the basis can
represent at either budget. The floor flag is then a statement about the
basis, not a fault in the fit: density, detection and group size are
stable. A finer basis would cost roughly in proportion (humpback 68 min at
1600) and the length-scales may keep following it.

**Decision (2026-10-22):** keep the 800 budget; the floor flag is reported as
a resolution limit of the basis (ROADMAP_2026-10-21, section 6).

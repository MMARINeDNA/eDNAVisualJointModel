# Clustered sightings and the length-scale floor (2026-10-21)

ROADMAP_2026-10-21 section 3 item 1. Question: the real visual fits put the GP
length-scales at the prior floor (explainer, "Refit at twice the basis
budget"). Is that caused by clustered sightings, i.e. segment counts more
variable than Poisson, and does a negative-binomial term remove it?

Code: `R/check_count_dispersion.r` (simulate on the real effort, fit, floor
check), the simulator options in `R/functions_visual.R` (`nb_size`,
`segments`), and the `use_nb` option in `stan/hsgp_visual.stan` and
`stan/hsgp_joint.stan`. Outputs: `outputs/visual/count_dispersion/`
(git-ignored). Seed 1 only; 4 chains x 500/500; basis budget 800.

## Set-up

- **Design:** the real 35°N effort (1,641 segments with their lengths and
  bottom depths, `data/processed/muri_visual.rds`). Its length-scale floor at
  the 800-function budget is **42 km (east-west) x 65 km (north-south)**.
  The explainer's 75 x 106 km floor was for the earlier fits on all 2,779
  segments, before the 35°N restriction.
- **Truth:** the simulator's field, with length-scales 150 x 250 km (well
  above the floor), so a fit that follows the floor is wrong. PWSD's mean
  density was raised (log 0.07) to give about the real number of detections.
- **Clustering:** each segment's expected groups multiplied by a
  gamma(k, k) draw, independently between segments, so groups per segment are
  negative binomial with size k (k = Inf: Poisson). Detection thinning keeps
  it negative binomial with the same k.
- **Calibration** (mgcv GAMs of real group counts on a 2-D smooth and bottom
  depth): humpback's overdispersion depends on how smooth the field is
  allowed to be. With about 38 effective df, NB θ = 1.1 and NB is preferred
  by 16 AIC; with about 71 df, the Poisson fits (Pearson dispersion 0.6).
  PWSD (89 groups): θ about 1, no AIC preference. Raw variance/mean of
  segment counts: humpback 2.5, PWSD 1.7.

## Results

Length-scales are posterior medians (km, east-west x north-south); "flags" is
`floor_check()` on both axes; R² is the field at the segments against the
truth.

| Species | Sightings | var/mean | Model | Length-scales | Flags | R² | φ (median, 95%) |
|---|---|---|---|---|---|---|---|
| humpback | Poisson | 1.20 | Poisson | 171 x 207 | 0 | 0.89 | |
| humpback | Poisson | 1.20 | NB | 173 x 212 | 0 | 0.89 | ~100 [10, 3e4] |
| humpback | NB k = 1 | 1.60 | Poisson | 130 x 254 | 0 | 0.83 | |
| humpback | NB k = 1 | 1.60 | NB | 133 x 274 | 0 | 0.82 | 1.01 [0.72, 1.59] |
| humpback | NB k = 0.3 | 3.35 | Poisson | **48 x 92** | **2** | **0.65** | |
| humpback | NB k = 0.3 | 3.35 | NB | 131 x 163 | 0 | 0.83 | 0.31 [0.24, 0.40] |
| PWSD | Poisson | 1.44 | Poisson | 169 x 276 | 0 | 0.85 | |
| PWSD | Poisson | 1.44 | NB | 167 x 275 | 0 | 0.85 | ~21 [2.3, 6e3] |
| PWSD | NB k = 1 | 1.71 | Poisson | 165 x 214 | 0 | 0.70 | |
| PWSD | NB k = 1 | 1.71 | NB | 154 x 212 | 0 | 0.70 | 0.78 [0.36, 2.34] |
| PWSD | NB k = 0.3 | 2.34 | Poisson | **80 x 139** | **2** | 0.78 | |
| PWSD | NB k = 0.3 | 2.34 | NB | 115 x 267 | 0 | 0.79 | 0.30 [0.17, 0.60] |

Real 35°N data:

| Species | var/mean | Model | Length-scales | Flags | GP SD | Mean density (/km²) | φ (median, 90%) |
|---|---|---|---|---|---|---|---|
| humpback | 2.5 | Poisson | 33 x 47 | 2 | 2.23 | 0.0226 | |
| humpback | 2.5 | NB | 35 x 50 | 2 | 2.12 | 0.0229 | 2.9 [1.8, 5.6] |
| PWSD | 1.7 | Poisson | 42 x 72 | 2 | 2.20 | 0.094 | |
| PWSD | 1.7 | NB | 46 x 77 | 2 | 2.17 | 0.098 | 1.1 [0.5, 8.6] |

No fit had more than 9 divergences (PWSD, NB model on Poisson data; 1 for the
Poisson model on the same data); max R-hat ≤ 1.015.

## What this shows

1. **Strong clustering does produce floor-following**, and the NB term
   removes it. With k = 0.3 (variance/mean 3.4 for humpback, above the real
   2.5) the Poisson model's length-scales fall to the floor and its field
   accuracy drops (humpback R² 0.65 against 0.83 with the NB model); the NB
   model recovers the length-scales and k.
2. **Moderate clustering (k = 1) does not**, with either model.
3. **The NB model costs nothing on Poisson data:** φ goes large and the field
   is the same.
4. **On the real data the NB term changes almost nothing.** It finds moderate
   overdispersion (humpback φ about 3, PWSD about 1, the range where the
   simulations show no floor-following), and the length-scales stay at the
   floor. So the real floor-following is **not** explained by clustering that
   is independent between segments. The real fields appear to vary on scales
   of 30–50 km, or the clustering is correlated over several 5 km segments,
   which the model cannot tell apart from fine-scale density.

Caveats: one seed; one truth; clustering independent between segments (an NB
term cannot absorb clustering correlated over several segments; that is
short-range field structure).

## Consequences

- Turn the NB option on for both cetaceans in the real fits and the tier-1
  replicates: it protects against clustering at no cost when there is none.
- The real floor-following needs a separate answer: a larger basis budget for
  the visual and joint real fits (the floor follows the budget; on the old
  domain doubling it lowered the posterior with the floor), or accepting it,
  since density estimates have not depended on it.

## Update (2026-10-27): what the floor is, and what is not ruled out

- The floor is the Riutort-Mayol et al. (2023) HSGP rule: a basis of m
  functions per axis, boundary factor c, half-range S, represents
  length-scales down to about 1.75 c S / m (1.75 × 1.5 × 345 / 22 ≈ 41 km
  east-west, 1.75 × 1.5 × 885 / 36 ≈ 65 km north-south at budget 800). Their
  check (the length-scale estimate should sit well above that limit; if not,
  raise m) failed at both budgets, so keeping budget 800 is a computational
  compromise, not evidence that the approximation is adequate.
- Ruled out: clustering independent between 5 km segments. Not ruled out:
  clustering correlated over several segments (10–30 km aggregations), or
  real fine-scale density; one survey pass cannot separate the two.
  Humpback's NB size rising from 2.9 to 8.3 when the basis doubles shows the
  variation is correlated between neighbouring segments.
- Planned tests (ROADMAP_2026-10-21, section 6): an along-track AR(1) segment
  effect alongside the GP, and a humpback fit at budget 3200.

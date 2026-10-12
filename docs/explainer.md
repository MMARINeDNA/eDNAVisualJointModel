---
title: "Joint eDNA and visual survey modelling: what we built and what we tested"
subtitle: "An explainer for the eDNAVisualJointModel repository (draft for review)"
date: "2026-10-11"
toc: true
toc-depth: 3
number-sections: true
---

> **Status of this draft.** Every figure in Sections 2 to 4 was produced from the
> code in this repository (`docs/explainer/`). **Section 5 (the joint model)
> summarises results recorded in the repository by the Phase 5 work; I did not
> re-run them for this document**, and the joint model has not yet been fitted to
> the real data. Items that need your input are marked **[REVIEW]**.

> **Out of date: set-up changed on 2026-10-12** (HSGP4eDNA v1.4,
> eDNAVisualJointModel #82). Every analysis is now restricted to the study
> area **north of 35°N**, one 690 × 1770 km box shared by the simulations and
> the real data. The changes:
>
> | | This document | Current set-up |
> |---|---|---|
> | Simulation box | 500 × 1270 km, about 37.7–49.2°N | 690 × 1770 km, 35°N northwards |
> | Simulated length-scales | 50 km E–W × 300 km N–S | 50 km E–W × 150 km N–S |
> | Simulated transects | 25 (1250 segments) | 35 (2415 segments) |
> | Basis budget | 400 functions | 800 functions |
> | Real visual effort | 2779 segments (12,967 km), from 27.4°N | 1641 segments (7,561 km), north of 35°N |
> | Real-data domain | 1360 × 2600 km | 690 × 1770 km |
>
> **Every number and figure below predates these changes** and will be
> regenerated after the re-runs listed in ROADMAP Phase 6.0. The real eDNA
> bottles (all north of 35°N) are unchanged; the real eDNA fit's numbers still
> change through the new domain and budget.

# What this document is

The aim of the project is to estimate where marine mammals and fish are, using two
very different kinds of data at once:

- **Visual line-transect surveys**: observers on a ship count groups of animals
  and record how far off the trackline each was.
- **Environmental DNA (eDNA)**: water bottles are filtered and the DNA of each
  species is counted by qPCR (hake only) and by metabarcoding (all species).

Both are noisy and indirect measures of the same thing, the density of animals in
space. The modelling idea is to give each species **one latent density surface**
and attach a separate observation model to it for each data source.

This document walks through everything we have built and tested, in the order a
reader would meet it:

1. the model in one page (this section);
2. how test data are simulated, for both data types;
3. the eDNA-only spatial model, on simulated and real data;
4. the visual (distance sampling) spatial model, on simulated and real data;
5. the joint model, on simulated data (real data not yet run);
6. limitations, and how to reproduce the figures.

## The model in one page

For species $s$ at location $x = (X, Y)$ with bottom depth $Z$:

$$
\log \lambda_s(x) \;=\; \mu_s \;+\; f_s(X, Y) \;+\; h_s(Z_\text{bathy})
\qquad \text{(animals per km}^2\text{)}
$$

- $\mu_s$ is an intercept.
- $f_s$ is a **Gaussian process** over position, one per species, with separate
  length-scales east-west and north-south (anisotropic). It is computed with a
  **Hilbert-space approximation (HSGP)**: the GP is written as a weighted sum of a
  fixed set of sine basis functions, so the sampler works with a few hundred
  coefficients instead of an $N \times N$ covariance matrix.
- $h_s$ is a smooth effect of **bottom depth**, a natural spline with 4 degrees of
  freedom.

Each data source then has its own observation model on top of $\lambda_s$:

| Data | Observation model |
|---|---|
| eDNA, qPCR (hake only) | Expected copies in a reaction $=\lambda_s \times$ conversion factor $\times$ litres filtered $\times$ water-column offset $\times$ bottle noise $\times$ aliquot fraction. Detection is $1-\exp(-\kappa \cdot \text{copies})$. Detected replicates have a Ct that falls linearly in the log of copies. |
| eDNA, metabarcoding (all species) | Zero-inflated beta-binomial on read counts, with read shares proportional to each species' copies relative to a "junk" background of all non-target DNA. |
| Visual line transects | Detected groups per segment are Poisson with mean $\lambda_s / E[\text{group size}] \times 2 L \times \text{ESW}$, where ESW is the effective strip half-width implied by the detection function. Perpendicular distances follow the detection function. |

Everything downstream is about whether this structure can be fitted reliably, and
whether its parts behave when they meet real data.

## Which model variants exist

The upstream eDNA repository (HSGP4eDNA) contains several spatial variants: a 2-D
GP over (X, Y) only; a 3-D GP with bottom depth as a third axis; and a 2-D GP plus
a **spline on bottom depth**. This document focuses on the last one, the 2-D HSGP
plus spline, because it is the cheap structure used by the visual and joint models.
The eDNA tests shown here use the scenario with **species-specific water-column
depth preferences** (`bathysp_depthspref`): bottles are taken at six depths, and
each species' eDNA concentration peaks at a different depth.

## Versions and where the code lives

| Piece | Where |
|---|---|
| eDNA model, HSGP functions, eDNA simulators | `external/HSGP4eDNA` (submodule, v1.3, `fecf3c4`) |
| Visual model | `stan/hsgp_visual.stan`, `stan/include/visual_functions.stan` |
| Visual and joint simulators | `R/functions_visual.R`, `R/functions_joint.R` |
| Joint model | `stan/hsgp_joint.stan`, `R/fit_joint.r`, `R/compare_joint.r` |
| Real-data loader | `R/01_load_muri.r`, notes in `notes/muri_data_formatting.md` |
| Figure scripts for this document | `docs/explainer/` |

# Data simulation

Every model in the repository is first tested on data simulated from a known truth,
so that "does the model recover the truth?" has an answer.

## One latent field per species

Three species are simulated, chosen to span how informative each data source is:

| Species | eDNA | Visual | Baseline density | GP marginal SD |
|---|---|---|---|---|
| Pacific hake | qPCR and metabarcoding | not surveyed | 15 per km² | 2.0 |
| Humpback whale | metabarcoding only | yes (400 to 1000 detections) | 0.01 per km² | 1.0 |
| Pacific white-sided dolphin (PWSD) | metabarcoding only | yes (40 to 120 detections) | 0.02 per km² | 1.3 |

All three share the same GP length-scales: 50 km east-west and 300 km north-south,
which mimics an ocean where structure runs along the coast.

![Spatial GP component of the latent log-density for each species (one draw).](figures/sim_field_maps.png)

The bottom-depth effect $h_s$ and the water-column preference are separate,
smooth, known functions:

![Left: the bottom-depth spline effect used as truth (arbitrary test curves, not ecological estimates). Right: the water-column offset on eDNA by sampling depth. Hake and humpback eDNA peak subsurface; PWSD is surface-associated.](figures/sim_depth_effects.png)

## Simulating eDNA observations

Given the density at each bottle, the simulator generates data in the same order
that nature does:

1. **Bottle copies.** Expected copies = density $\times$ a conversion factor
   (hake 10, humpback 200, PWSD 110) $\times$ litres filtered $\times$ the known
   water-column offset. The realised copy count is negative binomial around that
   mean ($\theta = 10$), which is **between-bottle noise**.
2. **Aliquot.** A small fraction of the bottle (2 $\mu$L of a 100 $\mu$L elution)
   goes into each PCR. The copies that land in it are binomial.
3. **qPCR (hake).** Detection with probability $1 - e^{-\kappa c}$, $\kappa = 0.85$,
   then a Ct that falls linearly in $\ln c$ with measurement noise.
4. **Metabarcoding.** Reads are split among the three species and a **junk**
   category in proportion to their copies; junk varies from bottle to bottle.

![Simulated eDNA observations. Left: qPCR detection against expected copies in the aliquot (red: the true hurdle). Middle: Ct against copies (red: the true standard curve). Right: share of metabarcoding replicates with any reads, by sampling depth.](figures/sim_edna_obs.png)

The whales sit near the detection limit: about 0.04 to 0.08 expected copies in a
typical aliquot, so only about 10% of aliquots contain any whale DNA. That is the
reason eDNA alone constrains the whale fields only weakly, and one reason the
visual data matter.

## Simulating visual surveys

The same field is observed along 25 east-west transects split into 10 km segments.
In each segment the number of animal groups in the strip is Poisson; each group is
at a uniform random distance out to the truncation distance (5 km) and is detected
with probability given by the detection function.

- **Humpback:** half-normal detection with scale 2.5 km, groups of 1 to 8.
- **PWSD:** half-normal with scale 1.5 km at the mean group size, wider for larger
  groups, groups of 2 to 452. Group sizes are resampled from the real survey's
  pools (`data/grpsz/`).
- **Sea state (optional):** detection falls with Beaufort, with a spatial trend in
  sea state so that ignoring it biases the field (details below).

![Detected perpendicular distances against the true detection function (top), and the empirical group-size pools the simulator resamples (bottom).](figures/sim_visual_obs.png)

### How sea state is simulated

With `sea_state = TRUE` each segment gets a Beaufort value, and that value changes the
detection scale of every group seen from it:

$$
b = 3 + 1.0\left(\frac{2Y}{Y_\max} - 1\right) + 0.75\left(1 - \frac{2X}{X_\max}\right) + \varepsilon,
\qquad \varepsilon \sim N(0, 0.5),\quad b \in [0, 6]
$$

$$
\sigma \;=\; \sigma_\text{det}\; e^{\beta_\text{size}(s - \bar s)}\; e^{\beta_\text{bft}(b - 3)},
\qquad g(x) = \exp\!\left(-\frac{x^2}{2\sigma^2}\right)
$$

- **The trend** makes the north (large $Y$) and the west (small $X$, the offshore side)
  rougher, with segment-level noise on top. The mean is 3 and the trend spans about
  $\pm 1.75$ Beaufort units.
- **The effect** is a log-linear multiplier on the detection scale, equal to 1 at
  Beaufort 3. The truths are $\beta_\text{bft} = -0.08$ (humpback) and $-0.15$
  (PWSD) per Beaufort unit.
- **Density does not depend on sea state.** Only the chance of seeing a group does.
  Because sea state is correlated with position, a model that ignores it reads the
  rough north and west as low density; the validation measured field shifts of up to
  0.58 log units across the survey.

![How sea state is simulated and what it does to detectability: the chain from Beaufort to detection (top), the expected sea-state surface (bottom left), detection functions at Beaufort 1, 3 and 5 (bottom middle), and the average detection probability in a 5 km strip (bottom right).](figures/sim_sea_state.png)

Averaged over a 5 km strip, going from Beaufort 1 to 5 lowers the detection
probability from 0.67 to 0.52 for humpback and from 0.50 to 0.28 for PWSD (at the
reference group size). The simulator uses the same log-linear form that the model
fits, so a fit to simulated data has the right structure by construction; real sea
state is also not a smooth spatial trend, and rough days are often not surveyed.

## The joint simulator

`simulate_joint()` draws **one** field over the union of eDNA station and transect
segment locations, then generates both data types from it. Each half has exactly
the shape of its single-source simulator's output, so every single-source tool
works on it unchanged.

![One joint simulation: 200 eDNA stations and 25 transects on the same latent field.](figures/sim_design.png)

The check that the two halves really see the same field (from the Phase 4
validation, two seeds):

| Check | Result |
|---|---|
| Field correlation between stations and segments within 10 km | 0.995 to 0.997 |
| eDNA-only fit on the eDNA half: field R² (hake / humpback / PWSD) | 0.98 / 0.73 / 0.68 and 0.98 / 0.58 / 0.65 |
| Visual-only fit on the visual half: field R² (humpback / PWSD) | 0.90 / 0.54 and 0.87 / 0.58 |

PWSD is weak in both sources, which is where a joint fit should help most.

# The eDNA-only spatial model

## Structure

The model is `stan/hsgp_2d_bathysp.stan` in HSGP4eDNA. Its main ingredients, beyond
the observation model above, are:

- **Design-based GP priors.** Priors on the length-scales and marginal SD are
  built from the **sample locations only**, never from observed counts. The
  length-scale prior is inverse-gamma with 1% of its mass below a floor
  $\ell_\min$ (the larger of the median sample spacing and the shortest scale a
  400-function basis can represent) and 1% above the domain extent. The basis size
  is chosen to match $\ell_\min$.
- **Between-bottle noise.** A negative-binomial multiplier on each bottle's expected
  copies, with one spread $\theta$ shared across species. Leaving it out makes the
  GP absorb bottle noise as spurious small-scale structure.
- **Junk noise.** The junk background also varies between samples, with its own
  spread $\theta_\text{junk}$. Without it the species' noise term absorbs junk
  variation and $\theta$ reads about 7 against a true 10.
- **Known water-column offset.** Depth preferences enter as a fixed offset per
  sample and species; the model does not estimate them.

![eDNA-only model. Each species' GP hyperparameters, intercept and bottom-depth spline set its density field. With fixed conversion factors, volumes and water-column offsets, the field sets each bottle's expected copies. The bottle's species copies (spread θ, shared by all species) and its junk copies (spread θ_junk) are unobserved negative-binomial draws. Its qPCR and metabarcoding replicates are observed and share those copies.](figures/model_edna.svg)

## Test on simulated data

**Setup.** A `bathysp_depthspref`-style simulation (seed 303) with 40 stations and 6
sampling depths (0 to 500 m): 235 bottles (samples deeper than the sea floor are
dropped), with 705 qPCR and 705 metabarcoding replicates. The basis has 144
functions (12 by 12), set by the design-based prior rule with a budget of 150.
Three chains of 400 + 400 iterations.

| | |
|---|---|
| Runtime (3 chains) | 84 min |
| Divergences | 0 |
| Transitions at the maximum treedepth (10) | 468 of 1200 (39%; 400 in chain 1, 68 in chain 2, 0 in chain 3) |
| max Rhat (hyperparameters, $\theta$) / min bulk ESS | 1.02 / 179 |
| E-BFMI per chain | 0.53, 0.36, 0.65 |

> **[REVIEW]** This is a **reduced version** of the reference scenario (the validated
> run has 200 stations, 1185 bottles, a 400-function basis, and takes about 15
> hours; a first attempt at 80 stations was projected to take about 20 hours on
> this shared machine). The sampler is not clean: 39% of transitions were
> truncated at treedepth 10 and the minimum ESS is only 179. The results are
> indicative; the 1185-bottle numbers in the table below are the reference.

### Recovery of the latent field and hyperparameters

![GP hyperparameters: posterior (blue, 95% interval) against truth (red x). The grey band is the central 95% of the prior, so a posterior well inside it shows the data, not the prior, set the value.](figures/edna_sim_hyper.png)

![Latent log-density at the 40 sampling stations: posterior mean (90% interval) against truth.](figures/edna_sim_field.png)

![Bottom-depth effect: posterior (blue band) against truth (red dashed).](figures/edna_sim_spline.png)

| | Hake | Humpback | PWSD |
|---|---|---|---|
| Field R² at the stations | 0.99 | 0.86 | 0.24 |
| 90% intervals covering the truth | 88% | 88% | 78% |
| Hyperparameters covered (of 4) | 4 | 3 (intercept $\mu$ misses) | 3 ($\ell_x$ misses) |

Ten of twelve hyperparameter intervals contain the truth. Hake is recovered very
well (R² 0.99). The humpback field is recovered with the right shape (R² 0.86). PWSD
is poorly recovered (R² 0.24): its eDNA is surface-associated and almost never
detected below 50 m, so 40 stations carry little information. The two misses are
both small: the humpback intercept's interval is [$-4.56$, $-2.85$] against a
truth of $-4.61$, and PWSD's east-west length-scale is 96 km [56, 163] against 50.
Intercepts are partly confounded with the realised field's mean, which is a known
property of this model (HSGP4eDNA limitation L4). The whale intercepts are
also pulled up toward the prior mean (see below).

### Prior and posterior predictive checks

The **prior predictive** check asks whether the priors, pushed through the model,
produce data of a sensible scale. I draw parameters from the priors, build the
field at the real bottle locations (an exact GP draw, which the HSGP approximates),
add bottle noise, and look at the implied eDNA signal.

![Prior predictive. Top: the typical (median) expected copies per aliquot across bottles, prior draws against the simulated truth (red). Bottom: the hake qPCR detection rate implied by the priors, with the observed value (red) and the posterior prediction (blue dashed).](figures/edna_sim_prior_pc.png)

- **Hake:** the truth sits in the middle of the prior predictive, and the implied
  qPCR detection rate spans 0.4 to 1.0 (5th to 95th percentile), containing the
  observed 0.95.
- **Whales:** the truth sits in the **lower tail** (humpback) or **below** (PWSD) of
  what the priors imply. The default intercept prior, normal(2, 1.5) on log
  density, was set for hake. The posterior intercepts for the whales are
  correspondingly biased high (humpback $-3.9$ against $-4.6$; PWSD $-3.0$ against
  $-3.9$), although the field shapes are recovered.

The **posterior predictive** check asks whether data replicated from the fitted
model look like the data we saw.

![Posterior predictive check for metabarcoding: share of replicates with reads, by species and sampling depth. Blue: model-implied probability (median, 90% over posterior draws). Red x: observed.](figures/edna_sim_ppc_mb.png)

![Posterior predictive check for qPCR: detection rate by depth (left) and observed against predicted Ct (right).](figures/edna_sim_ppc_qpcr.png)

- **Metabarcoding.** The depth profiles of the share of replicates with reads, which
  are the depth preferences, are reproduced for all three species: hake is nearly
  always positive, humpback peaks at 150 m, PWSD falls to zero below 150 m.
- **qPCR.** Detection rates by depth follow the observed shape, but the model
  predicts 1 to 2 percentage points fewer detections than observed at every
  depth, a small, consistent under-prediction. Ct is predicted with an RMSE of
  0.63 cycles; the largest errors are at the highest Cts, where the observed
  values are lower than predicted.

> **[REVIEW] Upstream gap.** The metabarcoding posterior predictive draws in
> HSGP4eDNA's generated quantities (`pp_mb_reads`) leave out the zero-inflation
> term that the likelihood contains, so they predict reads for essentially every
> replicate. The metabarcoding check above therefore computes the probability of
> any reads directly from posterior draws, using the likelihood's own formula.
> Adding the zero-inflation step to the generated quantities in HSGP4eDNA would fix
> this at the source; I did not change that repository.

### What HSGP4eDNA's validation adds

For context, the upstream validation of the same structure (3 chains of 800 + 800,
v1.3 defaults):

| Scenario fitted with 2-D HSGP + spline | N | Runtime | Divergences | max Rhat | Field R² (hake / humpback / PWSD) | Hyperparameters covered |
|---|---|---|---|---|---|---|
| `bathysp_surface` (matched) | 200 | 49 min | 0 | 1.01 | 0.98 / 0.59 / 0.74 | 12 of 12 |
| `bathysp_depthspref` (matched) | 1185 | 924 min | 0 | 1.01 | 1.00 / 0.94 / 0.76 | 12 of 12 |
| `bathygp_surface` (mis-specified: depth inside a GP) | 200 | 123 min | 0 | 1.30 | 0.70 / 0.81 / 0.76 | n/a, not converged |

The last row matters: when the truth is not a spline in bottom depth, the model
**says so**. Chains disagree (Rhat 1.30), the hake length-scales fall below the
prior floor, and $\theta$ collapses to 2.8. Mis-specification is detectable without
knowing the truth.

## Test on real data

**What is available.** The MURI eDNA data have 551 bottles at 182 stations and 19
distinct sampling depths (172 at the surface, then mostly 50, 150, 300 and 500 m).
There are 2663 hake qPCR replicates and 1709 metabarcoding replicates.

**What is not.** The qPCR file holds copies per $\mu$L, not Ct, and there is no
standard curve (slope, intercept, detection efficiency $\kappa$, Ct noise). The
model's qPCR part is built on Ct, so **it cannot be fitted yet**. Volumes filtered,
aliquot and elution volume are also missing, and the eDNA files carry no dates.

**What I did.** The model was fitted to metabarcoding only: all qPCR rows were
dropped except one non-detect (the model's data block requires at least one), the
water-column offset was set to 1 for all bottles, and the conversion factors and
volumes are the simulator's placeholders. So the fit says **where relative
density is high**, not how many animals there are.

The fit used 2 chains of 400 + 400 iterations and a basis of 399 functions (7 by
57; the design prior stretches the basis along the coast because the 551 bottles
lie in a long thin band).

| | |
|---|---|
| Bottles / metabarcoding replicates | 551 / 1709 |
| Runtime (2 chains) | 247 min |
| Divergences / max-treedepth hits | 0 / 0 |
| max Rhat / min bulk ESS | 1.02 / 80 |
| E-BFMI per chain | 0.41, 0.49 |
| Bottle-noise $\theta$ (median, 90% interval) | 0.24 [0.21, 0.26] |
| Junk-noise $\theta_\text{junk}$ | 0.65 [0.57, 0.74] |

![Posterior mean latent density at the stations, by species. Units are nominal (placeholder conversion factors and volumes), so only relative patterns are meaningful.](figures/edna_real_map.png)

![GP hyperparameters on the real bottles: posterior against the design-based prior.](figures/edna_real_hyper.png)

![Posterior predictive check: share of metabarcoding replicates with reads, by species and sampling depth (binned). Blue: model-implied probability. Red x: observed.](figures/edna_real_ppc_mb.png)

**What the checks say.**

- **The model reproduces the observed read frequencies.** For all three species and
  all depth classes the observed share of replicates with reads is inside, or within
  about a percentage point of, the predicted interval. The one visible miss is
  humpback at 300 m: none of 234 replicates had reads, against about 3.7 expected.
  This is an in-sample check, and with free per-bottle noise terms it is easy to
  pass.
- **The noise terms are very large.** $\theta$ = 0.24 means a coefficient of
  variation of about 2 in the bottle-to-bottle multiplier; on matched simulations
  $\theta$ is about 9 (truth 10). HSGP4eDNA documents a $\theta$ far below its
  matched value as a sign that the noise terms are absorbing structure the model
  lacks (its mis-specified example read 2.8).
- **The length-scales are short and at the floor.** The along-coast length-scale is
  about 34 to 43 km against a prior floor of 27 km, and `floor_check()` flags
  all six species-axis combinations (the posterior puts 16 to 65% of its mass within
  1.25 $\ell_\min$, against 3 to 6% under the prior).
- **Sampling is thin.** The minimum bulk ESS of 80 and E-BFMI near 0.4 come from
  only two chains.

> **[REVIEW] How to read this fit.** I would not interpret the map or the
> hyperparameters yet. Several things the data have and the model does not, each
> of which could push variation into the noise terms and the short length-scales:
> the 19 sampling depths are treated as one (offset 1), whereas the simulations
> show how strongly depth shapes eDNA; the metabarcoding replicates were run at
> template dilutions from 1 to 100-fold; the qPCR information (the most
> informative data for hake) is absent; and the conversion factors and volumes are
> placeholders. The honest summary is that the model **runs on the real bottles
> and reproduces their read frequencies**, and that the diagnostics point to
> structure it does not yet represent.

> **[REVIEW]** Open decisions before the real eDNA fit is more than a sanity check:
> where the Ct values and standard curve are; whether volumes filtered vary by
> bottle (the model takes one value); and how to treat the 19 sampling depths
> (a modelling choice: offset, covariate, or separate surfaces).

# The visual (distance sampling) spatial model

## Background: what a line-transect analysis does

On a line transect, animals close to the line are nearly always seen and animals
far away are often missed. A **detection function** $g(x)$ describes the
probability of detection at perpendicular distance $x$. Integrating it gives the
**effective strip half-width** (ESW): a strip of that half-width with perfect
detection would see the same number of animals. Density follows from counts
divided by $2 \times L \times \text{ESW}$.

A traditional analysis (for example Len Thomas's CalCOFI analysis, which we
followed) makes these choices, and so does our model:

| Choice | Traditional practice | Our model |
|---|---|---|
| Truncation | Drop the largest 5 to 10% of distances, per species | 95th percentile per species |
| Key function | Half-normal or hazard-rate, by AIC, with caution for species that bow-ride | Half-normal; **hazard-rate option** (used for humpback on real data) |
| Covariates | Group size, Beaufort, swell, observer platform; chosen by AIC | Group size and Beaufort; chosen by AIC on this data set |
| Size bias | Larger groups are seen farther away; correct for it | Group-size distribution is modelled and size-bias corrected |
| g(0) | Often corrected for animals missed on the line | **Assumed 1** (a limitation) |
| Spatial model | Often none | 2-D HSGP + bottom-depth spline on log density |

## Structure of the visual model

The model is `stan/hsgp_visual.stan`, one species at a time:

- **Field:** $\log\lambda = \mu + f(X,Y) + B(Z)\beta$, the same structure and units
  as the eDNA model, so the two can share it.
- **Encounters:** segment counts are Poisson with mean
  $\lambda / E[s] \times 2L \times \text{ESW}$.
- **Distances:** truncated half-normal or hazard-rate detection with
  $\log\sigma = \log\sigma_0 + \beta_\text{size}(s - \bar s) + \beta_\text{bft}(\text{Beaufort} - 3)$.
- **Group sizes:** zero-truncated negative binomial (humpback) or log-normal
  (PWSD), with the size-bias correction.
- **Priors:** the same design-based GP priors as the eDNA model, built from
  segment locations only.

![Visual model, one species. The density field, expected group size and effective strip half-width set the Poisson number of groups detected on each segment. For each detected group, the group size and perpendicular distance are observed; the distance follows the truncated detection function, whose scale depends on group size and optionally Beaufort. There is no per-segment latent draw.](figures/model_visual.svg)

## Group size: where it comes from, how it enters, and how it is biased

Group size matters twice in a line-transect model. It sets how many animals sit
behind each detected *group*, and (when larger groups are easier to see) it changes
how far out a group is detected. Both go through one parametric distribution of
group sizes and the effective strip half-width.

### Where the group sizes come from

| Stage | Source |
|---|---|
| **Real data** | The observers' best estimate of group size (`gsspbest`) for each *detected* group: on effort, with a perpendicular distance, within the truncation distance |
| **Simulator truth** | The empirical pools in `data/grpsz/` (all of the survey's sightings of the species: 683 humpback, 110 PWSD), resampled **as the population** of groups present |
| **Model** | A parametric population distribution fitted to the detected sizes: zero-truncated negative binomial (humpback) or log-normal (PWSD), on sizes $1,\dots,S_\max$ |

### How group size enters the model

Let $p_k$ be the population probability of a group of size $k$, and let
$\sigma_k$ be the detection scale for such a group in sea state $b$:

$$
E[s] \;=\; \sum_{k=1}^{S_\max} k\,p_k ,
\qquad
\sigma_k \;=\; \sigma_0\, e^{\beta_\text{size}(k-\bar s)}\, e^{\beta_\text{bft}(b-3)} ,
\qquad
\text{ESW}(\sigma_k) \;=\; \int_0^w g(x \mid \sigma_k)\,dx .
$$

The strip width is **averaged over the population of group sizes**, not evaluated at
the average size:

$$
\overline{\text{ESW}}(b) \;=\; \sum_{k=1}^{S_\max} p_k \,\text{ESW}\big(\sigma_k(b)\big).
$$

Segment counts then follow

$$
E[\text{detected groups}_j] \;=\; \frac{\lambda_j}{E[s]} \;\times\; 2L_j \;\times\; \overline{\text{ESW}}(b_j),
$$

where $\lambda_j$ is animal density (animals per km²), so $\lambda_j / E[s]$ is the
density of *groups*. Each detected group contributes its distance $x_i$ and size
$s_i$ through

$$
p(x_i, s_i \mid \text{detected}) \;=\; \frac{g(x_i \mid s_i)\, f(s_i)}{\overline{\text{ESW}}(b_i)},
$$

so the sizes of detected groups inform both the population distribution $f$ and the
size effect $\beta_\text{size}$.

**Segments with no detections.** $\overline{\text{ESW}}$ does not use a segment's own
detections. It depends only on the segment's covariate level (its binned Beaufort,
where used) and on the global parameters ($\sigma_0$, $\beta_\text{size}$,
$\beta_\text{bft}$, the hazard-rate shape, the group-size parameters). All segments
at the same sea-state level share the same value, whether or not anything was seen,
and a zero count is evidence of low density in proportion to that expected count.

**Without a size effect on detection** ($\beta_\text{size}$ off), $\text{ESW}(\sigma_k)$
is the same for every $k$, group size drops out of the strip width, and
$E[s]$ is simply the mean of the detected sizes (the model has no way to know
detected groups are unrepresentative).

**With a size effect**, larger groups have a wider strip, so the *detected* groups
have a mean size above the population mean:

$$
\text{mean size of detected groups} \;=\; \frac{\sum_k k\,p_k\,\text{ESW}_k}{\sum_k p_k\,\text{ESW}_k} \;>\; E[s],
\qquad
\text{effective strip width per animal} \;=\; \frac{\overline{\text{ESW}}}{E[s]} .
$$

The model estimates $E[s]$ from the detected sizes by undoing that weighting.

*Implementation.* The sum over $k$ is binned (sizes 1 to 30 individually, then 10%
log bins; it agrees with the full sum to about $10^{-5}$), the hazard-rate integral
is a 400-cell midpoint rule, and the half-normal integral is closed-form.

### Is it biased?

Yes, in several ways. The first three are the important ones.

| Source of bias | What happens | Handled? |
|---|---|---|
| **Size-biased detection** | Larger groups are seen farther out, so detected groups over-represent large ones | Only through $\beta_\text{size}$, and only if that term is in the model and its log-linear form is adequate (an extrapolation for the largest dolphin schools) |
| **No size effect assumed** | If $\beta_\text{size}$ is off, $E[s]$ equals the detected mean by construction | Not corrected; use the size effect where the data support it |
| **Simulator truth is a detected sample** | The pools are made of detected groups, then used as the population | Simulation tests consistency, not realism |
| **Observer error** | Best estimates of large schools are rough and rounded; "groups" can be subgroups | Not modelled; enters $E[s]$ one-for-one |
| **One distribution everywhere** | Mean group size is global, so spatial variation in group size is absorbed into the density surface | Not modelled |
| **Tail of the log-normal** | The PWSD support runs to 1000 against a largest observed group of 452 | Sensitivity check below |
| **Missing sightings** | 73 on-effort sightings have no perpendicular distance and are excluded | If missingness relates to size, biased |

Because density is group density times $E[s]$, any bias in $E[s]$ is a bias in
animal density of the same proportion.

### What the real-data fits show

Four fits on the current set-up (north of 35°N, 1641 segments, 792-function basis,
2 chains of 500 + 500, no divergences, max Rhat at most 1.01):

| | Humpback, no size effect | Humpback, with size effect | PWSD, $S_\max$ = 1000 | PWSD, $S_\max$ = 500 |
|---|---|---|---|---|
| Detections | 548 | 548 | 89 | 89 |
| Mean size of detected groups | 1.70 | 1.70 | 46.6 | 46.6 |
| Mean size within 1 km of the line | 1.56 (n = 180) | 1.56 | 37.4 (n = 39) | 37.4 |
| Model population mean $E[s]$ (90% interval) | 1.70 [1.63, 1.78] | 1.69 [1.61, 1.78] | 37.0 [29.7, 47.5] | 36.5 [29.2, 45.1] |
| $\beta_\text{size}$ (90% interval) | at prior (off) | 0.029 [$-0.043$, 0.101] | 0.0041 [0.0009, 0.0128] | 0.0039 [0.0011, 0.0123] |
| Mean animal density (per km²) | 0.0227 | 0.0225 | 0.093 | 0.093 |
| Mean detection probability | 0.44 | 0.44 | 0.51 | 0.50 |
| Runtime (2 chains) | 13 min | 18 min | 24 min | 22 min |

- **PWSD has a clear size bias, and the model corrects it.** The detected groups
  average 46.6 animals, but the population mean is estimated at 37.0, a correction
  of about $-21\%$. An independent check agrees: groups seen within 1 km of the
  line, where detection is nearly certain regardless of size, average 37.4. Beyond
  3 km the mean is 74 (n = 13). Without the correction PWSD density would be about
  26% too high.
- **The tail does not matter much.** Cutting the support from 1000 to 500 changes
  $E[s]$ by $-1.3\%$ and density by $-0.5\%$.
- **Humpback has little size bias to correct.** Groups are mostly singles. With the
  size effect on, $\beta_\text{size}$ is 0.029 with an interval that includes zero,
  $E[s]$ moves by $-1\%$ and density by $-0.8\%$. The data are not informative about
  it (548 detections, most of size 1 or 2); the near-line mean is a little lower
  than the model's (1.56 against 1.69), within sampling noise at this sample size.

> **Correction.** The AIC selection for humpback picked a group-size effect on
> detection ($\Delta$AIC $-3.2$), and the notes and the real-data tables in this
> document describe humpback as "hazard-rate, group size". The earlier real-data
> humpback fits did **not** include it: the loader never switched
> `use_size_covar` on for humpback (the simulator's humpback truth has none). The
> loader now sets it, and the processed objects are regenerated. The comparison
> above shows the omission changed humpback's mean group size and density by about
> 1%, so the earlier numbers stand to that accuracy. (All the real-data results
> above this point in the document also predate the 35°N set-up.)

## Test on simulated data

**Setup.** One survey (seed 1): 25 transects of 10 km segments, 549 humpback and
78 PWSD detections within 5 km. Basis of 392 functions (28 by 14, set by the
design-based prior rule); 3 chains of 800 + 800 iterations.

| | Humpback | PWSD |
|---|---|---|
| Detections | 549 | 78 |
| Runtime (3 chains) | 4.9 min | 14.4 min |
| Divergences / max-treedepth hits | 0 / 0 | 0 / 0 |
| max Rhat / min bulk ESS | 1.004 / 313 | 1.005 / 778 |
| Latent-field R² at the segments | 0.92 | 0.70 |

![Humpback: posterior mean log-density against the truth at the 1250 segments.](figures/visual_sim_humpback_field.png)

![PWSD: posterior mean log-density against the truth at the 1250 segments.](figures/visual_sim_pwsd_field.png)

**Recovery of the parameters** (posterior median and 95% interval):

| Quantity | Humpback truth | Humpback posterior | PWSD truth | PWSD posterior |
|---|---|---|---|---|
| Detection scale $\sigma$ (km) | 2.5 | 2.63 [2.40, 2.93] | 1.5 | 1.20 [1.02, 1.46] ✗ |
| Mean density (per km²) | 0.0136 | 0.0121 [0.0107, 0.0137] | 0.061 | 0.082 [0.060, 0.111] |
| $\ell_x$ (km) | 50 | 43 [30, 82] | 50 | 56 [34, 83] |
| $\ell_y$ (km) | 300 | 274 [206, 364] | 300 | 286 [159, 463] |
| GP marginal SD | 1.0 | 0.97 [0.73, 1.44] | 1.3 | 1.37 [0.90, 2.14] |
| Mean group size | 1.72 | 1.71 [1.64, 1.80] | 42.7 | 40.3 [31.5, 52.9] |

Eleven of twelve intervals contain the truth. The miss is PWSD's detection scale,
just outside its interval (1.46 against a true 1.5), which also happened in one of
the five replicates in the repository's validation study.

![Humpback: GP hyperparameters, posterior against prior (grey) and truth (red x).](figures/visual_sim_humpback_hyper.png)

![Humpback: posterior predictive counts per transect, fitted detection function, and prior against posterior predictive total.](figures/visual_sim_humpback_ppc.png)

![PWSD: the same three checks.](figures/visual_sim_pwsd_ppc.png)

**Reading the checks.**

- **Counts.** The observed number of detections per transect line sits inside the
  replicated 90% interval for nearly every line (about 90% by construction, so a
  few misses are expected).
- **Detection function.** The fitted curve follows the histogram of detected
  distances, and the bars nearest the trackline are the tallest, as a
  half-normal predicts. (Histograms in this document start their first bin at
  zero. An earlier version let the plotting default centre a bin on zero, which
  halves the first bar and looks like a deficit of detections near the line; the
  counts themselves, 87 in the first 0.5 km for humpback against 91 expected,
  were never short.) For PWSD (78 detections) the first two bars sit above the
  curve, within the noise of a small sample.
- **Prior against posterior.** The priors alone imply anywhere from tens to
  hundreds of thousands of detected groups, so they are weak. The data then
  collapse that to a narrow range around what was observed.
- **Hyperparameters.** All posteriors sit well inside the prior band, so the
  data, not the prior, determine them. The east-west length-scale ($\ell_x$) is
  the least certain parameter for both species.

### Replicated validation

The repository's multi-replicate study (5 simulated surveys, 4 chains of 1000 +
1000, design-based priors) gave:

| Species | Detections | Divergences | max Rhat | Field R² | Coverage of density, $\sigma$, group size | Coverage of GP SD, $\ell_x$, $\ell_y$ |
|---|---|---|---|---|---|---|
| Humpback | 549 to 917 | 0 in 5 fits | at most 1.01 | 0.87 to 0.93 | 5/5, 5/5, 5/5 | 5/5, 5/5, 5/5 |
| PWSD | 41 to 100 | 1 each in 2 fits | at most 1.00 | 0.60 to 0.79 | 5/5, 4/5, 5/5 | 5/5, 4/5, 5/5 |

With sea state simulated (3 replicates), fitting without the Beaufort covariate
reads poor sighting conditions as low density: the fitted field shifts by up to
0.58 log units across the survey. The covariate recovers its true effect in 3 of 3
fits.

## Real data: the MURI line-transect survey

**Data.** The 2018 CCES cruise: 2779 effort segments of about 5 km (12,967 km of
effort), with sightings of humpback and PWSD. *(Now restricted to north of
35°N: 1641 segments, 7,561 km. The fits below used all the effort; see the
note at the top.)* Perpendicular distances come from
the survey's own records (`SITEINFO_CCE_2018_Pp_Lo_Mn.csv`; km).

**Decisions made along the way** (all in `notes/muri_data_formatting.md`):

| Decision | Choice |
|---|---|
| Truncation | 95th percentile per species: humpback 6.45 km, PWSD 4.64 km (548 and 89 detections kept) |
| Detection key | Humpback hazard-rate (about 20 AIC points better than half-normal); PWSD half-normal |
| Covariates | Humpback: group size (selected by AIC; see the correction in Section 4.3). PWSD: group size and Beaufort (visibility dropped by choice) |
| Domain | Bounding box of all effort and eDNA, 1360 by 2600 km |
| Sightings without a distance | Excluded (69 humpback, 4 PWSD) |

![Humpback hazard-rate vs half-normal: where the choice matters. Real data distances for both species, with the Distance-package fits.](figures/visual_key_choice.png)

### Fit to the real data

Both species were fitted with the spatial model on all 2779 segments (basis of 391
functions, 2 chains of 500 + 500, design-based priors built from the segment
locations). These are single fits, shorter than the validation runs, and should
be read as a working demonstration rather than a final analysis.

| | Humpback | PWSD |
|---|---|---|
| Detection function | hazard-rate, group size | half-normal, group size and Beaufort |
| Detections used | 548 | 89 |
| Runtime (2 chains) | 25 min | 31 min |
| Divergences / max-treedepth hits | 0 / 0 | 0 / 0 |
| max Rhat / min bulk ESS | 1.008 / 387 | 1.004 / 459 |
| Mean density (animals per km², 90% interval) | 0.0135 [0.011, 0.016] | 0.056 [0.041, 0.075] |
| Mean group size | 1.71 [1.64, 1.79] | 37.9 [30.1, 47.6] |
| Detection-function shape | $b$ = 0.96 [0.77, 1.19] | Beaufort effect $-0.12$ per unit [$-0.29$, $0.006$] |
| $\ell_x$, $\ell_y$ (km, medians) | 80, 113 | 86, 128 |
| GP marginal SD (median) | 3.6 | 2.9 |

![Humpback: posterior mean density along the effort, with detections (red circles).](figures/visual_real_humpback_map.png)

![PWSD: posterior mean density along the effort, with detections (red circles).](figures/visual_real_pwsd_map.png)

![Humpback on real data: posterior predictive counts per 100 km latitude band, fitted detection function, and prior against posterior predictive total.](figures/visual_real_humpback_ppc.png)

![PWSD on real data: the same three checks.](figures/visual_real_pwsd_ppc.png)

**Reading the checks.**

- **Counts.** Observed counts per 100 km latitude band fall inside the replicated
  90% intervals almost everywhere for both species, including the sharp contrast
  between the empty south (Y below about 800 km, southern California and further
  south) and the dense regions near central California (Y of about 1000 km) and
  from southern Oregon to Vancouver Island (Y of about 1700 to 2400 km).
- **Detection function.** For humpback the fitted hazard-rate curve passes through
  the tallest bars near the trackline and then falls away with a long tail,
  which is the shape the half-normal could not reproduce (see the detection-key
  comparison above); the
  90% band is widest near zero, where few detections constrain it. For PWSD the
  fit tracks the histogram.
- **Prior against posterior.** As on simulated data, the priors are weak and the
  data narrow the predicted total to the observed value.

![GP hyperparameters on real data (humpback): posterior against prior.](figures/visual_real_humpback_hyper.png)

![GP hyperparameters on real data (PWSD): posterior against prior.](figures/visual_real_pwsd_hyper.png)

**Refit at twice the basis budget.** The check below flagged the first fits, so both
species were refitted with `M_max` = 800 (a basis of 792 functions; the floor
$\ell_\min$ falls to 75 km east-west and 106 km north-south). Neither fit had divergences
(max Rhat 1.005 and 1.004).

| | Humpback, 391 functions | Humpback, 792 functions | PWSD, 391 functions | PWSD, 792 functions |
|---|---|---|---|---|
| Runtime | 25 min | 42 min | 31 min | 54 min |
| Mean density (per km², 90% interval) | 0.0135 [0.0113, 0.0162] | 0.0135 [0.0112, 0.0161] | 0.056 [0.041, 0.075] | 0.057 [0.041, 0.077] |
| Mean group size | 1.71 | 1.71 | 37.9 | 37.7 |
| Detection shape / Beaufort effect | $b$ = 0.96 | $b$ = 0.97 | $-0.12$ | $-0.10$ |
| $\ell_x$, $\ell_y$ (km, medians) | 80, 113 | 58, 85 | 86, 128 | 64, 93 |
| Floor $\ell_\min$ (km) | 106, 150 | 75, 106 | 106, 150 | 75, 106 |
| Floor check flags | both axes | both axes | both axes | both axes |
| GP marginal SD (median) | 3.6 | 3.6 | 2.9 | 2.8 |

![Humpback at the doubled basis: GP hyperparameters, posterior against prior.](figures/visual_real800_humpback_hyper.png)

![PWSD at the doubled basis: GP hyperparameters, posterior against prior.](figures/visual_real800_pwsd_hyper.png)

> **[REVIEW] The length-scale posteriors track the floor, and density does not
> depend on it.** At both budgets `floor_check()` flags both axes of both species:
> the posterior length-scales (58 to 128 km) lie below the floor ($\ell_\min$ = 75 to
> 150 km), and doubling the basis lowered the floor and the posterior with it
> (humpback $\ell_x$ 80 to 58 km). So a larger basis alone does not make the flag
> go away. Density, detection and group-size estimates did not change (mean
> humpback density 0.0135 both times; PWSD 0.056 and 0.057). Two things are
> consistent with this pattern, and I have not separated them: real
> sightings are more clustered than a smooth GP with Poisson counts allows, so
> the GP uses short length-scales to chase hot spots (the same mechanism that
> HSGP4eDNA documents for unmodelled bottle noise); or the real field really
> has structure on scales of tens of kilometres. The visual model has no
> overdispersion term for segment counts. The marginal SD is also above the
> prior's central 95% (upper bound 2.7) in all four fits.

> **[REVIEW]** The real-data effort extends well beyond the eDNA region: eDNA
> bottles cover about 20% of the width and 45% of the height of the domain.
> *(On the 35°N study domain they cover about 39% of the width and 66% of the
> height.)* A
> visual-only fit is informative everywhere; a joint fit will only borrow
> strength in the eDNA region.

# The joint eDNA and visual spatial model

The joint model (`stan/hsgp_joint.stan`, driven by `R/fit_joint.r`) was added to
`main` after the rest of this document was written, in Phase 5 of the roadmap.
This section summarises what the repository records about it. The numbers are
from `README.md` and `ROADMAP.md`, produced by that work.

## Structure

Each species has **one** latent field, a 2-D HSGP plus a bottom-depth spline in
animals per km², exactly as in the single-source models. Two observation models
share it:

- **eDNA:** HSGP4eDNA's v1.3 model, run unchanged: qPCR for hake, zero-inflated
  beta-binomial metabarcoding for all species, and the between-bottle and junk
  noise terms.
- **Visual:** the line-transect model for humpback and PWSD, with the detection
  key (half-normal or hazard-rate) and per-species truncation of the visual
  model. A species map ties the cetaceans' eDNA fields to their line-transect
  data.

![Joint model. One density field per species feeds both observation models. The eDNA side is the HSGP4eDNA model above (bottle and junk noise, replicates), with conversion factors fixed or, optionally, free for the cetaceans. The visual side is the line-transect model for humpback and PWSD. The likelihood switches drop either side without changing anything else.](figures/model_joint.svg)

Data switches turn either likelihood off, and optionally **free the cetaceans'
eDNA conversion factors**. The eDNA-only, visual-only and joint fits therefore use
identical priors, basis and starting values, and differ only in the data.

**Exactness checks** (recorded in the repository): with the visual likelihood off
the log density equals HSGP4eDNA's `hsgp_2d_bathysp.stan` up to a constant; with
the eDNA likelihood off it equals the sum of the two single-species visual fits,
for both half-normal and hazard-rate detection. The spread over random draws is
at most 3e-11.

## Results on simulated data

Four fits on the same simulated dataset (`simulate_joint()`, seeds 1 and 2, 3
chains of 800 + 800). Cells are the field R² over all 200 stations and 1250
segments, with the mean 95% interval width of log density in brackets.

| Species | Seed | eDNA only | Visual only | Joint | Joint, conversion factor free |
|---|---|---|---|---|---|
| Hake | 1 | 0.987 (1.07) | not surveyed | 0.987 (1.09) | 0.987 (1.08) |
| Hake | 2 | 0.989 (1.27) | not surveyed | 0.988 (1.26) | |
| Humpback | 1 | 0.696 (2.08) | 0.897 (1.25) | **0.922 (1.15)** | 0.923 (1.15) |
| Humpback | 2 | 0.579 (2.40) | 0.869 (1.31) | **0.895 (1.26)** | |
| PWSD | 1 | 0.659 (2.54) | 0.533 (3.57) | **0.661** (2.63) | 0.665 (2.79) |
| PWSD | 2 | 0.751 (2.26) | 0.574 (3.47) | **0.809 (2.03)** | |

![Joint against single-source field recovery, seed 1: posterior mean against true log density, by species (rows) and fit (columns). This figure is from the repository's Phase 5 comparison.](figures/joint_compare_s1.png)

- **Hake is unchanged** by the joint fit, as it should be: it has no visual data.
- **Humpback improves** over the better single source in both seeds, in R² and in
  interval width.
- **PWSD** has the best R² in both seeds and narrower intervals in seed 2. In seed 1
  (37 sightings) its interval is slightly wider than eDNA-only: the visual data
  pull the along-coast length-scale from 435 to 199 km (truth 300), and a shorter
  length-scale means more uncertainty between data points.
- **Conversion factors are identified.** eDNA alone cannot identify them (only
  density times conversion enters its likelihood). With the visual data, and a
  deliberately off-centre prior, the joint model recovers them: humpback 5.33
  [5.09, 5.57] (truth 5.30) and PWSD 4.94 [4.46, 5.42] (truth 4.70), on the log
  scale.
- **One joint run did not converge** (seed 2, first attempt, max Rhat 1.66): one
  chain sat in a mode with very large junk noise. A rerun converged. This is
  an intermittent mixing risk, and the repository's advice is to check Rhat on
  $\theta_\text{junk}$ in every joint fit.
- **Cost.** eDNA-only and joint fits take 2 to 3 hours; visual-only 12 to 15 minutes.

> **[REVIEW] What is not yet shown here.**
> (1) Prior and posterior predictive checks for the joint model, in the style of
> Sections 3 and 4: not produced. (2) The joint model has not been fitted to the
> real MURI data. (3) The joint model's visual part does not yet include the
> **Beaufort covariate** that the PWSD real-data fit uses, and the real eDNA fit
> showed the diagnostics that need attention first (very large noise terms,
> length-scales at the floor, no Ct values, one flat water-column offset).
> Item (3) is the main reason I would treat a real-data joint fit as a
> demonstration rather than an analysis.

## What the real-data objects provide

`R/01_load_muri.r` writes an eDNA object and a visual object with the same
structure as the simulators' output, on a common domain and UTM origin, so the
same formatters apply to the joint model. Two data gaps carry over:
no qPCR Ct values (so the eDNA likelihood would use metabarcoding only, as in
Section 3.3) and no volumes filtered, which make the absolute conversion factor
meaningless until they are supplied. The visual data supply the absolute density.

# Limitations and open issues

## Data

- **No Ct or standard curve** for the real qPCR data, so the qPCR part of the eDNA
  model cannot be used. Volumes (filtered, aliquot, elution) are placeholders.
- **eDNA files have no date or cruise**, and the visual and eDNA surveys are
  different cruises. Both are treated as one snapshot of the same field.
- **Dilutions.** 1162 qPCR replicates were run at dilutions of 0.1 to 0.5 and were
  not corrected; metabarcoding replicates also carry template dilutions of 1 to
  100-fold that the model ignores.
- **Perpendicular distance missing** for 73 on-effort sightings, which are excluded.

## Models

- **g(0) = 1** in the visual model; there is no availability or perception bias
  correction.
- **Beaufort is the segment mean**, not the value at each sighting.
- **Unidentified animals** (baleen whales, unidentified dolphins) are not
  apportioned to species.
- **Water-column depth** is a known fixed offset in the eDNA model; on real data it
  is set to 1 because the real profiles are unknown.
- **Projection.** Effort to about -114.5 degrees longitude is far from UTM zone 10's
  central meridian, so distances there are distorted by a few percent.
- **Sampling cost.** The reference eDNA fit takes about 15 hours; real-data fits of
  both models are long, and the fits shown here use fewer chains and iterations
  than the validation runs.

# Reproducing this document

Run from the repository root (a CmdStan installation and the submodule are needed):

```bash
Rscript docs/explainer/05_fig_simulation.R     # simulation figures (needs the eDNA sim object)
Rscript docs/explainer/01_fit_edna_sim.R       # eDNA fit on simulated data
Rscript docs/explainer/06_fig_edna_sim.R       # eDNA figures
Rscript docs/explainer/02_fit_visual_sim.R     # visual fit on simulated data
SP=humpback Rscript docs/explainer/03_fit_visual_real.R   # real visual fit (also SP=pwsd)
SP=humpback M_MAX=800 TAG=real800 Rscript docs/explainer/03_fit_visual_real.R   # refit at twice the basis budget
Rscript docs/explainer/08_fig_key_choice.R     # detection-key comparison on the real distances
Rscript docs/explainer/04_fit_edna_real.R      # real eDNA fit (metabarcoding only)
Rscript docs/explainer/07_fig_visual.R         # visual figures
quarto render docs/explainer.md --to html --embed-resources
```

Fits are saved under `outputs/explainer/` (git-ignored); the figures are in
`docs/figures/`.

From HSGP4eDNA v1.4 the default basis budget is 800, so the first real visual
fit already uses about 800 functions. The "twice the budget" refit would be
`M_MAX=1600 TAG=real1600`. The real eDNA fit (`04`) also uses the default
budget, and will take roughly twice as long.

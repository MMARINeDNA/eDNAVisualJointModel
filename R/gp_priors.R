# =============================================================================
# gp_priors.R
#
# Reproducible, design-based priors for the HSGP latent field.
#
#   pr <- gp_prior_from_design(coords, domain_range, coord_scale)
#
# PRINCIPLE: the priors depend only on the survey DESIGN (sample locations and
# domain extent) and the compute budget, never on the observed outcomes
# (detections, counts, reads). So they are fixed before the data are seen,
# and two analysts with the same design get the same priors. Species surveyed
# with the same design get identical priors; the likelihood separates them.
#
# Length-scales get an inverse-gamma prior per axis, with `tail` mass below a
# lower bound ell_min and `tail` mass above an upper bound ell_max. This is the
# boundary-avoiding construction of Betancourt's "Robust Gaussian Processes in
# Stan" case study. The light left tail suppresses length-scales the design
# cannot resolve (the ell -> 0 overfitting collapse seen in HSGP4eDNA's basis
# sweep); the heavy right tail stays permissive up to the domain size.
#
#   ell_max_d = domain extent along axis d: longer features are
#               indistinguishable from the intercept.
#   ell_min_d = the larger of two floors:
#     design  - the median gap between distinct sample coordinates along d
#               (transect spacing, segment length; rarely binds for scattered
#               stations);
#     compute - the shortest scale an HSGP basis of at most M_max functions
#               can represent (hsgp_basis_rule(), Riutort-Mayol et al.). If
#               the design floor needs a bigger basis, both axes are scaled
#               up by a common factor (this minimises the largest inflation).
#
# The marginal SD gets a gamma prior with `tail` mass below sigma_bounds[1]
# (no mode at 0, avoiding the low-sigma trap of v3.2) and above
# sigma_bounds[2] (on the log-density scale; 3 means +/- 2 sigma spans ~400x).
#
# The returned basis `m` is hsgp_basis_rule() at ell_min, so the basis can
# represent every length-scale the prior supports.
#
# Requires hsgp_basis_rule() from HSGP4eDNA's R/functions.R (sourced by
# R/functions_visual.R).
# =============================================================================

# HSGP basis budget: the most basis functions (prod(m)) a single species field
# may use. This is a COMPUTE budget, not a statistical quantity. It sets the
# `compute` floor on length-scales, so lowering it raises the shortest scale
# the prior (and the model) can represent. Chosen 2026-10-07 so the floor stays
# below the 50 km simulation length-scale for every simulated species/data
# type (23-33 km); check runtime and floor_check() before changing it.
GP_BASIS_BUDGET <- 400L

# Gamma(shape, rate) with P(X < lo) = tail and P(X > hi) = tail.
gamma_from_quantiles <- function(lo, hi, tail = 0.01) {
  stopifnot(lo > 0, hi > lo)
  ratio <- function(a) qgamma(1 - tail, a) / qgamma(tail, a) - hi / lo
  a <- uniroot(ratio, c(0.05, 1e4), tol = 1e-10)$root
  c(shape = a, rate = qgamma(tail, a) / lo)
}

# Inverse-gamma(shape, scale) with P(X < lo) = tail and P(X > hi) = tail.
# If X ~ IG(a, b) then 1/X ~ Gamma(a, rate = b), so this is the gamma fit on
# 1/X with the bounds swapped.
inv_gamma_from_quantiles <- function(lo, hi, tail = 0.01) {
  g <- gamma_from_quantiles(1 / hi, 1 / lo, tail)
  c(shape = unname(g["shape"]), scale = unname(g["rate"]))
}

# Median gap between distinct sorted coordinate values (km).
axis_spacing <- function(x) {
  u <- sort(unique(round(x, 6)))
  if (length(u) < 2) return(0)
  median(diff(u))
}

# coords       n x D matrix of sample locations (km; one row per sample /
#              segment, detections or not)
# domain_range length-D domain extent (km)
# coord_scale  length-D normalisation half-range (km); the Stan models sample
#              length-scales in units of coord_scale
gp_prior_from_design <- function(coords, domain_range, coord_scale,
                                 tail = 0.01, M_max = GP_BASIS_BUDGET, c = 1.5,
                                 sigma_bounds = c(0.25, 3)) {
  coords <- as.matrix(coords)
  D <- ncol(coords)
  stopifnot(length(domain_range) == D, length(coord_scale) == D)

  ell_max <- domain_range
  floor_design <- apply(coords, 2, axis_spacing)
  ell_min <- floor_design

  # Compute floor: scale ell_min up (keeping its aspect) until the basis that
  # represents it fits within M_max.
  basis_at <- function(ell) hsgp_basis_rule(ell, coord_scale, c = c)
  f <- 1
  while (prod(basis_at(ell_min * f)) > M_max) f <- f * 1.01
  floor_compute <- ell_min * f
  ell_min <- floor_compute

  binding <- if (f > 1) "compute" else "design"
  ig <- t(vapply(seq_len(D), function(d)
    inv_gamma_from_quantiles(ell_min[d] / coord_scale[d],
                             ell_max[d] / coord_scale[d], tail), numeric(2)))
  sg <- gamma_from_quantiles(sigma_bounds[1], sigma_bounds[2], tail)

  list(
    # Stan data (length-scales in normalised units: ell / coord_scale)
    prior_gp_l_shape = unname(ig[, "shape"]),
    prior_gp_l_scale = unname(ig[, "scale"]),
    prior_gp_sigma_shape = unname(sg["shape"]),
    prior_gp_sigma_rate  = unname(sg["rate"]),
    HSGP_M = basis_at(ell_min),
    # Provenance
    table = data.frame(axis = seq_len(D), floor_design = floor_design,
                       floor_compute = floor_compute,
                       ell_min = ell_min, ell_max = ell_max, binding = binding,
                       m = basis_at(ell_min)),
    coord_scale = coord_scale,
    settings = list(tail = tail, M_max = M_max, c = c,
                    sigma_bounds = sigma_bounds, n_locations = nrow(unique(coords)))
  )
}

# -----------------------------------------------------------------------------
# floor_check(): truth-free check that the basis budget is not limiting a fit.
# For each axis, compares the posterior probability that ell < 1.25 * ell_min
# with the same probability under the prior. Posterior mass piling up against
# the floor (much more than the prior puts there) means the data want shorter
# length-scales than the budget allows -> raise GP_BASIS_BUDGET and refit.
#   ell_draws  draws x D matrix of length-scales (km)
#   prior      the gp_prior_from_design() result used for the fit
# -----------------------------------------------------------------------------
floor_check <- function(ell_draws, prior, near = 1.25, flag_ratio = 3) {
  ell_draws <- as.matrix(ell_draws)
  t <- prior$table
  out <- lapply(seq_len(ncol(ell_draws)), function(d) {
    thr <- near * t$ell_min[d]
    # prior P(ell < thr): ell_norm ~ IG(shape, scale) with ell = ell_norm * coord_scale
    p_prior <- 1 - pgamma(prior$coord_scale[d] / thr, prior$prior_gp_l_shape[d],
                          rate = prior$prior_gp_l_scale[d])
    p_post  <- mean(ell_draws[, d] < thr)
    data.frame(axis = d, ell_min = t$ell_min[d], post_median = median(ell_draws[, d]),
               p_post_near_floor = p_post, p_prior_near_floor = p_prior,
               flag = p_post > max(flag_ratio * p_prior, 0.05))
  })
  do.call(rbind, out)
}

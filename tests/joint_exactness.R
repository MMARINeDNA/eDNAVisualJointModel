# =============================================================================
# tests/joint_exactness.R
#
# Exactness checks for stan/hsgp_joint.stan's per-species visual options (the
# sea-state and visibility covariates and negative-binomial segment counts), on log densities
# at random parameter values (cmdstanr model methods):
#
#   1. Negative-binomial counts and visibility off for every species (eDNA +
#      visual on, Beaufort off): the log density equals the joint model at
#      BASE_REF (before these options) exactly.
#   2. eDNA off; PWSD with the Beaufort and visibility covariates and NB
#      counts, humpback with none; data with sea state and visibility: the log
#      density equals humpback's and PWSD's single-species visual models
#      (stan/hsgp_visual.stan; PWSD with use_bft, use_vis and use_nb) up to a
#      constant. Only the cetaceans' field and
#      visual parameters vary between points; the visual models get the joint
#      model's field priors and coordinates.
#
#   Rscript tests/joint_exactness.R
# Env: BASE_REF (git ref of the joint model before the NB option; default
#      origin/main), N_POINTS (default 10), TOL (default 1e-8)
# Exit status non-zero on failure. Models are compiled in a temporary directory.
# =============================================================================

suppressMessages({library(cmdstanr); library(posterior)})
source("R/fit_joint.r")
source("R/fit_visual.r")
options(cmdstanr_warn_inits = FALSE)

BASE_REF <- Sys.getenv("BASE_REF", "origin/main")
N_POINTS <- as.integer(Sys.getenv("N_POINTS", "10"))
TOL      <- as.numeric(Sys.getenv("TOL", "1e-8"))

# ---- Data (simulated before any model methods are loaded) -------------------
sim0 <- simulate_joint(seed = 2L, n_stations = 120L, n_transects = 12L)
sim1 <- simulate_joint(seed = 2L, n_stations = 120L, n_transects = 12L, sea_state = TRUE,
                       visibility = TRUE)
d_new <- format_stan_data_joint(sim0, M_max = 64L)
stopifnot(all(d_new$use_bft_covar == 0L), d_new$n_bft == 1L, all(d_new$use_nb == 0L),
          all(d_new$use_vis_covar == 0L))
d_old <- d_new[setdiff(names(d_new), c("use_nb", "prior_nb_inv_sqrt_sd", "use_vis_covar", "vis_level",
                                       "beta_vis_prior_mean", "beta_vis_prior_sd"))]
d_vis <- format_stan_data_joint(sim1, use_edna = 0L, bft_species = "pwsd", vis_species = "pwsd",
                                nb_species = "pwsd", M_max = 64L)
stopifnot(all(d_vis$use_bft_covar == c(0L, 1L)), all(d_vis$use_nb == c(0L, 1L)),
          all(d_vis$use_vis_covar == c(0L, 1L)))
HSGP_M <- apply(d_vis$INDICES, 2, max)
vd <- lapply(c(humpback = FALSE, pwsd = TRUE), function(ub) NULL)
vd$humpback <- format_stan_data_visual(sim1$visual, "humpback", HSGP_M = HSGP_M, use_bft = FALSE,
                                       use_vis = FALSE, use_nb = 0L)
vd$pwsd     <- format_stan_data_visual(sim1$visual, "pwsd",     HSGP_M = HSGP_M, use_bft = TRUE,
                                       use_vis = TRUE, use_nb = 1L)
for (k in seq_along(vd)) {                     # the joint model's field priors and coordinates
  s <- d_vis$lt_species_idx[k]
  vd[[k]] <- modifyList(vd[[k]], list(
    coords = d_vis$seg_coords, coord_scale = d_vis$coord_scale, L_hsgp = d_vis$L_hsgp,
    B_bathy = d_vis$B_bathy_seg, prior_beta_bathy_sig = d_vis$prior_beta_bathy_sig,
    prior_gp_sigma_shape = d_vis$prior_gp_sigma_shape, prior_gp_sigma_rate = d_vis$prior_gp_sigma_rate,
    prior_gp_l_shape = d_vis$prior_gp_l_shape, prior_gp_l_scale = d_vis$prior_gp_l_scale,
    prior_mu_sp_mu = d_vis$prior_mu_sp_mu[s], prior_mu_sp_sig = d_vis$prior_mu_sp_sig[s],
    prior_nb_inv_sqrt_sd = d_vis$prior_nb_inv_sqrt_sd))
  stopifnot(identical(vd[[k]]$INDICES, d_vis$INDICES))
}
stopifnot(identical(vd$pwsd$seg_bft, as.integer(d_vis$seg_bft)),
          isTRUE(all.equal(as.numeric(vd$pwsd$bft_level), as.numeric(d_vis$bft_level))),
          isTRUE(all.equal(as.numeric(vd$pwsd$vis_level), as.numeric(d_vis$vis_level))))

# ---- Models with methods, in a temporary directory --------------------------
tmp <- file.path(tempdir(), "joint_exactness"); dir.create(file.path(tmp, "old"), recursive = TRUE)
old_src <- system2("git", c("show", sprintf("%s:stan/hsgp_joint.stan", BASE_REF)), stdout = TRUE)
writeLines(old_src, file.path(tmp, "old", "hsgp_joint.stan"))
invisible(file.copy(JOINT_STAN, file.path(tmp, "hsgp_joint.stan"), overwrite = TRUE))
invisible(file.copy(VISUAL_STAN, file.path(tmp, "hsgp_visual.stan"), overwrite = TRUE))
inc <- normalizePath(JOINT_INCLUDE)
mk <- function(f) cmdstan_model(f, include_paths = inc, compile_model_methods = TRUE,
                                force_recompile = TRUE, quiet = TRUE)
# A fit object to carry the model methods: one fixed-parameter draw at the inits
carrier <- function(mod, data, init) {
  f <- mod$sample(data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
                  init = list(init), refresh = 0, show_messages = FALSE)
  f$init_model_methods(); f
}
lp <- function(f, pars) f$log_prob(f$unconstrain_variables(pars), jacobian = TRUE)
drop_empty <- function(p) p[lengths(p) > 0]

set.seed(7)
jitter <- function(p, sp_rows, lt_rows) {
  # perturb the field of species sp_rows and the visual parameters of lt_rows
  p$mu_sp[sp_rows]       <- p$mu_sp[sp_rows] + rnorm(length(sp_rows), 0, 0.3)
  p$gp_sigma[sp_rows]    <- p$gp_sigma[sp_rows] * exp(rnorm(length(sp_rows), 0, 0.2))
  p$gp_l_raw[sp_rows, ]  <- p$gp_l_raw[sp_rows, ] * exp(rnorm(2 * length(sp_rows), 0, 0.2))
  p$z_beta[sp_rows, ]    <- p$z_beta[sp_rows, ] + rnorm(length(p$z_beta[sp_rows, ]), 0, 0.3)
  p$beta_bathy[sp_rows, ] <- p$beta_bathy[sp_rows, ] + rnorm(length(p$beta_bathy[sp_rows, ]), 0, 0.3)
  p$log_sigma[lt_rows]   <- p$log_sigma[lt_rows] + rnorm(length(lt_rows), 0, 0.1)
  p$beta_size[lt_rows]   <- p$beta_size[lt_rows] + rnorm(length(lt_rows), 0, 0.005)
  if (length(p$beta_bft)) p$beta_bft <- p$beta_bft + rnorm(length(p$beta_bft), 0, 0.05)
  if (length(p$beta_vis)) p$beta_vis <- p$beta_vis + rnorm(length(p$beta_vis), 0, 0.05)
  if (length(p$inv_sqrt_phi_nb)) p$inv_sqrt_phi_nb <- array(exp(rnorm(length(p$inv_sqrt_phi_nb), log(0.7), 0.5)))
  p$mu_s[lt_rows]        <- p$mu_s[lt_rows] * exp(rnorm(length(lt_rows), 0, 0.1))
  p$phi_s[lt_rows]       <- p$phi_s[lt_rows] * exp(rnorm(length(lt_rows), 0, 0.1))
  p$mu_log_s[lt_rows]    <- p$mu_log_s[lt_rows] + rnorm(length(lt_rows), 0, 0.1)
  p$sigma_log_s[lt_rows] <- p$sigma_log_s[lt_rows] * exp(rnorm(length(lt_rows), 0, 0.1))
  p
}
fail <- FALSE

# ---- 1. NB and visibility off: identical to the joint model at BASE_REF ------
m_new <- mk(file.path(tmp, "hsgp_joint.stan")); m_old <- mk(file.path(tmp, "old", "hsgp_joint.stan"))
p0 <- drop_empty(joint_init(d_new)())
f_new <- carrier(m_new, d_new, p0); f_old <- carrier(m_old, d_old, p0)
d1 <- vapply(seq_len(N_POINTS), function(k) {
  p <- jitter(p0, 1:d_new$S, seq_len(d_new$L))
  lp(f_new, p) - lp(f_old, p)
}, numeric(1))
cat(sprintf("1. NB and visibility off vs %s: max |difference| over %d points = %.2e\n",
            BASE_REF, N_POINTS, max(abs(d1))))
if (max(abs(d1)) > TOL) { cat("   FAIL\n"); fail <- TRUE }

# ---- 2. eDNA off, PWSD with the covariate: the two visual models -----------
m_vis <- mk(file.path(tmp, "hsgp_visual.stan"))
pj <- drop_empty(joint_init(d_vis)())
to_vis <- function(p, k) {
  s <- d_vis$lt_species_idx[k]
  drop_empty(list(mu_sp = p$mu_sp[s], gp_sigma = array(p$gp_sigma[s]),
       gp_l_raw = matrix(p$gp_l_raw[s, ], 1), z_beta = matrix(p$z_beta[s, ], 1),
       beta_bathy = array(p$beta_bathy[s, ]),
       log_sigma = p$log_sigma[k], beta_size = p$beta_size[k],
       beta_bft = if (d_vis$use_bft_covar[k] == 1L) array(p$beta_bft[sum(d_vis$use_bft_covar[1:k])]),
       beta_vis = if (d_vis$use_vis_covar[k] == 1L) array(p$beta_vis[sum(d_vis$use_vis_covar[1:k])]),
       inv_sqrt_phi_nb = if (d_vis$use_nb[k] == 1L) array(p$inv_sqrt_phi_nb[sum(d_vis$use_nb[1:k])]),
       mu_s = p$mu_s[k], phi_s = p$phi_s[k], mu_log_s = p$mu_log_s[k], sigma_log_s = p$sigma_log_s[k]))
}
f_j  <- carrier(m_new, d_vis, pj)
f_vh <- carrier(m_vis, vd$humpback, to_vis(pj, 1)); f_vp <- carrier(m_vis, vd$pwsd, to_vis(pj, 2))
d2 <- vapply(seq_len(N_POINTS), function(k) {
  p <- jitter(pj, d_vis$lt_species_idx, seq_len(d_vis$L))
  lp(f_j, p) - lp(f_vh, to_vis(p, 1)) - lp(f_vp, to_vis(p, 2))
}, numeric(1))
cat(sprintf("2. eDNA off, PWSD with Beaufort, visibility and NB vs the visual models: spread of the difference over %d points = %.2e\n",
            N_POINTS, diff(range(d2))))
if (diff(range(d2)) > TOL) { cat("   FAIL\n"); fail <- TRUE }

if (fail) { cat("\nJOINT EXACTNESS FAILED\n"); quit(status = 1) }
cat("\nJOINT EXACTNESS PASSED\n")

# =============================================================================
# tests/joint_exactness.R
#
# Exactness checks for stan/hsgp_joint.stan's options (the visual sea-state and
# visibility covariates and negative-binomial segment counts; the eDNA dilution
# offsets, per-bottle volume, estimated qPCR curves, Ct-SD bound and
# sampling-depth smooth ported from HSGP4eDNA v1.8-v1.12), on log densities at
# random parameter values (cmdstanr model methods):
#
#   1. Every option off (eDNA + visual on, Beaufort off, no dilutions, one
#      volume, the fixed qPCR curve, no depth smooth): the log density equals
#      the joint model at BASE_REF (before these options) exactly.
#   2. eDNA off; PWSD with the Beaufort and visibility covariates and NB
#      counts, humpback with none; data with sea state and visibility: the log
#      density equals humpback's and PWSD's single-species visual models
#      (stan/hsgp_visual.stan; PWSD with use_bft, use_vis and use_nb) up to a
#      constant. Only the cetaceans' field and
#      visual parameters vary between points; the visual models get the joint
#      model's field priors and coordinates.
#   3. Visual off; the eDNA options all on (dilutions, per-bottle volumes,
#      estimated per-plate curves from standards, the conditional Ct SD, the
#      depth smooth): the log density equals HSGP4eDNA's own
#      stan/hsgp_2d_bathysp.stan on the same data up to a constant. The eDNA
#      parameters vary between points; the visual parameters are fixed.
#
#   Rscript tests/joint_exactness.R
# Env: BASE_REF (git ref of the joint model before the NB and visibility
#      options; default 0a60f03, the commit before that work), N_POINTS
#      (default 10), TOL (default 1e-8)
# Exit status non-zero on failure. Models are compiled in a temporary directory.
# =============================================================================

suppressMessages({library(cmdstanr); library(posterior)})
source("R/fit_joint.r")
source("R/fit_visual.r")
options(cmdstanr_warn_inits = FALSE)

BASE_REF <- Sys.getenv("BASE_REF", "0a60f03")
N_POINTS <- as.integer(Sys.getenv("N_POINTS", "10"))
TOL      <- as.numeric(Sys.getenv("TOL", "1e-8"))

# ---- Data (simulated before any model methods are loaded) -------------------
sim0 <- simulate_joint(seed = 2L, n_stations = 120L, n_transects = 12L)
sim1 <- simulate_joint(seed = 2L, n_stations = 120L, n_transects = 12L, sea_state = TRUE,
                       visibility = TRUE)
d_new <- format_stan_data_joint(sim0, M_max = 64L)
stopifnot(all(d_new$use_bft_covar == 0L), d_new$n_bft == 1L, all(d_new$use_nb == 0L),
          all(d_new$use_vis_covar == 0L))
stopifnot(all(d_new$qpcr_log_dilution == 0), all(d_new$mb_log_dilution == 0),
          all(d_new$log_vol_bottle == d_new$log_vol_filtered), d_new$qpcr_curve_model == 0L,
          d_new$ct_sd_cond == 0L, d_new$K_zdepth == 0L)
d_old <- d_new[setdiff(names(d_new), c("use_nb", "prior_nb_inv_sqrt_sd", "use_vis_covar", "vis_level",
                                       "beta_vis_prior_mean", "beta_vis_prior_sd",
                                       "qpcr_log_dilution", "mb_log_dilution", "log_vol_bottle",
                                       "qpcr_curve_model", "N_plate", "qpcr_plate_idx", "N_std",
                                       "std_plate_idx", "std_detect", "std_ct", "std_log_copies",
                                       "prior_alpha_ct", "prior_beta_ct", "beta_ct_min",
                                       "prior_gamma0_ct", "prior_gamma1_ct", "prior_sigma0_ct_sig",
                                       "prior_logit_kappa", "ct_sd_cond", "K_zdepth", "B_zdepth",
                                       "prior_beta_zdepth_sig"))]
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

# eDNA options on (check 3): five sampling depths, and dilutions, volumes,
# plates and standards added to the eDNA object (values need not be realistic)
set.seed(3)
sim3 <- simulate_joint(seed = 3L, n_stations = 60L, n_transects = 8L, sample_depths = c(0, 50, 150, 300, 500))
o3 <- sim3$edna$observed; nq <- length(o3$qpcr_sample_idx); nm <- length(o3$mb_sample_idx)
o3$qpcr_log_dilution <- log(sample(c(1, 0.2), nq, replace = TRUE))
o3$mb_log_dilution   <- -log(sample(c(1, 5, 10), nm, replace = TRUE))
o3$qpcr_plate_idx    <- sample(1:4, nq, replace = TRUE)
std <- expand.grid(rep = 1:3, K = 10^(0:4), plate = 1:5)
std_mu <- 38 - 1.45 * log(2 * std$K)
std_det <- as.integer(runif(nrow(std)) < 1 - exp(-0.9 * 2 * std$K))
o3$qpcr_std_plate_idx <- std$plate; o3$qpcr_std_copies_ul <- std$K; o3$qpcr_std_detect <- std_det
o3$qpcr_std_ct <- ifelse(std_det == 1L, std_mu + rnorm(nrow(std), 0, 0.3), NA)
sim3$edna$observed <- o3
sim3$edna$design$samples$vol_filtered <- runif(sim3$edna$meta$N, 1.2, 2.5)
d3 <- format_stan_data_joint(sim3, use_visual = 0L, qpcr_curve = "estimated", zdepth = "smooth",
                             ct_sd = "conditional", M_max = 64L)
stopifnot(d3$qpcr_curve_model == 1L, d3$ct_sd_cond == 1L, d3$K_zdepth > 0L, d3$N_plate == 5L,
          any(d3$qpcr_log_dilution != 0), any(d3$mb_log_dilution != 0), sd(d3$log_vol_bottle) > 0)
# HSGP4eDNA's data for the same eDNA object and options, with the joint model's
# field (coordinates, basis, priors, bottom-depth basis) and noise settings
e3 <- format_stan_data_gp2d_bathysp(sim3$edna, prediction_n = 2L, mb_phi = "empirical",
                                    qpcr_curve = "estimated", zdepth = "smooth")$stan_data
e3 <- modifyList(e3, d3[c("coords", "INDICES", "M", "L_hsgp", "coord_scale", "prior_gp_l_shape",
                          "prior_gp_l_scale", "prior_gp_sigma_shape", "prior_gp_sigma_rate",
                          "B_bathy", "K_bathy", "prior_mu_sp_mu", "prior_mu_sp_sig",
                          "prior_beta_bathy_sig", "log_zsample_effect", "use_bottle_noise",
                          "bottle_noise_dist", "bottle_theta_shared", "prior_cv_bottle_sig",
                          "use_junk_noise", "junk_noise_dist", "prior_cv_junk_sig", "ct_sd_cond")])
e3$m_hsgp <- as.array(apply(d3$INDICES, 2, max)); e3$prior_gp_l_family <- 1L
e3$pred_coords <- d3$coords[seq_len(e3$N_pred), , drop = FALSE]
e3$B_bathy_pred <- d3$B_bathy[seq_len(e3$N_pred), , drop = FALSE]
shared3 <- c("qpcr_sample_idx", "qpcr_detect", "qpcr_ct", "mb_sample_idx", "mb_reads", "mb_total",
             "qpcr_log_dilution", "mb_log_dilution", "log_vol_bottle", "qpcr_plate_idx", "std_plate_idx",
             "std_detect", "std_ct", "std_log_copies", "B_zdepth", "log_conv_factor", "log_vol_filtered")
stopifnot(all(vapply(shared3, function(k) isTRUE(all.equal(as.vector(as.matrix(e3[[k]])),
                                                           as.vector(as.matrix(d3[[k]])))), logical(1))))

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

# ---- 1. Every option off: identical to the joint model at BASE_REF ----------
m_new <- mk(file.path(tmp, "hsgp_joint.stan")); m_old <- mk(file.path(tmp, "old", "hsgp_joint.stan"))
p0 <- drop_empty(joint_init(d_new)())
f_new <- carrier(m_new, d_new, p0); f_old <- carrier(m_old, d_old, p0)
d1 <- vapply(seq_len(N_POINTS), function(k) {
  p <- jitter(p0, 1:d_new$S, seq_len(d_new$L))
  lp(f_new, p) - lp(f_old, p)
}, numeric(1))
cat(sprintf("1. every option off vs %s: max |difference| over %d points = %.2e\n",
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

# ---- 3. Visual off, eDNA options on: HSGP4eDNA's eDNA model ------------------
tmp_e <- file.path(tmp, "edna"); dir.create(tmp_e)
invisible(file.copy(file.path(HSGP4EDNA_DIR, "stan", "hsgp_2d_bathysp.stan"),
                    file.path(tmp_e, "hsgp_2d_bathysp.stan"), overwrite = TRUE))
m_edna <- cmdstan_model(file.path(tmp_e, "hsgp_2d_bathysp.stan"),
                        include_paths = normalizePath(file.path(HSGP4EDNA_DIR, "stan", "include")),
                        compile_model_methods = TRUE, force_recompile = TRUE, quiet = TRUE)
p3 <- drop_empty(joint_init(d3)())
edna_pars <- c("mu_sp", "gp_sigma", "gp_l_raw", "z_beta", "beta_bathy", "beta0_phi", "gamma0_phi",
               "gamma1_phi", "g_raw", "cv_bottle", "g_junk_raw", "cv_junk", "alpha_ct_plate",
               "beta_ct_plate", "gamma0_ct_est", "gamma1_ct_est", "sigma0_ct_est", "logit_kappa_est",
               "beta_zdepth")
to_edna <- function(p) drop_empty(p[intersect(edna_pars, names(p))])
jitter3 <- function(p) {
  p <- jitter(p, 1:d3$S, integer(0))
  p$beta0_phi <- p$beta0_phi + rnorm(1, 0, 0.3); p$gamma0_phi <- p$gamma0_phi * exp(rnorm(1, 0, 0.2))
  p$gamma1_phi <- p$gamma1_phi * exp(rnorm(1, 0, 0.2))
  if (length(p$g_raw)) p$g_raw <- p$g_raw + rnorm(length(p$g_raw), 0, 0.2)
  if (length(p$g_junk_raw)) p$g_junk_raw <- p$g_junk_raw + rnorm(length(p$g_junk_raw), 0, 0.2)
  p$alpha_ct_plate <- p$alpha_ct_plate + rnorm(length(p$alpha_ct_plate), 0, 0.3)
  p$beta_ct_plate  <- pmax(p$beta_ct_plate + rnorm(length(p$beta_ct_plate), 0, 0.02), d3$beta_ct_min + 1e-3)
  p$gamma0_ct_est  <- p$gamma0_ct_est + rnorm(1, 0, 0.2); p$gamma1_ct_est <- p$gamma1_ct_est + rnorm(1, 0, 0.1)
  p$sigma0_ct_est  <- p$sigma0_ct_est * exp(rnorm(1, 0, 0.2))
  p$logit_kappa_est <- p$logit_kappa_est + rnorm(1, 0, 0.3)
  p$beta_zdepth <- p$beta_zdepth + rnorm(length(p$beta_zdepth), 0, 0.3)
  p
}
f_j3 <- carrier(m_new, d3, p3); f_e3 <- carrier(m_edna, e3, to_edna(p3))
d3d <- vapply(seq_len(N_POINTS), function(k) {
  p <- jitter3(p3)
  lp(f_j3, p) - lp(f_e3, to_edna(p))
}, numeric(1))
cat(sprintf("3. visual off, eDNA options on vs HSGP4eDNA's hsgp_2d_bathysp.stan: spread of the difference over %d points = %.2e\n",
            N_POINTS, diff(range(d3d))))
if (diff(range(d3d)) > TOL) { cat("   FAIL\n"); fail <- TRUE }

if (fail) { cat("\nJOINT EXACTNESS FAILED\n"); quit(status = 1) }
cat("\nJOINT EXACTNESS PASSED\n")

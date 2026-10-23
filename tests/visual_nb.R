# =============================================================================
# tests/visual_nb.R
#
# Checks the negative-binomial segment counts of stan/hsgp_visual.stan
# (use_nb), on log densities at random parameter values (cmdstanr model
# methods) and on one recovery fit:
#
#   1. use_nb = 0: the log density equals the visual model at BASE_REF
#      (before the option) exactly.
#   2. use_nb = 1: the per-segment log-likelihood (log_lik_seg) equals R's
#      dnbinom(size = phi_nb, mu) at the model's expected detected groups.
#   3. Recovery (RECOVERY=1, the default): NB groups per segment (size 1)
#      simulated on a small design; the true phi_nb lies in the 99% interval.
#      The field's length-scales (150 x 400 km) are long enough for the small
#      basis: a field the basis cannot represent reads as extra overdispersion
#      and pulls phi_nb down.
#
#   Rscript tests/visual_nb.R
# Env: BASE_REF (the visual model before the NB and visibility options;
#      default 0a60f03), N_POINTS (default 10), TOL (1e-8),
#      RECOVERY (1), CHAINS (4)
# Exit status non-zero on failure. Models are compiled in a temporary directory.
# =============================================================================

suppressMessages({library(cmdstanr); library(posterior)})
source("R/fit_visual.r")
options(cmdstanr_warn_inits = FALSE)

BASE_REF <- Sys.getenv("BASE_REF", "0a60f03")
N_POINTS <- as.integer(Sys.getenv("N_POINTS", "10"))
TOL      <- as.numeric(Sys.getenv("TOL", "1e-8"))

# ---- Data -------------------------------------------------------------------
gp <- default_gp_params(); gp$humpback$mu <- log(0.03)
gp$humpback$lx <- 150; gp$humpback$ly <- 400
sim <- simulate_visual(seed = 3L, n_transects = 12L, gp_params = gp,
                       nb_size = c(humpback = 1, pwsd = 1))
d0 <- format_stan_data_visual(sim, "humpback", M_max = 64L, use_nb = 0L)
d1 <- modifyList(d0, list(use_nb = 1L))
d_old <- d0[setdiff(names(d0), c("use_nb", "prior_nb_inv_sqrt_sd", "use_vis_covar", "vis_level",
                                 "beta_vis_prior_mean", "beta_vis_prior_sd"))]

# ---- Models with methods, in a temporary directory --------------------------
tmp <- file.path(tempdir(), "visual_nb"); dir.create(file.path(tmp, "old", "include"), recursive = TRUE)
show <- function(path) system2("git", c("show", sprintf("%s:%s", BASE_REF, path)), stdout = TRUE)
writeLines(show("stan/hsgp_visual.stan"), file.path(tmp, "old", "hsgp_visual.stan"))
writeLines(show("stan/include/visual_functions.stan"), file.path(tmp, "old", "include", "visual_functions.stan"))
invisible(file.copy(VISUAL_STAN, file.path(tmp, "hsgp_visual.stan"), overwrite = TRUE))
hsgp_inc <- normalizePath(file.path(HSGP4EDNA_DIR, "stan", "include"))
mk <- function(f, inc) cmdstan_model(f, include_paths = inc, compile_model_methods = TRUE,
                                     force_recompile = TRUE, quiet = TRUE)
m_new <- mk(file.path(tmp, "hsgp_visual.stan"), c(normalizePath("stan/include"), hsgp_inc))
m_old <- mk(file.path(tmp, "old", "hsgp_visual.stan"), c(file.path(tmp, "old", "include"), hsgp_inc))
carrier <- function(mod, data, init) {
  f <- mod$sample(data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
                  init = list(init), refresh = 0, show_messages = FALSE)
  f$init_model_methods(); f
}
lp <- function(f, pars) f$log_prob(f$unconstrain_variables(pars), jacobian = TRUE)
drop_empty <- function(p) p[lengths(p) > 0]
set.seed(7)
jitter <- function(p) {
  p$mu_sp <- p$mu_sp + rnorm(1, 0, 0.3); p$gp_sigma <- p$gp_sigma * exp(rnorm(1, 0, 0.2))
  p$gp_l_raw <- p$gp_l_raw * exp(rnorm(2, 0, 0.2)); p$z_beta <- p$z_beta + rnorm(length(p$z_beta), 0, 0.3)
  p$beta_bathy <- p$beta_bathy + rnorm(length(p$beta_bathy), 0, 0.3)
  p$log_sigma <- p$log_sigma + rnorm(1, 0, 0.1); p$mu_s <- p$mu_s * exp(rnorm(1, 0, 0.1))
  if (length(p$inv_sqrt_phi_nb)) p$inv_sqrt_phi_nb <- array(exp(rnorm(1, log(0.7), 0.5)))
  p
}
fail <- FALSE

# ---- 1. use_nb = 0 equals the model at BASE_REF ------------------------------
p0 <- drop_empty(visual_init(d0)())
f_new <- carrier(m_new, d0, p0); f_old <- carrier(m_old, d_old, p0)
diff1 <- vapply(seq_len(N_POINTS), function(k) { p <- jitter(p0); lp(f_new, p) - lp(f_old, p) }, numeric(1))
cat(sprintf("1. use_nb = 0 vs %s: max |difference| over %d points = %.2e\n", BASE_REF, N_POINTS, max(abs(diff1))))
if (max(abs(diff1)) > TOL) { cat("   FAIL\n"); fail <- TRUE }

# ---- 2. use_nb = 1: log_lik_seg equals dnbinom -------------------------------
p1 <- drop_empty(visual_init(d1)())
f1 <- carrier(m_new, d1, p1)
diff2 <- vapply(seq_len(N_POINTS), function(k) {
  p <- jitter(p1)
  v <- f1$constrain_variables(f1$unconstrain_variables(p), transformed_parameters = TRUE,
                              generated_quantities = TRUE)
  mu <- exp(v$log_lambda - log(v$mean_group_size) + log(2 * d1$seg_l * v$esw_bft[d1$seg_bft]))
  max(abs(v$log_lik_seg - dnbinom(d1$seg_count, size = v$phi_nb, mu = mu, log = TRUE)))
}, numeric(1))
cat(sprintf("2. use_nb = 1: max |log_lik_seg - dnbinom| over %d points = %.2e\n", N_POINTS, max(diff2)))
if (max(diff2) > 1e-6) { cat("   FAIL\n"); fail <- TRUE }

# ---- 3. Recovery of phi_nb ---------------------------------------------------
if (as.integer(Sys.getenv("RECOVERY", "1")) == 1L) {
  fit <- m_new$sample(data = d1, chains = as.integer(Sys.getenv("CHAINS", "4")),
                      parallel_chains = as.integer(Sys.getenv("CHAINS", "4")),
                      iter_warmup = 500, iter_sampling = 500, adapt_delta = 0.9, max_treedepth = 12,
                      seed = 11, init = visual_init(d1), refresh = 0, show_messages = FALSE)
  q <- quantile(fit$draws("phi_nb", format = "draws_matrix"), c(0.005, 0.5, 0.995))
  cnt <- d1$seg_count
  cat(sprintf("3. recovery: %d segments, %d groups (var/mean %.2f); phi_nb median %.2f, 99%% [%.2f, %.2f] (truth 1); divergences %d\n",
              d1$n_seg, sum(cnt), var(cnt) / mean(cnt), q[2], q[1], q[3],
              sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)))
  if (q[1] > 1 || q[3] < 1) { cat("   FAIL\n"); fail <- TRUE }
}

if (fail) { cat("\nVISUAL NB CHECKS FAILED\n"); quit(status = 1) }
cat("\nVISUAL NB CHECKS PASSED\n")

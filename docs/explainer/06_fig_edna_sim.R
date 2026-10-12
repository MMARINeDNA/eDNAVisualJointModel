# Figures for the eDNA-only model on simulated data (bathysp_depthspref, 80 stations).
# Needs outputs/explainer/edna_sim{,_standata,_fit}.rds from 01_fit_edna_sim.R.
source("docs/explainer/00_common.R")
source(file.path(HSGP4EDNA_DIR, "R/functions.R"))
sim <- readRDS(file.path(EXPL_OUT, "edna_sim.rds"))
sd  <- readRDS(file.path(EXPL_OUT, "edna_sim_standata.rds"))
fit <- readRDS(file.path(EXPL_OUT, "edna_sim_fit.rds"))
S <- 3; N <- sd$N; spn <- names(SP_LAB)
samp <- sim$design$samples
lvf <- log(sd$vol_aliquot / 100)

draws <- as_draws_matrix(fit$draws())
pick <- function(prefix) draws[, startsWith(colnames(draws), prefix), drop = FALSE]
mat_idx <- function(v) { m <- str_match(colnames(v), "\\[(\\d+),(\\d+)\\]"); cbind(as.integer(m[, 2]), as.integer(m[, 3])) }

# ---- diagnostics summary -----------------------------------------------------------
diag <- fit$diagnostic_summary(quiet = TRUE)
hyper_vars <- c(grep("^gp_l\\[|^gp_sigma|^mu_sp|^theta", colnames(draws), value = TRUE))
sm <- summarise_draws(draws[, hyper_vars], rhat, ess_bulk, ess_tail)
saveRDS(list(divergences = sum(diag$num_divergent), treedepth = sum(diag$num_max_treedepth),
             ebfmi = diag$ebfmi, max_rhat = max(sm$rhat, na.rm = TRUE),
             min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE), M = sd$M, m_hsgp = sd$m_hsgp, N = sd$N,
             N_qpcr = sd$N_qpcr_long, N_mb = sd$N_mb_long,
             runtime_min = max(fit$time()$chains$total) / 60),
        file.path(EXPL_OUT, "edna_sim_diag.rds"))

# ---- hyperparameter recovery vs prior -----------------------------------------------
gpt <- sim$truth$gp_params
cs  <- sd$coord_scale
rec <- map_dfr(seq_len(S), function(s) {
  get <- function(v) as.numeric(draws[, v])
  tibble(species = spn[s], param = c("lx (km)", "ly (km)", "gp_sigma", "mu"),
         truth = c(gpt[[s]]$lx, gpt[[s]]$ly, gpt[[s]]$sigma, gpt[[s]]$mu),
         post = list(get(sprintf("gp_l[%d,1]", s)), get(sprintf("gp_l[%d,2]", s)),
                     get(sprintf("gp_sigma[%d]", s)), get(sprintf("mu_sp[%d]", s)))) |>
    mutate(med = map_dbl(post, median), lo = map_dbl(post, quantile, .025), hi = map_dbl(post, quantile, .975))
})
qig <- function(p, shape, scale) 1 / qgamma(1 - p, shape, rate = scale)
prior_rng <- tibble(param = c("lx (km)", "ly (km)", "gp_sigma", "mu"),
  lo = c(qig(.025, sd$prior_gp_l_shape[1], sd$prior_gp_l_scale[1]) * cs[1],
         qig(.025, sd$prior_gp_l_shape[2], sd$prior_gp_l_scale[2]) * cs[2],
         qgamma(.025, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate),
         qnorm(.025, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)),
  hi = c(qig(.975, sd$prior_gp_l_shape[1], sd$prior_gp_l_scale[1]) * cs[1],
         qig(.975, sd$prior_gp_l_shape[2], sd$prior_gp_l_scale[2]) * cs[2],
         qgamma(.975, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate),
         qnorm(.975, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)))
rec$param <- factor(rec$param, c("lx (km)", "ly (km)", "gp_sigma", "mu")); prior_rng$param <- factor(prior_rng$param, levels(rec$param))
rec$species <- factor(rec$species, spn, SP_LAB)
p1 <- ggplot(rec, aes(species)) +
  geom_rect(data = prior_rng, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey88") +
  geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac") +
  geom_point(aes(y = truth), colour = "#b2182b", shape = 4, size = 3, stroke = 1.3) +
  facet_wrap(~param, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "GP hyperparameters: posterior (blue, 95% CI) vs truth (red x)",
       subtitle = "grey band: central 95% of the design-based prior") + theme_expl() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_fig(p1, "edna_sim_hyper", 12, 4)
rec |> mutate(covered = truth >= lo & truth <= hi) |> select(species, param, truth, med, lo, hi, covered) |>
  saveRDS(file.path(EXPL_OUT, "edna_sim_hyper_table.rds"))

# ---- latent field recovery ---------------------------------------------------------------
ll <- pick("log_lambda["); ll <- ll[, !startsWith(colnames(ll), "log_lambda_edna") & !startsWith(colnames(ll), "log_lambda_pred"), drop = FALSE]
ix <- mat_idx(ll)
pm <- matrix(NA, N, S); pq <- array(NA, c(N, S, 2))
for (j in seq_len(ncol(ll))) { pm[ix[j, 1], ix[j, 2]] <- mean(ll[, j])
  pq[ix[j, 1], ix[j, 2], ] <- quantile(ll[, j], c(.05, .95)) }
truef <- log(sim$truth$lambda_true_si)
fdf <- map_dfr(seq_len(S), function(s) tibble(species = factor(spn[s], spn, SP_LAB), truth = truef[, s],
        post = pm[, s], lo = pq[, s, 1], hi = pq[, s, 2], station = samp$station)) |>
  group_by(species, station) |> slice(1) |> ungroup()   # density is the same for every depth at a station
r2 <- fdf |> group_by(species) |> summarise(r2 = cor(truth, post)^2, cover90 = mean(truth >= lo & truth <= hi), n = n())
p2 <- ggplot(fdf, aes(truth, post)) + geom_abline(colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#2166ac", alpha = 0.7, size = 0.25) +
  geom_text(data = r2, aes(x = -Inf, y = Inf, label = sprintf("R2 = %.2f\n90%% band covers %.0f%%", r2, 100 * cover90)),
            hjust = -0.1, vjust = 1.2, size = 3.3, inherit.aes = FALSE) +
  facet_wrap(~species, scales = "free") +
  labs(x = "true log density", y = "posterior mean log density (90% interval)",
       title = "Latent field recovery at the sampling stations") + theme_expl()
save_fig(p2, "edna_sim_field", 12, 4.3)
saveRDS(r2, file.path(EXPL_OUT, "edna_sim_field_r2.rds"))

# ---- bottom-depth effect ------------------------------------------------------------------
zg <- seq(min(samp$Z_bathy), max(samp$Z_bathy), length.out = 80)
Bg <- bathy_spline_basis(zg, setup = sim$truth$bathy_spline)$B
bb <- pick("beta_bathy[")
sp_curves <- map_dfr(seq_len(S), function(s) {
  cols <- sprintf("beta_bathy[%d,%d]", s, seq_len(sd$K_bathy))
  eff <- Bg %*% t(as.matrix(bb[, cols]))
  tibble(species = factor(spn[s], spn, SP_LAB), Z = zg, med = apply(eff, 1, median),
         lo = apply(eff, 1, quantile, .05), hi = apply(eff, 1, quantile, .95),
         truth = as.numeric(Bg %*% BATHY_BETA_TRUE[spn[s], seq_len(sd$K_bathy)]))
})
# both curves are centred against the same basis; compare after removing the mean over the grid
sp_curves <- sp_curves |> group_by(species) |> mutate(across(c(med, lo, hi, truth), ~ . - mean(truth))) |> ungroup()
p3 <- ggplot(sp_curves, aes(Z)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#2166ac", alpha = 0.25) +
  geom_line(aes(y = med), colour = "#2166ac", linewidth = 1) + geom_line(aes(y = truth), colour = "#b2182b", linetype = 2, linewidth = 1) +
  facet_wrap(~species) + labs(x = "bottom depth (m)", y = "h(Z_bathy), centred",
  title = "Bottom-depth spline: posterior (blue, 90%) vs truth (red dashed)") + theme_expl()
save_fig(p3, "edna_sim_spline", 12, 3.8)

# ---- posterior predictive: metabarcoding positives by species x depth -----------------------
# The model's generated quantities (pp_mb_reads) omit the zero-inflation term of the likelihood, so
# the probability that a replicate has any reads is computed here, exactly as in the likelihood:
#   P(reads > 0) = (1 - exp(-aliquot copies)) * (1 - BetaBinomial(0 | total, a, b)).
le_all <- pick("log_lambda_edna[")
li_all <- mat_idx(le_all)
ord <- order(li_all[, 2], li_all[, 1])
le_arr <- array(as.matrix(le_all)[, ord], c(nrow(le_all), N, S + 1))             # draws x N x (S+1)
b0 <- as.numeric(draws[, "beta0_phi"]); g0 <- as.numeric(draws[, "gamma0_phi"]); g1 <- as.numeric(draws[, "gamma1_phi"])
log1mexp <- function(x) log(-expm1(x))
lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
mb_idx <- sd$mb_sample_idx; mb_tot <- sd$mb_total
dsel <- round(seq(1, nrow(le_all), length.out = 200))
ppos <- matrix(NA_real_, length(dsel), length(mb_idx) * S)
for (k in seq_along(dsel)) {
  d <- dsel[k]
  for (r in seq_along(mb_idx)) {
    ll <- le_arr[d, mb_idx[r], ] + lvf                                        # K log aliquot copies
    llK <- ll - log1mexp(-exp(ll))
    lnu <- llK - llK[S + 1]; lpi <- lnu - lse(lnu)
    lphi <- b0[d] + lse(ll) + pmax(g0[d] - g1[d] * llK, 0)
    a <- exp(lpi + lphi); bb <- exp(log1mexp(lpi) + lphi)
    p_bb0 <- exp(lbeta(a, mb_tot[r] + bb) - lbeta(a, bb))                     # BB(0 | total)
    ppos[k, (r - 1) * S + 1:S] <- (1 - exp(-exp(ll[1:S]))) * (1 - p_bb0[1:S])
  }
}
mb_d <- samp$Z_sample[mb_idx]
obs_reads <- sd$mb_reads
ppc_mb <- map_dfr(seq_len(S), function(s) {
  cols <- (seq_along(mb_idx) - 1) * S + s
  map_dfr(sort(unique(mb_d)), function(dp) {
    r <- which(mb_d == dp)
    pr <- rowMeans(ppos[, cols[r], drop = FALSE])
    tibble(species = factor(spn[s], spn, SP_LAB), depth = dp, obs = mean(obs_reads[r, s] > 0),
           lo = quantile(pr, .05), med = median(pr), hi = quantile(pr, .95))
  })
})
p4 <- ggplot(ppc_mb, aes(factor(depth))) + geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac") +
  geom_point(aes(y = obs), colour = "#b2182b", shape = 4, size = 3, stroke = 1.3) +
  facet_wrap(~species) + labs(x = "sampling depth (m)", y = "share of MB replicates with reads",
  title = "Posterior predictive check: metabarcoding positives",
  subtitle = "blue: model-implied probability of reads (median, 90% over posterior draws); red x: observed") + theme_expl()
save_fig(p4, "edna_sim_ppc_mb", 12, 3.8)
saveRDS(ppc_mb, file.path(EXPL_OUT, "edna_sim_ppc_mb.rds"))

# ---- posterior predictive: qPCR detection by depth, and Ct --------------------------------------
le <- pick("log_lambda_edna[")
li <- mat_idx(le)
hake_cols <- which(li[, 2] == 1); hake_rows <- li[hake_cols, 1]
le_h <- as.matrix(le[, hake_cols])[, order(hake_rows), drop = FALSE]    # draws x N
p_det <- 1 - exp(-sd$kappa * exp(le_h + lvf))                           # P(detect) per sample
qi <- sd$qpcr_sample_idx; qd <- samp$Z_sample[qi]
ppc_q <- map_dfr(sort(unique(qd)), function(d) {
  r <- which(qd == d)
  rate <- rowMeans(p_det[, qi[r], drop = FALSE])
  tibble(depth = d, obs = mean(sd$qpcr_detect[r]), lo = quantile(rate, .05), med = median(rate), hi = quantile(rate, .95), n = length(r))
})
p5a <- ggplot(ppc_q, aes(factor(depth))) + geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac") +
  geom_point(aes(y = obs), colour = "#b2182b", shape = 4, size = 3, stroke = 1.3) +
  labs(x = "sampling depth (m)", y = "qPCR detection rate (hake)", title = "qPCR detection by depth",
       subtitle = "blue: predicted (90%); red x: observed") + theme_expl()
cop <- exp(le_h + lvf)                                               # copies in the reaction
det_i <- which(sd$qpcr_detect == 1)
ct_pred <- apply(cop[, qi[det_i], drop = FALSE], 2, function(x) median(sd$alpha_ct - sd$beta_ct * log(pmax(x, 1))))
ctdf <- tibble(obs = sd$qpcr_ct[det_i], pred = ct_pred)
p5b <- ggplot(ctdf, aes(pred, obs)) + geom_point(alpha = 0.4, size = 0.9) + geom_abline(colour = "#b2182b") +
  labs(x = "predicted Ct (posterior median)", y = "observed Ct", title = "Ct: observed vs predicted",
       subtitle = sprintf("RMSE %.2f cycles", sqrt(mean((ctdf$obs - ctdf$pred)^2)))) + theme_expl()
save_fig(p5a + p5b, "edna_sim_ppc_qpcr", 11, 4)

# ---- prior predictive ---------------------------------------------------------------------------
# Draw parameters from the model's priors and push them through the same likelihood inputs:
# exact anisotropic GP at the bottles (the HSGP approximates it), spline, NB bottle noise.
set.seed(7)
coords_km <- as.matrix(samp[, c("X", "Y")])
B <- sd$B_bathy
nd <- 200
pp <- vector("list", nd)
for (d in seq_len(nd)) {
  cv <- abs(rnorm(1, 0, sd$prior_cv_bottle_sig)); th <- 1 / cv^2
  lg <- matrix(NA, N, S)
  for (s in 1:S) {
    gsig <- rgamma(1, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate)
    ell  <- c(1 / rgamma(1, sd$prior_gp_l_shape[1], rate = sd$prior_gp_l_scale[1]),
              1 / rgamma(1, sd$prior_gp_l_shape[2], rate = sd$prior_gp_l_scale[2])) * cs
    f <- as.numeric(MASS::mvrnorm(1, rep(0, N), aniso_cov_2d(coords_km, gsig, ell[1], ell[2])))
    mu <- rnorm(1, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)
    h <- as.numeric(B %*% rnorm(sd$K_bathy, 0, sd$prior_beta_bathy_sig))
    g <- rgamma(N, th, th)
    lg[, s] <- mu + f + h + sd$log_zsample_effect[, s] + sd$log_conv_factor[s] + sd$log_vol_filtered + log(g)
  }
  pp[[d]] <- lg + lvf              # log expected copies in the aliquot, per sample x species
}
pr_cop <- map_dfr(seq_len(nd), function(d) tibble(draw = d, species = rep(spn, each = N), l10 = as.numeric(pp[[d]]) / log(10)))
tru_cop <- tibble(species = rep(spn, each = N), l10 = as.numeric(log(sim$truth$lambda_true_si) +
  sd$log_zsample_effect + rep(sd$log_conv_factor[1:S], each = N) + sd$log_vol_filtered + lvf) / log(10))
post_cop <- tibble(species = rep(spn, each = N),
                   l10 = as.numeric(sapply(1:S, function(s) { j <- which(li[, 2] == s); r <- li[j, 1]
                     colMeans(as.matrix(le[, j]))[order(r)] })) / log(10) + lvf / log(10))
prior_sum <- pr_cop |> group_by(species, draw) |> summarise(q = list(quantile(l10, c(.1, .5, .9))), .groups = "drop") |>
  mutate(lo = map_dbl(q, 1), med = map_dbl(q, 2), hi = map_dbl(q, 3)) |>
  mutate(species = factor(species, spn, SP_LAB))
obs_q <- tru_cop |> group_by(species) |> summarise(lo = quantile(l10, .1), med = median(l10), hi = quantile(l10, .9)) |>
  mutate(species = factor(species, spn, SP_LAB))
p6a <- ggplot(prior_sum, aes(med)) + geom_histogram(bins = 30, fill = "grey75", colour = "white") +
  geom_vline(data = obs_q, aes(xintercept = med), colour = "#b2182b", linewidth = 1) +
  facet_wrap(~species, scales = "free_y") + labs(x = "median log10 expected copies per aliquot across bottles", y = "prior draws",
  title = "Prior predictive: typical eDNA signal per bottle", subtitle = "red: the simulated truth") + theme_expl()
# prior predictive of the hake qPCR detection rate (the single most informative summary)
rate_pp <- map_dbl(pp, function(lg) mean(1 - exp(-sd$kappa * exp(lg[qi, 1]))))
p6b <- ggplot(tibble(rate = rate_pp), aes(rate)) + geom_histogram(bins = 25, fill = "grey75", colour = "white") +
  geom_vline(xintercept = mean(sd$qpcr_detect), colour = "#b2182b", linewidth = 1) +
  geom_vline(xintercept = median(rowMeans(p_det[, qi, drop = FALSE])), colour = "#2166ac", linewidth = 1, linetype = 2) +
  labs(x = "hake qPCR detection rate", y = "prior draws", title = "Prior predictive: hake qPCR detection rate",
       subtitle = "red: observed; blue dashed: posterior prediction") + theme_expl()
save_fig(p6a / p6b, "edna_sim_prior_pc", 10, 7)
saveRDS(list(prior_cover_truth = mean(obs_q$med >= tapply(prior_sum$med, prior_sum$species, quantile, .01) &
                                       obs_q$med <= tapply(prior_sum$med, prior_sum$species, quantile, .99)),
             rate_range = quantile(rate_pp, c(.05, .5, .95)), obs_rate = mean(sd$qpcr_detect)),
        file.path(EXPL_OUT, "edna_sim_prior_pc.rds"))
cat("done\n")

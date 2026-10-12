# Figures and numbers for the eDNA model on the real MURI bottles (metabarcoding only).
# SRC=real (default) uses outputs/explainer/edna_real_*; SRC=sim runs the same code on the simulated fit
# (a test of this script). Needs the fit from 04_fit_edna_real.R.
source("docs/explainer/00_common.R")
source(file.path(HSGP4EDNA_DIR, "R/functions.R"))
SRC <- Sys.getenv("SRC", "real")
sd  <- readRDS(file.path(EXPL_OUT, sprintf("edna_%s_standata.rds", SRC)))
fit <- readRDS(file.path(EXPL_OUT, sprintf("edna_%s_fit.rds", SRC)))
samp <- if (SRC == "real") readRDS("data/processed/muri_edna.rds")$design$samples else
  readRDS(file.path(EXPL_OUT, "edna_sim.rds"))$design$samples
S <- 3; N <- sd$N; spn <- names(SP_LAB); lvf <- log(sd$vol_aliquot / 100)
draws <- as_draws_matrix(fit$draws())
pick <- function(prefix) draws[, startsWith(colnames(draws), prefix), drop = FALSE]
mat_idx <- function(v) { m <- str_match(colnames(v), "\\[(\\d+),(\\d+)\\]"); cbind(as.integer(m[, 2]), as.integer(m[, 3])) }
tag <- if (SRC == "real") "edna_real" else "edna_realtest"

# ---- diagnostics ------------------------------------------------------------------------------
dg <- fit$diagnostic_summary(quiet = TRUE)
hv <- grep("^gp_l\\[|^gp_sigma|^mu_sp|^theta", colnames(draws), value = TRUE)
sm <- summarise_draws(draws[, hv], rhat, ess_bulk)
th <- summarise_draws(draws[, grep("^theta", colnames(draws), value = TRUE)], median,
                      q05 = ~quantile(.x, .05), q95 = ~quantile(.x, .95))
diag <- list(div = sum(dg$num_divergent), td = sum(dg$num_max_treedepth), ebfmi = dg$ebfmi,
             max_rhat = max(sm$rhat, na.rm = TRUE), min_ess = min(sm$ess_bulk, na.rm = TRUE),
             runtime_min = max(fit$time()$chains$total) / 60, N = N, M = sd$M, m_hsgp = sd$m_hsgp,
             N_mb = sd$N_mb_long, theta = as.data.frame(th))

# ---- hyperparameters against the prior + floor check ----------------------------------------------
cs <- sd$coord_scale
qig <- function(p, shape, scale) 1 / qgamma(1 - p, shape, rate = scale)
hp <- map_dfr(seq_len(S), function(s) tibble(species = factor(spn[s], spn, SP_LAB),
  param = factor(c("lx (km)", "ly (km)", "gp_sigma", "mu"), c("lx (km)", "ly (km)", "gp_sigma", "mu")),
  med = c(median(draws[, sprintf("gp_l[%d,1]", s)]), median(draws[, sprintf("gp_l[%d,2]", s)]),
          median(draws[, sprintf("gp_sigma[%d]", s)]), median(draws[, sprintf("mu_sp[%d]", s)])),
  lo = c(quantile(draws[, sprintf("gp_l[%d,1]", s)], .025), quantile(draws[, sprintf("gp_l[%d,2]", s)], .025),
         quantile(draws[, sprintf("gp_sigma[%d]", s)], .025), quantile(draws[, sprintf("mu_sp[%d]", s)], .025)),
  hi = c(quantile(draws[, sprintf("gp_l[%d,1]", s)], .975), quantile(draws[, sprintf("gp_l[%d,2]", s)], .975),
         quantile(draws[, sprintf("gp_sigma[%d]", s)], .975), quantile(draws[, sprintf("mu_sp[%d]", s)], .975))))
prior_rng <- tibble(param = factor(c("lx (km)", "ly (km)", "gp_sigma", "mu"), levels(hp$param)),
  lo = c(qig(.025, sd$prior_gp_l_shape[1], sd$prior_gp_l_scale[1]) * cs[1], qig(.025, sd$prior_gp_l_shape[2], sd$prior_gp_l_scale[2]) * cs[2],
         qgamma(.025, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate), qnorm(.025, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)),
  hi = c(qig(.975, sd$prior_gp_l_shape[1], sd$prior_gp_l_scale[1]) * cs[1], qig(.975, sd$prior_gp_l_shape[2], sd$prior_gp_l_scale[2]) * cs[2],
         qgamma(.975, sd$prior_gp_sigma_shape, sd$prior_gp_sigma_rate), qnorm(.975, sd$prior_mu_sp_mu, sd$prior_mu_sp_sig)))
p1 <- ggplot(hp, aes(species)) +
  geom_rect(data = prior_rng, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey88") +
  geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac") + facet_wrap(~param, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "eDNA model on the MURI bottles: GP hyperparameters",
       subtitle = "blue: posterior 95%; grey: central 95% of the design-based prior") + theme_expl() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_fig(p1, paste0(tag, "_hyper"), 12, 4)
gpp <- attr(sd, "gp_prior")
floor_tab <- NULL
if (!is.null(gpp)) floor_tab <- map_dfr(seq_len(S), function(s)
  cbind(species = spn[s], floor_check(cbind(draws[, sprintf("gp_l[%d,1]", s)], draws[, sprintf("gp_l[%d,2]", s)]), gpp)))

# ---- posterior mean density map (log10 animals per km2, nominal units) -----------------------------
ll <- pick("log_lambda["); ll <- ll[, !startsWith(colnames(ll), "log_lambda_edna") & !startsWith(colnames(ll), "log_lambda_pred"), drop = FALSE]
ix <- mat_idx(ll); pm <- matrix(NA, N, S)
for (j in seq_len(ncol(ll))) pm[ix[j, 1], ix[j, 2]] <- mean(ll[, j])
mp <- map_dfr(seq_len(S), function(s) tibble(species = factor(spn[s], spn, SP_LAB), X = samp$X, Y = samp$Y,
                                             station = samp$station, l10 = pm[, s] / log(10))) |>
  group_by(species, station) |> summarise(X = first(X), Y = first(Y), l10 = mean(l10), .groups = "drop")
p2 <- ggplot(mp, aes(X, Y, colour = l10)) + geom_point(size = 2) + coord_equal() + facet_wrap(~species) +
  scale_colour_viridis_c(name = "log10\ndensity\n(nominal)") +
  labs(x = "X (km)", y = "Y (km)", title = "Posterior mean latent density at the MURI stations",
       subtitle = "nominal units: placeholder conversion factors and volumes (see text)") + theme_expl()
save_fig(p2, paste0(tag, "_map"), 12, 6)

# ---- posterior predictive: share of replicates with reads, by species x depth ----------------------------
le_all <- pick("log_lambda_edna["); li <- mat_idx(le_all)
ord <- order(li[, 2], li[, 1])
le_arr <- array(as.matrix(le_all)[, ord], c(nrow(le_all), N, S + 1))
b0 <- as.numeric(draws[, "beta0_phi"]); g0 <- as.numeric(draws[, "gamma0_phi"]); g1 <- as.numeric(draws[, "gamma1_phi"])
log1mexp <- function(x) log(-expm1(x)); lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
mb_idx <- sd$mb_sample_idx; mb_tot <- sd$mb_total
dsel <- round(seq(1, nrow(le_all), length.out = 100))
ppos <- matrix(NA_real_, length(dsel), length(mb_idx) * S)
for (k in seq_along(dsel)) { d <- dsel[k]
  for (r in seq_along(mb_idx)) {
    ll_r <- le_arr[d, mb_idx[r], ] + lvf; llK <- ll_r - log1mexp(-exp(ll_r))
    lnu <- llK - llK[S + 1]; lpi <- lnu - lse(lnu); lphi <- b0[d] + lse(ll_r) + pmax(g0[d] - g1[d] * llK, 0)
    a <- exp(lpi + lphi); bb <- exp(log1mexp(lpi) + lphi)
    ppos[k, (r - 1) * S + 1:S] <- (1 - exp(-exp(ll_r[1:S]))) * (1 - exp(lbeta(a, mb_tot[r] + bb) - lbeta(a, bb)))[1:S]
  } }
dep <- samp$Z_sample[mb_idx]
dep_bin <- cut(dep, c(-1, 1, 75, 200, 400, 1000), labels = c("0", "50", "150", "300", "500"))
ppc <- map_dfr(seq_len(S), function(s) { cols <- (seq_along(mb_idx) - 1) * S + s
  map_dfr(levels(dep_bin), function(b) { r <- which(dep_bin == b); if (!length(r)) return(NULL)
    pr <- rowMeans(ppos[, cols[r], drop = FALSE])
    tibble(species = factor(spn[s], spn, SP_LAB), depth = factor(b, levels(dep_bin)), n = length(r), obs = mean(sd$mb_reads[r, s] > 0),
           lo = quantile(pr, .05), med = median(pr), hi = quantile(pr, .95)) }) })
p3 <- ggplot(ppc, aes(depth)) + geom_pointrange(aes(y = med, ymin = lo, ymax = hi), colour = "#2166ac") +
  geom_point(aes(y = obs), colour = "#b2182b", shape = 4, size = 3, stroke = 1.3) + facet_wrap(~species) +
  labs(x = "sampling depth (m, binned)", y = "share of MB replicates with reads",
       title = "Posterior predictive check: metabarcoding positives, MURI bottles",
       subtitle = "blue: model-implied probability (median, 90%); red x: observed") + theme_expl()
save_fig(p3, paste0(tag, "_ppc_mb"), 12, 3.8)
saveRDS(list(diag = diag, floor = floor_tab, hyper = hp, ppc = ppc), file.path(EXPL_OUT, paste0(tag, "_results.rds")))
cat("done\n")

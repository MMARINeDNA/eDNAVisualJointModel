# =============================================================================
# fit_qpcr_standards.r
#
# Hake qPCR standard curves, one per plate, estimated from the gBlocks
# standards in the qPCR lab's qPCR_2019A.rds (stan/qpcr_standards.stan, after
# the qPCR vignette's multi-plate model and priors).
#
#   Rscript R/fit_qpcr_standards.r
# Env: RAW_DIR (folder with qPCR_2019A.rds; default data), INHIBIT (standards
#      with inhibition_rate at or above this are dropped, as in the vignette;
#      default 0.75), CHAINS, WARMUP, SAMPLE (default 4 x 1000/1000)
# Out: data/processed/qpcr_standard_curves.rds  per-plate alpha_ct, beta_ct
#        (posterior mean, sd, 95% interval), the shared gamma0, gamma1,
#        sigma0, kappa, 1000 joint posterior draws, diagnostics and checks;
#        read by R/01_load_muri.r
#      outputs/qpcr/standard_curves.{png,csv}, outputs/qpcr/checks.txt
# =============================================================================

suppressMessages({library(cmdstanr); library(posterior); library(dplyr); library(ggplot2)})

env_int <- function(k, d) as.integer(Sys.getenv(k, d))
RAW_DIR <- Sys.getenv("RAW_DIR", "data")
INHIBIT <- as.numeric(Sys.getenv("INHIBIT", "0.75"))
OUT <- "outputs/qpcr"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
VOL_REACTION <- 2      # uL of template per reaction (qPCR vignette: 2 uL for standards and samples)

A  <- readRDS(file.path(RAW_DIR, "qPCR_2019A.rds"))
st <- as_tibble(A$standards) |> filter(species == "hake")
n_all <- nrow(st)
st <- st |> filter(is.na(inhibition_rate) | inhibition_rate < INHIBIT)
plates <- sort(unique(st$qPCR))
st <- st |> mutate(plate = match(qPCR, plates), detect = as.integer(Ct_bin == 1),
                   ct = ifelse(detect == 1, Ct, 0))
stopifnot(all(st$ct[st$detect == 1] > 0), all(st$copies_ul > 0))
cat(sprintf("hake standards: %d wells, %d dropped as inhibited (>= %.2f), %d used on %d plates; %d non-detects\n",
            n_all, n_all - nrow(st), INHIBIT, nrow(st), length(plates), sum(st$detect == 0)))

sdat <- list(N = nrow(st), N_plate = length(plates), plate = st$plate, detect = st$detect,
           ct = st$ct, K = st$copies_ul, vol_reaction = VOL_REACTION,
           beta0_prior = c(40, 5), beta1_prior = c(-1 / log(2), 0.20), beta1_max = -1 / log(2),
           gamma0_prior = c(0, 2), gamma1_prior = c(-0.2, 0.5), logit_kappa_prior = c(1, 3))

mod <- cmdstan_model("stan/qpcr_standards.stan")
set.seed(1)
fit <- mod$sample(data = sdat, chains = env_int("CHAINS", "4"), parallel_chains = env_int("CHAINS", "4"),
                  iter_warmup = env_int("WARMUP", "1000"), iter_sampling = env_int("SAMPLE", "1000"),
                  seed = 1, refresh = 0, show_messages = FALSE,
                  init = function() list(beta0 = rnorm(sdat$N_plate, 38, 1),
                                         beta1 = rep(-1.6, sdat$N_plate) + rnorm(sdat$N_plate, 0, 0.02),
                                         sigma0 = 0.3, gamma0 = 0, gamma1 = -0.2, logit_kappa = 2))

q95 <- function(x) quantile(x, c(0.025, 0.975))
per_plate <- function(v) fit$summary(v, mean, sd, q95, "rhat", "ess_bulk") |>
  mutate(plate = plates[as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable))])
pa <- per_plate("alpha_ct"); pb <- per_plate("beta_ct"); pe <- per_plate("efficiency")
curves <- tibble(plate = plates,
                 alpha_ct = pa$mean, alpha_ct_sd = pa$sd, alpha_ct_q025 = pa$`2.5%`, alpha_ct_q975 = pa$`97.5%`,
                 beta_ct = pb$mean, beta_ct_sd = pb$sd, beta_ct_q025 = pb$`2.5%`, beta_ct_q975 = pb$`97.5%`,
                 efficiency = pe$mean, n_standards = as.integer(table(factor(st$qPCR, plates))))
shared <- fit$summary(c("gamma0", "gamma1", "sigma0", "kappa"), mean, sd, q95, "rhat", "ess_bulk")
diag <- fit$diagnostic_summary(quiet = TRUE)
max_rhat <- max(fit$summary(c("beta0", "beta1", "gamma0", "gamma1", "sigma0", "logit_kappa"), "rhat")$rhat)

# ---- Checks -------------------------------------------------------------------
# 1. predictive: share of detected standards whose Ct lies in its 95% predictive interval
rep <- fit$draws("ct_rep", format = "draws_matrix")
pi <- apply(rep, 2, quantile, c(0.025, 0.975))
pos <- st$detect == 1
cover <- mean(st$ct[pos] >= pi[1, pos] & st$ct[pos] <= pi[2, pos])
# 2. detection by concentration: observed against the fitted 1 - exp(-kappa * 2K)
kap <- shared$mean[shared$variable == "kappa"]
det_tab <- st |> group_by(copies_ul) |>
  summarise(n = n(), observed = mean(detect), fitted = 1 - exp(-kap * VOL_REACTION * first(copies_ul)), .groups = "drop")
# 3. negatives: hake detections in field negatives and extraction / PCR controls
ctl <- as_tibble(A$controls) |> filter(species == "hake")
neg <- ctl |> group_by(type = coalesce(field.negative.type, type)) |>
  summarise(wells = n(), detected = sum(Ct > 0, na.rm = TRUE), max_copies_ul = suppressWarnings(max(copies_ul, na.rm = TRUE)),
            .groups = "drop")

res <- list(curves = curves, shared = shared,
            draws = as_draws_df(fit$draws(c("alpha_ct", "beta_ct", "gamma0", "gamma1", "sigma0", "kappa"))) |>
              (\(d) d[round(seq(1, nrow(d), length.out = 1000)), ])(),
            meta = list(source = "qPCR_2019A.rds (hake standards)", vol_reaction = VOL_REACTION,
                        inhibit_cutoff = INHIBIT, n_standards = nrow(st), n_dropped_inhibited = n_all - nrow(st),
                        model = "stan/qpcr_standards.stan", fitted = Sys.time()),
            diagnostics = list(divergences = sum(diag$num_divergent), treedepth = sum(diag$num_max_treedepth),
                               max_rhat = max_rhat, ebfmi = diag$ebfmi),
            checks = list(ct_coverage_95 = cover, detection = det_tab, negatives = neg))
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "data/processed/qpcr_standard_curves.rds")
write.csv(curves, file.path(OUT, "standard_curves.csv"), row.names = FALSE)

sink(file.path(OUT, "checks.txt"), split = TRUE)
cat(sprintf("divergences %d | treedepth %d | max Rhat %.3f\n", res$diagnostics$divergences,
            res$diagnostics$treedepth, max_rhat))
cat("shared parameters:\n"); print(as.data.frame(shared), digits = 3)
cat(sprintf("per-plate intercept alpha_ct %.2f-%.2f, slope beta_ct %.3f-%.3f, efficiency %.2f-%.2f\n",
            min(curves$alpha_ct), max(curves$alpha_ct), min(curves$beta_ct), max(curves$beta_ct),
            min(curves$efficiency), max(curves$efficiency)))
cat(sprintf("detected standards inside their 95%% predictive interval: %.3f\n", cover))
cat("detection by concentration:\n"); print(as.data.frame(det_tab), digits = 3)
cat("hake in negatives and controls:\n"); print(as.data.frame(neg))
sink()

# ---- Figure -------------------------------------------------------------------
grid <- expand.grid(plate = plates, K = 10^seq(0, 5, length.out = 30)) |>
  left_join(curves, by = "plate") |> mutate(ct = alpha_ct - beta_ct * log(VOL_REACTION * K))
g <- ggplot(st |> filter(detect == 1), aes(copies_ul, ct)) +
  geom_line(data = grid, aes(K, ct, group = plate), colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.7, alpha = 0.6, colour = "#0072B2") +
  geom_point(data = st |> filter(detect == 0), aes(copies_ul, 41), shape = 4, colour = "#D55E00") +
  scale_x_log10() + labs(x = "standard concentration (copies per uL)", y = "Ct",
                         title = sprintf("Hake qPCR standards, %d plates: per-plate fitted curves (grey)", length(plates)),
                         subtitle = "Crosses at Ct 41: non-detects") + theme_bw()
ggsave(file.path(OUT, "standard_curves.png"), g, width = 7, height = 5, dpi = 130)
cat("Saved data/processed/qpcr_standard_curves.rds and", OUT, "\n")

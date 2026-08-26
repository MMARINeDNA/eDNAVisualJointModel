# =============================================================================
# 07_plot_basis_sweep_vS1.R
#
# Summarise and plot the HSGP basis-count sweep produced by
# 06_basis_sweep_vS1.r. Aggregates over simulation seeds and produces:
#
#   outputs/whale_edna_output_vS1/basis_sweep/basis_sweep_summary.pdf
#   outputs/whale_edna_output_vS1/basis_sweep/basis_sweep_summary.csv
#
# Works on partial results (safe to run mid-sweep).
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(patchwork)
})

OUTPUT_DIR <- "outputs/whale_edna_output_vS1/basis_sweep"
res <- readRDS(file.path(OUTPUT_DIR, "sweep_results.rds"))
sp_common <- c("Pacific hake", "Humpback whale", "P. white-sided dolphin")

# order configs by M (then label), so H/I (same M) stay distinct on a
# categorical x-axis instead of overplotting
lab_levels <- res$diag %>% distinct(label, M) %>% arrange(M, label) %>% pull(label)
labf <- function(x) factor(x, levels = lab_levels)
xlab_M <- res$diag %>% distinct(label, M) %>%
  mutate(tick = sprintf("%s\n(%d)", label, M)) %>% arrange(match(label, lab_levels))

# -----------------------------------------------------------------------------
# Aggregate across seeds
# -----------------------------------------------------------------------------
field_ag <- res$field %>% group_by(label, M, species) %>%
  summarise(R2 = mean(R2), rmse = mean(rmse), bias = mean(bias),
            n = n(), .groups = "drop") %>%
  mutate(species = sp_common[species], label = labf(label))

rec_ag <- res$recovery %>% group_by(label, M, species, param) %>%
  summarise(mean = mean(mean), true = first(true),
            coverage = mean(covered), n = n(), .groups = "drop") %>%
  mutate(species = sp_common[species], label = labf(label))

diag_ag <- res$diag %>% group_by(label, M) %>%
  summarise(runtime_s = mean(runtime_s), max_rhat = max(max_rhat),
            n_divergent = sum(n_divergent), min_ess = min(min_ess),
            .groups = "drop") %>% mutate(label = labf(label))

xsc <- scale_x_discrete(limits = lab_levels, labels = xlab_M$tick)
thm <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
# 1. Latent-field recovery R^2 vs basis
p_field <- ggplot(field_ag, aes(labf(label), R2, colour = species, group = species)) +
  geom_line() + geom_point(size = 2) + xsc +
  scale_colour_viridis_d(end = 0.85) +
  labs(title = "Latent-field recovery vs basis count",
       subtitle = "posterior-mean log_lambda vs truth (mean over seeds)",
       x = "config (M = total basis)", y = expression(R^2)) + thm

# 2. Length-scale recovery (lx, ly) vs basis, with truth lines
ls_df <- rec_ag %>% filter(param %in% c("lx","ly")) %>%
  mutate(dim = recode(param, lx = "lx (cross-shore, true 50)",
                             ly = "ly (along-shore, true 300)"))
p_ls <- ggplot(ls_df, aes(labf(label), mean, colour = species, group = species)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed", colour = "grey40") +
  facet_wrap(~dim, scales = "free_y") + xsc +
  scale_colour_viridis_d(end = 0.85) +
  labs(title = "Length-scale recovery vs basis count",
       subtitle = "dashed = truth; overshoot at low M = oversmoothing, collapse at high M = overfit",
       x = "config (M)", y = "posterior mean length scale (km)") + thm

# 3. gp_sigma recovery
p_sig <- ggplot(rec_ag %>% filter(param == "gp_sigma"),
                aes(labf(label), mean, colour = species, group = species)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed", colour = "grey40") +
  xsc + scale_colour_viridis_d(end = 0.85) +
  labs(title = "Marginal SD (gp_sigma) recovery",
       x = "config (M)", y = "posterior mean gp_sigma") + thm

# 4. 95% CI coverage of truth, by parameter
p_cov <- ggplot(rec_ag, aes(labf(label), coverage, colour = param, group = param)) +
  geom_line() + geom_point(size = 2) + xsc +
  scale_colour_viridis_d(end = 0.9) + ylim(0, 1) +
  labs(title = "95% CI coverage of truth", x = "config (M)",
       y = "fraction of species x seeds covering truth") + thm

# 5. Cost + sampling health
p_time <- ggplot(diag_ag, aes(labf(label), runtime_s, group = 1)) +
  geom_line() + geom_point(size = 2) + xsc +
  labs(title = "Runtime", x = "config (M)", y = "seconds / fit") + thm
p_rhat <- ggplot(diag_ag, aes(labf(label), max_rhat, group = 1)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(yintercept = 1.01, linetype = "dashed", colour = "red") + xsc +
  labs(title = "Max Rhat (dashed = 1.01)", x = "config (M)", y = "max Rhat") + thm
p_div <- ggplot(diag_ag, aes(labf(label), n_divergent, group = 1)) +
  geom_line() + geom_point(size = 2) + xsc +
  labs(title = "Divergences (summed over seeds)", x = "config (M)", y = "# divergent") + thm

# -----------------------------------------------------------------------------
# Recommendation summary table (printed + CSV)
# -----------------------------------------------------------------------------
summary_tbl <- field_ag %>% filter(species == "Pacific hake") %>%
  select(label, M, hake_field_R2 = R2) %>%
  left_join(rec_ag %>% filter(species == "Pacific hake", param == "lx") %>%
              select(label, hake_lx_mean = mean, hake_lx_cov = coverage), by = "label") %>%
  left_join(diag_ag %>% select(label, runtime_s, max_rhat, n_divergent), by = "label") %>%
  arrange(M)
write_csv(summary_tbl, file.path(OUTPUT_DIR, "basis_sweep_summary.csv"))
cat("\n=== Basis-sweep summary (hake; true lx = 50) ===\n")
print(as.data.frame(summary_tbl), digits = 3, row.names = FALSE)

# -----------------------------------------------------------------------------
# Write PDF
# -----------------------------------------------------------------------------
pdf(file.path(OUTPUT_DIR, "basis_sweep_summary.pdf"), width = 11, height = 8.5, onefile = TRUE)
print(p_field)
print(p_ls)
print(p_sig / p_cov)
print((p_time | p_rhat | p_div) +
        plot_annotation(title = "Cost and sampling health vs basis count"))
invisible(dev.off())
cat(sprintf("\nWrote %s and .csv\n", file.path(OUTPUT_DIR, "basis_sweep_summary.pdf")))

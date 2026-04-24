# =============================================================================
# 00_pipeline_v4.R
#
# End-to-end pipeline for the V3 eDNA model. Run from the project root:
#   Rscript 00_pipeline_v4.r
#
# Executes, in order:
#   1. scripts/01_simulate_whale_edna_v4.r     — simulate eDNA data
#   2. scripts/02_plot_simulated_data_v4.r     — plot the simulated truth
#   3. scripts/03_format_stan_data_v4.r        — assemble the Stan data list
#   4. scripts/04_run_whale_edna_model_v4.r    — compile + fit Stan model
#   5. scripts/05_check_whale_edna_model_v4.r  — diagnostics + plots
#
# Outputs land in outputs/whale_edna_output_v4/ and outputs/.
# =============================================================================

cat("\n--- 00_pipeline_v4: simulate ---\n")
source("scripts/01_simulate_whale_edna_v4.r")

cat("\n--- 00_pipeline_v4: plot simulated data ---\n")
source("scripts/02_plot_simulated_data_v4.r")

cat("\n--- 00_pipeline_v4: format Stan data ---\n")
source("scripts/03_format_stan_data_v4.r")

cat("\n--- 00_pipeline_v4: run model ---\n")
source("scripts/04_run_whale_edna_model_v4.r")

cat("\n--- 00_pipeline_v4: check model ---\n")
source("scripts/05_check_whale_edna_model_v4.r")

cat("\n=== 00_pipeline_v4: done ===\n")

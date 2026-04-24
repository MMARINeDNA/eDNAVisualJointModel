# =============================================================================
# 00_pipeline_v1.R
#
# End-to-end pipeline for the V1 eDNA model. Run from the project root:
#   Rscript 00_pipeline_v1.r
#
# Executes, in order:
#   1. scripts/01_simulate_whale_edna_v1.r     — simulate eDNA data
#   2. scripts/03_format_stan_data_v1.r        — assemble the Stan data list
#   3. scripts/04_run_whale_edna_model_v1.r    — compile + fit Stan model
#   4. scripts/05_check_whale_edna_model_v1.r  — diagnostics + plots
#
# Outputs land in outputs/whale_edna_output_v1/.
# =============================================================================

cat("\n--- 00_pipeline_v1: simulate ---\n")
source("scripts/01_simulate_whale_edna_v1.r")

cat("\n--- 00_pipeline_v1: format Stan data ---\n")
source("scripts/03_format_stan_data_v1.r")

cat("\n--- 00_pipeline_v1: run model ---\n")
source("scripts/04_run_whale_edna_model_v1.r")

cat("\n--- 00_pipeline_v1: check model ---\n")
source("scripts/05_check_whale_edna_model_v1.r")

cat("\n=== 00_pipeline_v1: done ===\n")

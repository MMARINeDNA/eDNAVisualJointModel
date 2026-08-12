# =============================================================================
# helper_functions.R
#
# Shared helpers for the vS1 simulation/estimation pipeline.
# =============================================================================

# -----------------------------------------------------------------------------
# pred_summary()
#
# Summarise draws of one or more matrix-valued Stan quantities (e.g.
# log_lambda[N, S], log_lambda_pred[N_pred, S+1]) into per-element
# posterior mean / SD / 2.5% / 97.5%, in both wide and long form.
#
# Takes a cmdstanr CmdStanMCMC fit (cmdstanr is the canonical runner).
# For each named parameter it returns a list with:
#   Mean, SD, q025, q975 : data frames, one row per first index (idx),
#                          columns sp_1..sp_C for the second index, plus idx
#   Long                 : long form (idx, species, Mean, SD, q025, q975, param)
#
# Each parameter must be indexed as par[i, j] (a 2-D array / matrix), which
# is the case for every quantity summarised in the vS1 pipeline.
# -----------------------------------------------------------------------------
pred_summary <- function(fit, pars) {
  pred_out <- list()
  for (par in pars) {
    # as.numeric() strips quantile()'s inner names so the "q025"/"q975"
    # column names stick.
    smry <- posterior::summarise_draws(
      fit$draws(variables = par),
      mean = mean, sd = sd,
      q025 = ~as.numeric(stats::quantile(.x, 0.025)),
      q975 = ~as.numeric(stats::quantile(.x, 0.975))
    )

    # Parse "par[i,j]" -> integer i (row / idx) and j (column / species).
    m  <- stringr::str_match(smry$variable, "\\[(\\d+),(\\d+)\\]")
    ii <- as.integer(m[, 2])
    jj <- as.integer(m[, 3])
    R  <- max(ii)
    C  <- max(jj)

    to_wide <- function(vals) {
      mat <- matrix(NA_real_, R, C)
      mat[cbind(ii, jj)] <- vals
      d <- as.data.frame(mat)
      colnames(d) <- paste0("sp_", seq_len(C))
      d$idx <- seq_len(R)
      d
    }

    Long <- tibble::tibble(
      idx     = ii,
      species = paste0("sp_", jj),
      Mean    = smry$mean,
      SD      = smry$sd,
      q025    = smry$q025,
      q975    = smry$q975,
      param   = par
    ) |> dplyr::arrange(idx, species)

    pred_out[[par]] <- list(
      Mean = to_wide(smry$mean),
      SD   = to_wide(smry$sd),
      q025 = to_wide(smry$q025),
      q975 = to_wide(smry$q975),
      Long = Long
    )
  }
  pred_out
}

# =============================================================================
# plot_joint.r
#
# Figures from the joint model's saved fit summaries (fit_joint() results with
# a prediction grid and use_ppc = 1; the fits themselves are not needed):
#
#   plot_field_maps()    species (rows) x truth / eDNA only / visual only /
#                        joint (columns): posterior mean log density on the
#                        grid; with what = "sd", the posterior SD instead
#   grid_recovery()      R^2, RMSE, 90% coverage and interval width on the grid
#                        against truth_on_grid()
#   plot_ppc_visual()    per line-transect species: detections along the coast,
#                        detected distances, detected group sizes and segment
#                        count statistics, observed against predictive
#   plot_ppc_edna()      metabarcoding share with reads and qPCR detection by
#                        sampling depth; observed against predicted Ct
#   plot_prior_pred()    prior predictive totals against the observed
#
#   source("R/plot_joint.r")
#   res <- lapply(c(edna = ..., visual = ..., joint = ...), readRDS)
#   p <- plot_field_maps(res, grid, sim, truth = truth_on_grid(sim, grid))
# =============================================================================

suppressMessages({library(ggplot2); library(patchwork)})

PAL <- list(model = "#2a78d6", obs = "#eb6834", ink = "#3d3d3a", muted = "#8a8984",
            grid = "#e8e7e2", ramp = c("#f3f8fe", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                                       "#256abf", "#184f95", "#0d366b"))
MODE_LABEL <- c(truth = "Truth", edna = "eDNA only", visual = "Visual only", joint = "Joint",
                joint_conv = "Joint, conversion free")
SP_LABEL <- c(hake = "Hake", humpback = "Humpback", pwsd = "PWSD")

theme_joint <- function(base = 10) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = PAL$grid, linewidth = 0.3),
        axis.text = element_text(colour = PAL$muted), axis.title = element_text(colour = PAL$ink),
        plot.title = element_text(colour = PAL$ink, face = "bold", size = base + 1),
        plot.subtitle = element_text(colour = PAL$muted, size = base - 1),
        strip.text = element_text(colour = PAL$ink, face = "bold"),
        legend.title = element_text(colour = PAL$ink, size = base - 1),
        legend.text = element_text(colour = PAL$muted, size = base - 2))

# ---- Maps ----------------------------------------------------------------------
# res: named list of fit summaries (names among edna, visual, joint,
# joint_conv), each with $grid; grid: the prediction points (X, Y); obj: the
# simulate_joint() or real_joint_object() object (for the sample locations);
# truth: truth_on_grid() output or NULL (real data). what: "mean" or "sd".
plot_field_maps <- function(res, grid, obj, truth = NULL, what = c("mean", "sd"),
                            modes = intersect(c("edna", "visual", "joint"), names(res)),
                            title = NULL, subtitle = NULL, show_samples = TRUE) {
  what <- match.arg(what)
  sp <- names(res[[modes[1]]]$grid)
  lt <- obj$meta$lt_species
  cols <- c(if (!is.null(truth) && what == "mean") "truth", modes)
  st <- unique(obj$edna$design$samples[, c("X", "Y")])
  sg <- obj$visual$design$segments[, c("X", "Y")]
  spacing <- attr(grid, "spacing") %||% min(diff(sort(unique(grid$X))))
  val <- function(col, s) {
    if (col == "truth") return(truth$log_lambda[, s])
    if (col == "visual" && !(s %in% lt)) return(NULL)            # not surveyed
    res[[col]]$grid[[s]][[what]]
  }
  rows <- lapply(sp, function(s) {
    v <- unlist(lapply(cols, val, s = s))
    lim <- if (what == "mean") quantile(v, c(0.01, 0.99), na.rm = TRUE) else c(0, quantile(v, 0.99, na.rm = TRUE))
    panels <- lapply(seq_along(cols), function(j) {
      col <- cols[j]; z <- val(col, s)
      p <- ggplot() + coord_equal(expand = FALSE) + theme_joint(9) +
        theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid.major = element_blank(),
              panel.background = element_rect(fill = "#fafaf8", colour = NA),
              plot.title = element_text(size = 9.5, face = "plain", colour = PAL$ink, hjust = 0.5))
      if (is.null(z)) {
        p <- p + annotate("text", x = mean(range(grid$X)), y = mean(range(grid$Y)), label = "not\nsurveyed",
                          colour = PAL$muted, size = 3.2) +
          scale_x_continuous(limits = range(grid$X)) + scale_y_continuous(limits = range(grid$Y))
      } else {
        d <- data.frame(X = grid$X, Y = grid$Y, z = pmin(pmax(z, lim[1]), lim[2]))
        p <- p + geom_tile(data = d, aes(X, Y, fill = z), width = spacing, height = spacing) +
          scale_fill_gradientn(colours = PAL$ramp, limits = lim,
                               name = if (what == "mean") expression(atop("log density", (animals / km^2)))
                                      else expression(atop("posterior SD", (log~scale))))
        if (show_samples) {
          if (col %in% c("truth", "edna", "joint", "joint_conv"))
            p <- p + geom_point(data = st, aes(X, Y), shape = 21, size = 0.55, stroke = 0.25,
                                colour = "white", fill = PAL$ink)
          if (col %in% c("truth", "visual", "joint", "joint_conv") && s %in% lt)
            p <- p + geom_point(data = sg, aes(X, Y), size = 0.08, colour = PAL$obs, alpha = 0.6)
        }
      }
      p + labs(title = if (s == sp[1]) MODE_LABEL[[col]] else NULL,
               y = if (j == 1) SP_LABEL[[s]] else NULL) +
        theme(axis.title.y = if (j == 1) element_text(angle = 90, face = "bold", colour = PAL$ink, size = 10)
                             else element_blank())
    })
    wrap_plots(panels, nrow = 1) + plot_layout(guides = "collect")
  })
  wrap_plots(rows, ncol = 1) +
    plot_annotation(title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold", colour = PAL$ink, size = 12),
                                  plot.subtitle = element_text(colour = PAL$muted, size = 9)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Recovery on the grid against truth_on_grid() (90% intervals)
# (the visual-only fit's unsurveyed species, outside lt, are left out)
grid_recovery <- function(res, truth, lt = c("humpback", "pwsd")) {
  do.call(rbind, lapply(names(res), function(m) {
    g <- res[[m]]$grid; if (is.null(g)) return(NULL)
    do.call(rbind, lapply(names(g), function(s) {
      if (m == "visual" && !(s %in% lt)) return(NULL)
      tr <- truth$log_lambda[, s]; e <- g[[s]]
      data.frame(mode = m, species = s, R2_grid = cor(e$mean, tr)^2, rmse_grid = sqrt(mean((e$mean - tr)^2)),
                 cover90_grid = mean(tr >= e$q05 & tr <= e$q95), width90_grid = mean(e$q95 - e$q05))
    }))
  }))
}

# ---- Predictive checks: visual ---------------------------------------------------
pr_point <- function(p, colour = PAL$model) p + geom_pointrange(aes(ymin = lo, ymax = hi), colour = colour, size = 0.25, linewidth = 0.6)

plot_ppc_visual <- function(r, title = NULL) {
  pv <- r$ppc$visual
  rows <- lapply(names(pv), function(s) {
    v <- pv[[s]]
    b <- v$bands
    p1 <- ggplot(b, aes(Y_km)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = PAL$model, alpha = 0.2) +
      geom_line(aes(y = med), colour = PAL$model, linewidth = 0.7) +
      geom_point(aes(y = obs), colour = PAL$obs, shape = 4, size = 1.8, stroke = 0.9) +
      coord_flip() + labs(x = "Y, along the coast (km)", y = "groups detected per 100 km band",
                          title = sprintf("%s: detections along the coast", SP_LABEL[[s]])) + theme_joint(9)
    plots <- list(p1)
    if (!is.null(v$dist)) {
      d <- v$dist
      p2 <- ggplot() +
        geom_histogram(aes(x = d$obs, y = after_stat(density)), breaks = seq(0, d$w, length.out = 21),
                       fill = "#f4c6b3", colour = "white", linewidth = 0.3) +
        geom_ribbon(data = d$pred, aes(x, ymin = lo, ymax = hi), fill = PAL$model, alpha = 0.25) +
        geom_line(data = d$pred, aes(x, med), colour = PAL$model, linewidth = 0.8) +
        labs(x = "perpendicular distance (km)", y = "density", title = "detected distances") + theme_joint(9)
      z <- v$size
      p3 <- pr_point(ggplot(z, aes(class, med))) +
        geom_point(aes(y = obs), colour = PAL$obs, shape = 4, size = 2, stroke = 0.9) +
        labs(x = "group size", y = "share of detected groups",
             title = sprintf("detected group sizes (n = %d)", z$n_obs[1])) +
        theme_joint(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plots <- c(plots, list(p2, p3))
    }
    cs <- v$counts
    lab <- c(zero = "share of segments with none", disp = "variance / mean of segment counts")
    cl <- do.call(rbind, lapply(names(lab), function(k) data.frame(stat = lab[[k]], rep = cs$rep[[k]], obs = cs$obs[[k]])))
    p4 <- ggplot(cl, aes(rep)) + geom_histogram(bins = 25, fill = PAL$model, alpha = 0.5, colour = "white", linewidth = 0.2) +
      geom_vline(aes(xintercept = obs), colour = PAL$obs, linewidth = 0.9) +
      facet_wrap(~stat, scales = "free", ncol = 1) +
      labs(x = NULL, y = "replicates", title = "segment counts") + theme_joint(9)
    wrap_plots(c(plots, list(p4)), nrow = 1, widths = c(1, 1.1, 1.1, 1)[seq_len(length(plots) + 1)])
  })
  wrap_plots(rows, ncol = 1) +
    plot_annotation(title = title, subtitle = "Blue: predictive (median, 90% over draws); orange: observed",
                    theme = theme(plot.title = element_text(face = "bold", colour = PAL$ink, size = 12),
                                  plot.subtitle = element_text(colour = PAL$muted, size = 9)))
}

# ---- Predictive checks: eDNA ------------------------------------------------------
plot_ppc_edna <- function(r, title = NULL) {
  e <- r$ppc$edna
  mb <- e$mb; mb$species <- factor(SP_LABEL[mb$species], SP_LABEL)
  p1 <- pr_point(ggplot(mb, aes(factor(depth), med))) +
    geom_point(aes(y = obs), colour = PAL$obs, shape = 4, size = 2, stroke = 0.9) +
    facet_wrap(~species, nrow = 1) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = "sampling depth (m)", y = "share of replicates with reads", title = "metabarcoding") + theme_joint(9)
  plots <- list(p1)
  if (!is.null(e$qpcr)) {
    p2 <- pr_point(ggplot(e$qpcr, aes(factor(depth), med))) +
      geom_point(aes(y = obs), colour = PAL$obs, shape = 4, size = 2, stroke = 0.9) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "sampling depth (m)", y = "detection rate", title = "hake qPCR detection") + theme_joint(9)
    ct <- e$ct; inside <- mean(ct$obs >= ct$lo & ct$obs <= ct$hi)
    p3 <- ggplot(ct, aes(pred, obs)) +
      geom_linerange(aes(ymin = lo, ymax = hi), colour = PAL$model, alpha = 0.15, linewidth = 0.4) +
      geom_abline(colour = PAL$muted, linetype = 2) +
      geom_point(colour = PAL$obs, size = 0.7, alpha = 0.7) +
      labs(x = "predicted mean Ct", y = "observed Ct",
           title = sprintf("hake Ct (%.0f%% inside the 90%% predictive interval)", 100 * inside)) + theme_joint(9)
    plots <- c(plots, list(p2 + p3))
  }
  wrap_plots(plots, ncol = 1) +
    plot_annotation(title = title, subtitle = "Blue: predictive (median, 90% over draws); orange: observed",
                    theme = theme(plot.title = element_text(face = "bold", colour = PAL$ink, size = 12),
                                  plot.subtitle = element_text(colour = PAL$muted, size = 9)))
}

# ---- Prior predictive --------------------------------------------------------------
# r_prior: the fit_joint(mode = "prior") summary; r_post (optional): a fit
# whose predictive totals are drawn for comparison
plot_prior_pred <- function(r_prior, r_post = NULL, title = "Prior predictive") {
  plots <- list()
  pv <- r_prior$ppc$visual
  if (!is.null(pv)) {
    d <- do.call(rbind, lapply(names(pv), function(s) rbind(
      data.frame(species = SP_LABEL[[s]], src = "prior", total = pv[[s]]$counts$rep$total),
      if (!is.null(r_post$ppc$visual)) data.frame(species = SP_LABEL[[s]], src = "posterior",
                                                  total = r_post$ppc$visual[[s]]$counts$rep$total))))
    o <- data.frame(species = SP_LABEL[names(pv)], total = vapply(pv, function(v) v$counts$obs$total, 0))
    plots$vis <- ggplot(d, aes(total + 1, fill = src)) +
      geom_histogram(bins = 30, colour = "white", linewidth = 0.2, alpha = 0.7, position = "identity") +
      geom_vline(data = o, aes(xintercept = total + 1), colour = PAL$obs, linewidth = 0.9) +
      scale_x_log10() + scale_fill_manual(values = c(prior = PAL$muted, posterior = PAL$model)) +
      facet_wrap(~species, scales = "free") +
      labs(x = "groups detected in the survey + 1 (log scale)", y = "draws", fill = NULL,
           title = "total detections") + theme_joint(9)
  }
  pe <- r_prior$ppc$edna
  if (!is.null(pe$qpcr)) {
    q <- pe$qpcr
    plots$q <- pr_point(ggplot(q, aes(factor(depth), med)), PAL$muted) +
      geom_point(aes(y = obs), colour = PAL$obs, shape = 4, size = 2, stroke = 0.9) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "sampling depth (m)", y = "detection rate", title = "hake qPCR detection (prior: median, 90%)") + theme_joint(9)
    mb <- pe$mb; mb$species <- factor(SP_LABEL[mb$species], SP_LABEL)
    plots$mb <- pr_point(ggplot(mb, aes(factor(depth), med)), PAL$muted) +
      geom_point(aes(y = obs), colour = PAL$obs, shape = 4, size = 2, stroke = 0.9) +
      facet_wrap(~species, nrow = 1) + coord_cartesian(ylim = c(0, 1)) +
      labs(x = "sampling depth (m)", y = "share with reads", title = "metabarcoding (prior)") + theme_joint(9)
  }
  wrap_plots(plots, ncol = 1) +
    plot_annotation(title = title, subtitle = "Grey: prior predictive (median, 90% over draws); orange: observed",
                    theme = theme(plot.title = element_text(face = "bold", colour = PAL$ink, size = 12),
                                  plot.subtitle = element_text(colour = PAL$muted, size = 9)))
}

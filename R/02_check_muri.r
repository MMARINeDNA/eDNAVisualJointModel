# =============================================================================
# 02_check_muri.r
#
# Checks + summary table for the formatted MURI data (R/01_load_muri.r), and a
# smoke test that the existing Stan-data formatters run. No models are fitted.
# Run from the repo root:  source("R/02_check_muri.r")
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
source("R/functions_visual.R")

edna <- readRDS("data/processed/muri_edna.rds")
vis  <- readRDS("data/processed/muri_visual.rds")
sam  <- edna$design$samples; ob <- edna$observed; seg <- vis$design$segments
hdr  <- function(x) cat("\n====", x, "====\n")

hdr("formatters (full objects)")
fe <- format_stan_data_gp2d_bathysp(edna)
cat("eDNA: N =", fe$stan_data$N, " N_qpcr_long =", fe$stan_data$N_qpcr_long,
    " N_mb_long =", fe$stan_data$N_mb_long, " K_bathy =", fe$stan_data$K_bathy,
    " M =", fe$stan_data$M, "\n")
for (sp in c("humpback", "pwsd")) {
  fv <- format_stan_data_visual(vis, sp)
  cat(sp, ": n_seg =", fv$n_seg, " n =", fv$n, " sum(seg_count) =", sum(fv$seg_count),
      " M =", fv$M, " K_bathy =", fv$K_bathy, " coords range =",
      paste(round(range(fv$coords), 2), collapse = " .. "), "\n")
}

hdr("summary table")
mbpos <- function(sp) mean(ob$mb_reads[, sp] > 0)
cat(sprintf("eDNA   bottles %d | stations %d | qPCR reps %d | MB reps %d\n",
            nrow(sam), n_distinct(sam$station), length(ob$qpcr_detect), nrow(ob$mb_reads)))
cat(sprintf("       hake qPCR detection rate %.3f | MB reps with reads: hake %.3f humpback %.3f pwsd %.3f\n",
            mean(ob$qpcr_detect), mbpos("hake"), mbpos("humpback"), mbpos("pwsd")))
cat(sprintf("visual segments %d | transects(sections) %d | effort %.1f km | sightings humpback %d, pwsd %d\n",
            nrow(seg), n_distinct(seg$transect_id), sum(seg$seg_l),
            nrow(vis$observed$humpback$obs), nrow(vis$observed$pwsd$obs)))
print(vis$meta$sight_flow)

hdr("domain box")
cat("eDNA bottles outside box:", sum(!sam$in_domain), "of", nrow(sam),
    "| X range", paste(round(range(sam$X)), collapse = ".."),
    "Y range", paste(round(range(sam$Y)), collapse = ".."), "\n")
cat("segments outside box:", sum(!seg$in_domain), "of", nrow(seg), "(",
    round(sum(seg$seg_l[!seg$in_domain])), "km of", round(sum(seg$seg_l)), "km )\n")
print(seg |> mutate(side = case_when(X < 0 ~ "X<0 (west)", X > 500 ~ "X>500 (east)",
                                     Y < 0 ~ "Y<0 (south)", Y > 1270 ~ "Y>1270 (north)",
                                     TRUE ~ "inside")) |> count(side))
for (sp in names(vis$observed)) {
  o <- vis$observed[[sp]]$obs
  cat(sp, "sightings in segments outside box:", sum(!seg$in_domain[o$seg_id]), "of", nrow(o), "\n")
}

hdr("bathymetry")
cat("segments Z_bathy NA:", sum(is.na(seg$Z_bathy)), " range:",
    paste(round(range(seg$Z_bathy, na.rm = TRUE)), collapse = ".."), "\n")
cat("bottles  Z_bathy NA:", sum(is.na(sam$Z_bathy)), " range:",
    paste(round(range(sam$Z_bathy, na.rm = TRUE)), collapse = ".."), "\n")
cat("bottles with Z_sample > Z_bathy:", sum(sam$Z_sample > sam$Z_bathy, na.rm = TRUE), "\n")
# compare with the survey's own BATH field at sightings
sg <- read.csv("data/sightings.csv")
bath <- do.call(rbind, lapply(seq_len(nrow(sg)), function(i) {
  j <- tryCatch(jsonlite::fromJSON(sg$oceano[i]), error = function(e) NULL)
  data.frame(S2004 = if (!is.null(j$BATH$S2004)) -j$BATH$S2004 else NA_real_,
             ETOPO1 = if (!is.null(j$BATH$ETOPO1)) -j$BATH$ETOPO1 else NA_real_)
}))
elev <- terra::extract(terra::rast(marmap::as.raster(readRDS("data/processed/etopo2022_60s_cce.rds"))),
                       cbind(sg$longitude, sg$latitude), method = "bilinear")[[1]]
bath$mine <- ifelse(elev < 0, -elev, NA_real_)
ok <- complete.cases(bath) & bath$ETOPO1 > 0
cat("sightings n =", sum(ok), " cor(mine, BATH$ETOPO1) =", round(cor(bath$mine[ok], bath$ETOPO1[ok]), 4),
    " median |diff| (m) =", round(median(abs(bath$mine - bath$ETOPO1)[ok], na.rm = TRUE), 1),
    " | vs S2004:", round(median(abs(bath$mine - bath$S2004)[ok], na.rm = TRUE), 1), "\n")

hdr("perpendicular distance (km) by species, before truncation")
sg2 <- sg |> mutate(sp = suppressWarnings(as.integer(species_name)))
si <- read.csv("data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv")
for (nm in c(humpback = 76, pwsd = 22)) {
  d <- si$pdist[si$species == nm]
  cat(names(nm), "n =", length(d), " quantiles(50/80/90/95/99/100) =",
      paste(round(quantile(d, c(.5, .8, .9, .95, .99, 1)), 2), collapse = " / "),
      "| share > 5 km:", round(mean(d > 5), 3), "\n")
}

hdr("group sizes: data/grpsz pools vs real best group sizes (obs)")
for (sp in c("humpback", "pwsd")) {
  pool <- vis$truth$sp_params[[sp]]$group_size_pool; real <- vis$observed[[sp]]$obs$size
  cat(sprintf("%-8s pool n=%d mean=%.2f median=%g q90=%g max=%g | real n=%d mean=%.2f median=%g q90=%g max=%g\n",
              sp, length(pool), mean(pool), median(pool), quantile(pool, .9), max(pool),
              length(real), mean(real), median(real), quantile(real, .9), max(real)))
}

hdr("sighting conditions / mixed species")
for (sp in names(vis$observed)) {
  o <- vis$observed[[sp]]$obs
  cat(sp, ": mixed =", sum(o$mixed %in% TRUE), "; bft NA =", sum(is.na(o$bft)),
      "; vis NA =", sum(is.na(o$vis)), "\n")
}
print(table(seg$mode, seg$efftype)); print(table(seg$eswsides))
cat("effort dates:", format(range(seg$datetime_begin)), "\n")
cat("sighting dates:", format(range(c(vis$observed$humpback$obs$datetime,
                                      vis$observed$pwsd$obs$datetime))), "\n")

hdr("eDNA details")
print(table(Z_sample = sam$Z_sample))
cat("bottles with >1 qPCR sample id:", sum(sam$n_qpcr_sample_ids > 1), "\n")
print(table(qpcr_dilution = ob$qpcr_dilution))
print(table(mb_dilution = ob$mb_dilution))
print(table(reps_per_bottle_qpcr = table(ob$qpcr_sample_idx)))
print(table(reps_per_bottle_mb = table(ob$mb_sample_idx)))
cat("qPCR inhibition_rate: summary\n"); print(summary(ob$qpcr_inhibition))
cat("share of qPCR reps with inhibition_rate > 0.5:", round(mean(ob$qpcr_inhibition > 0.5), 3), "\n")
cat("MB total reads quantiles:\n"); print(quantile(ob$mb_total, c(0, .01, .05, .5, .95, 1)))
cat("MB reps with total < 1000:", sum(ob$mb_total < 1000), "\n")

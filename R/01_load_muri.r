# =============================================================================
# 01_load_muri.r
#
# Format the real MURI data into the same list structures the simulators
# return, so the existing formatters (format_stan_data_gp2d_bathysp(),
# format_stan_data_visual()) run on real data unchanged.
#
#   eDNA   <- simulate_bathysp() shape   -> data/processed/muri_edna.rds
#   visual <- simulate_visual()  shape   -> data/processed/muri_visual.rds
#
# Run from the repo root:  source("R/01_load_muri.r")
# Every assumption and every value that could not be filled is listed in
# notes/muri_data_formatting.md. NO model is fitted here.
#
# Inputs (data/):
#   hake_qPCR_MURI_df.csv          qPCR, hake, copies/uL per replicate
#   MV1_MURI_df.csv                metabarcoding reads, long (taxon x PCR rep)
#   effort.csv, sightings.csv      OBIS-SEAMAP dataset 2147 (SWFSC CCES 2018, cruise 1651)
#   SITEINFO_CCE_2018_Pp_Lo_Mn.csv per-sighting perpendicular distance (pdist, km)
# Bathymetry: NOAA ETOPO 2022 60-arc-second bedrock, downloaded once via marmap
#   and cached in data/processed/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(sf); library(data.table)
})
source("R/functions_visual.R")   # LT_DOMAIN, lt_species_params(), HSGP4eDNA functions

W_TRUNC      <- 5.0              # km; kept at the simulations' value until decided
UTM_EPSG     <- 32610L
X_OFFSET_M   <- 100000; Y_OFFSET_M <- 4180000
OUT_DIR      <- "data/processed"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- helpers ----------------------------------------------------------------
lonlat_to_model <- function(lon, lat) {
  pts <- sf::st_transform(sf::st_as_sf(data.frame(lon = lon, lat = lat),
                                       coords = c("lon", "lat"), crs = 4326), UTM_EPSG)
  xy <- sf::st_coordinates(pts)
  tibble(X_utm = xy[, 1], Y_utm = xy[, 2],
         X = (xy[, 1] - X_OFFSET_M) / 1000, Y = (xy[, 2] - Y_OFFSET_M) / 1000)
}
in_domain <- function(X, Y) X >= 0 & X <= LT_DOMAIN$X_km_max & Y >= 0 & Y <= LT_DOMAIN$Y_km_max

# Bottom depth (m, positive down), bilinear from ETOPO 2022 60s; NA on land.
BATHY_SOURCE <- "NOAA ETOPO 2022 v1 60-arc-second bedrock elevation (via marmap::getNOAA.bathy), bilinear"
bottom_depth <- function(lon, lat) {
  cache <- file.path(OUT_DIR, "etopo2022_60s_cce.rds")
  if (!file.exists(cache)) {
    b <- marmap::getNOAA.bathy(lon1 = -131, lon2 = -113, lat1 = 26, lat2 = 52,
                               resolution = 1, keep = FALSE)
    saveRDS(b, cache)
  }
  r <- terra::rast(marmap::as.raster(readRDS(cache)))
  elev <- terra::extract(r, cbind(lon, lat), method = "bilinear")[[1]]
  ifelse(elev < 0, -elev, NA_real_)
}

# =============================================================================
# eDNA
# =============================================================================
q  <- read.csv("data/hake_qPCR_MURI_df.csv", stringsAsFactors = FALSE)
mb <- data.table::fread("data/MV1_MURI_df.csv")

# One bottle = one location_id (one station x one depth; verified: each
# location_id has a single lat/lon/depth and is identical across both files).
stopifnot(all(tapply(q$depth, q$location_id, function(x) length(unique(x))) == 1))
bottle <- unique(mb[, .(location_id, lat, lon, depth)]) |> as_tibble() |>
  arrange(location_id)
stopifnot(!anyDuplicated(bottle$location_id), setequal(bottle$location_id, q$location_id))

qsamp <- q |> group_by(location_id) |>
  summarise(qpcr_sample_ids = paste(sort(unique(sample)), collapse = ";"),
            n_qpcr_sample_ids = n_distinct(sample), .groups = "drop")

sample_depths <- sort(unique(bottle$depth))
stn_key <- paste(round(bottle$lat, 5), round(bottle$lon, 5))
xy <- lonlat_to_model(bottle$lon, bottle$lat)

samples <- bottle |>
  transmute(location_id, lat, lon, Z_sample = as.numeric(depth),
            station = match(stn_key, unique(stn_key))) |>
  bind_cols(xy) |>
  mutate(Z_bathy   = bottom_depth(lon, lat),
         depth_idx = match(Z_sample, sample_depths),
         sample_id = row_number(),
         in_domain = in_domain(X, Y),
         cruise = NA_character_, datetime = as.POSIXct(NA)) |>   # not in the files
  left_join(qsamp, by = "location_id") |>
  select(station, X, Y, X_utm, Y_utm, Z_bathy, Z_sample, depth_idx, sample_id,
         location_id, qpcr_sample_ids, n_qpcr_sample_ids, lat, lon, datetime, cruise, in_domain)
N <- nrow(samples)

# ---- qPCR (hake): one entry per replicate row of hake_qPCR_MURI_df.csv ------
qpcr_sample_idx <- match(q$location_id, samples$location_id)
stopifnot(!anyNA(qpcr_sample_idx), all(q$detected == as.integer(!is.na(q$hake_copies_ul))))

# ---- metabarcoding: wide reads per PCR replicate -----------------------------
SP_TARGET <- c(hake = "Merluccius productus", humpback = "Megaptera novaeangliae",
               pwsd = "Lagenorhynchus obliquidens")
mbrep <- unique(mb[, .(pcr_replicate_id, location_id, Dilution, Rep)])[order(pcr_replicate_id)]
tot   <- mb[, .(total = sum(Nreads)), by = pcr_replicate_id]
tgt   <- dcast(mb[species %in% SP_TARGET], pcr_replicate_id ~ species, value.var = "Nreads", fill = 0L)
mbrep <- Reduce(function(a, b) merge(a, b, by = "pcr_replicate_id", all.x = TRUE),
                list(mbrep, tot, tgt))
mb_reads <- as.matrix(mbrep[, ..SP_TARGET]); colnames(mb_reads) <- names(SP_TARGET)
mb_reads <- cbind(mb_reads, junk = mbrep$total - rowSums(mb_reads))
storage.mode(mb_reads) <- "integer"
stopifnot(all(mb_reads >= 0L), all(rowSums(mb_reads) == mbrep$total))
# Genus-level "Lagenorhynchus" reads stay in junk (see notes).

edna <- list(
  meta = list(
    n_species = 3L,
    sp_names  = c(hake = "Merluccius_productus", humpback = "Megaptera_novaeangliae",
                  pwsd = "Lagenorhynchus_obliquidens"),
    sp_common = c("Pacific hake", "Humpback whale", "Pacific white-sided dolphin"),
    conv_factor  = c(hake = 10, humpback = 200, pwsd = 110, junk = 10),  # PLACEHOLDER (simulator)
    vol_filtered = 2.5,    # PLACEHOLDER (simulator): litres; real volume not in the data
    vol_aliquot  = 2,      # PLACEHOLDER (simulator): uL template; real value not in the data
    vol_elution  = NA_real_,  # uL; model assumes 100
    N = N, N_qpcr_long = nrow(q), N_mb_long = nrow(mbrep),
    X_km_max = LT_DOMAIN$X_km_max, Y_km_max = LT_DOMAIN$Y_km_max,
    placeholders = c("conv_factor", "vol_filtered", "vol_aliquot",
                     "truth$qpcr_params (all NA: no Ct / standard curve available)"),
    ct_available = FALSE,
    bathy_source = BATHY_SOURCE,
    sample_depths = sample_depths
  ),
  design = list(samples = samples),
  observed = list(
    qpcr_sample_idx = qpcr_sample_idx,
    qpcr_detect     = as.integer(q$detected),
    qpcr_ct         = rep(NA_real_, nrow(q)),     # Ct NOT available
    mb_sample_idx   = match(mbrep$location_id, samples$location_id),
    mb_reads        = mb_reads,
    mb_total        = as.integer(rowSums(mb_reads)),
    # extras (ignored by the formatter): per-replicate provenance
    qpcr_copies_ul = q$hake_copies_ul, qpcr_dilution = q$dilution,
    qpcr_inhibition = q$inhibition_rate, qpcr_plate = q$qPCR, qpcr_sample_id = q$sample,
    mb_pcr_replicate_id = mbrep$pcr_replicate_id, mb_dilution = mbrep$Dilution,
    mb_rep = mbrep$Rep
  ),
  truth = list(
    zsample_effect = matrix(1, N, 3, dimnames = list(NULL, names(SP_TARGET))),
    qpcr_params = list(kappa = NA_real_, alpha_ct = NA_real_, beta_ct = NA_real_,
                       gamma_0 = NA_real_, gamma_1 = NA_real_, sigma_0 = NA_real_),
    bathy_spline = NULL
  )
)
stopifnot(!anyNA(edna$observed$mb_sample_idx))

# =============================================================================
# Visual
# =============================================================================
SPECIES_CODE <- c(humpback = 76L, pwsd = 22L)   # SWFSC species codes (Mn, Lo)

ef <- read.csv("data/effort.csv", stringsAsFactors = FALSE) |>
  mutate(t_begin = as.POSIXct(datetime_begin, tz = "UTC"),
         t_end   = as.POSIXct(datetime_end,   tz = "UTC"),
         t_mid   = as.POSIXct(mdatetime,      tz = "UTC")) |>
  arrange(t_begin)
stopifnot(all(ef$ds_type == "lneff"))   # on-effort line-transect rows only

xy_seg <- lonlat_to_model(ef$mlon, ef$mlat)
segments <- bind_cols(
  tibble(seg_id = seq_len(nrow(ef)), transect_id = as.integer(factor(ef$section_id, unique(ef$section_id))),
         X = xy_seg$X, Y = xy_seg$Y, seg_l = ef$length_km,
         Z_bathy = bottom_depth(ef$mlon, ef$mlat)),
  tibble(X_utm = xy_seg$X_utm, Y_utm = xy_seg$Y_utm, lat = ef$mlat, lon = ef$mlon,
         datetime_begin = ef$t_begin, datetime_end = ef$t_end, datetime_mid = ef$t_mid,
         cruise = ef$cruise, section_id = ef$section_id, segnum = ef$segnum, row_id = ef$row_id,
         mode = ef$mode, efftype = ef$efftype, eswsides = ef$eswsides,
         speed_kph = ef$speed_kph, avgbft = ef$avgbft, avgvis = ef$avgvis,
         avgswellhght = ef$avgswellhght,
         in_domain = in_domain(xy_seg$X, xy_seg$Y))
)

# ---- sightings --------------------------------------------------------------
sg <- read.csv("data/sightings.csv", stringsAsFactors = FALSE) |>
  mutate(sp_code = suppressWarnings(as.integer(species_name)),
         t = as.POSIXct(datetime, tz = "UTC")) |>
  filter(sp_code %in% SPECIES_CODE)
si <- read.csv("data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv", stringsAsFactors = FALSE) |>
  transmute(sightno = snum, sp_code = species, pdist)
stopifnot(!anyDuplicated(si[, c("sightno", "sp_code")]))

# On effort := sighting time lies inside an effort segment's [begin, end].
seg_of <- vapply(sg$t, function(x) {
  w <- which(ef$t_begin <= x & x <= ef$t_end); if (length(w)) w[1] else NA_integer_
}, integer(1))
sg <- sg |> mutate(seg_id = seg_of) |>
  left_join(si, by = c("sightno", "sp_code")) |>
  mutate(size = as.integer(round(gsspbest)))

sight_flow <- sg |> group_by(sp_code) |>
  summarise(n_all = n(), n_on_effort = sum(!is.na(seg_id)),
            n_on_effort_with_pdist = sum(!is.na(seg_id) & !is.na(pdist)),
            n_on_effort_pdist_missing = sum(!is.na(seg_id) & is.na(pdist)),
            n_size_missing = sum(!is.na(seg_id) & !is.na(pdist) & is.na(size)),
            n_beyond_w = sum(!is.na(seg_id) & !is.na(pdist) & pdist > W_TRUNC),
            .groups = "drop")

make_obs <- function(code) {
  d <- sg |> filter(sp_code == code, !is.na(seg_id), !is.na(pdist), !is.na(size),
                    pdist <= W_TRUNC) |> arrange(seg_id, t)
  data.frame(seg_id = as.integer(d$seg_id), distance = d$pdist, size = d$size,
             sightno = d$sightno, datetime = d$t, lat = d$latitude, lon = d$longitude,
             bft = d$bft, vis = d$vis, swellhght = d$swellhght,
             mixed = d$mixed, nsp = d$nsp, row_id = d$row_id)
}
observed_vis <- lapply(SPECIES_CODE, function(code) {
  o <- make_obs(code)
  list(obs = o, seg_count = tabulate(o$seg_id, nbins = nrow(segments)))
})

visual <- list(
  meta = list(seed = NA, species = names(SPECIES_CODE), w = W_TRUNC,
              domain = LT_DOMAIN, bathy_source = BATHY_SOURCE,
              pdist_source = "data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv (pdist, km), joined on sightno x species",
              sight_flow = sight_flow),
  design = list(segments = as.data.frame(segments)),
  observed = observed_vis,
  truth = list(sp_params = lt_species_params()[names(SPECIES_CODE)], spline_setup = NULL)
)

# ---- in-domain subsets (same structure; seg_id renumbered) ------------------
subset_visual <- function(v) {
  seg <- v$design$segments; keep <- which(seg$in_domain)
  map <- setNames(seq_along(keep), keep)
  v$design$segments <- transform(seg[keep, ], seg_id = seq_along(keep))
  for (sp in names(v$observed)) {
    o <- v$observed[[sp]]$obs; o <- o[o$seg_id %in% keep, , drop = FALSE]
    o$seg_id <- as.integer(map[as.character(o$seg_id)])
    v$observed[[sp]] <- list(obs = o, seg_count = tabulate(o$seg_id, nbins = length(keep)))
  }
  v
}
subset_edna <- function(e) {
  s <- e$design$samples; keep <- which(s$in_domain)
  map <- setNames(seq_along(keep), keep)
  e$design$samples <- transform(s[keep, ], sample_id = seq_along(keep))
  ob <- e$observed
  qk <- ob$qpcr_sample_idx %in% keep; mk <- ob$mb_sample_idx %in% keep
  for (nm in grep("^qpcr_", names(ob), value = TRUE)) ob[[nm]] <- ob[[nm]][qk]
  for (nm in grep("^mb_", names(ob), value = TRUE))
    ob[[nm]] <- if (is.matrix(ob[[nm]])) ob[[nm]][mk, , drop = FALSE] else ob[[nm]][mk]
  ob$qpcr_sample_idx <- as.integer(map[as.character(ob$qpcr_sample_idx)])
  ob$mb_sample_idx   <- as.integer(map[as.character(ob$mb_sample_idx)])
  e$observed <- ob
  e$truth$zsample_effect <- e$truth$zsample_effect[keep, , drop = FALSE]
  e$meta$N <- length(keep); e$meta$N_qpcr_long <- length(ob$qpcr_sample_idx)
  e$meta$N_mb_long <- length(ob$mb_sample_idx)
  e
}

saveRDS(edna,   file.path(OUT_DIR, "muri_edna.rds"))
saveRDS(visual, file.path(OUT_DIR, "muri_visual.rds"))
saveRDS(subset_edna(edna),     file.path(OUT_DIR, "muri_edna_indomain.rds"))
saveRDS(subset_visual(visual), file.path(OUT_DIR, "muri_visual_indomain.rds"))
message("wrote ", OUT_DIR, "/muri_{edna,visual}{,_indomain}.rds")

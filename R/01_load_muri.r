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
# Study area: north of 35°N (LAT_MIN). Effort segments with a midpoint south of
# it, and their sightings, are dropped; all eDNA bottles are north of it. The
# domain box is HSGP4eDNA's STUDY_DOMAIN (checked against the data).
#
# Inputs. Private raw files (git-ignored) are read from RAW_DIR (env var;
# default data/); the tracked ones from data/:
#   Hake_2019_metadata.csv          [RAW_DIR] 2019 eDNA bottles: sample number,
#                                   station, Niskin, depth, volume, date, flags
#   M3_compiled_taxon_table_wide.csv [RAW_DIR] metabarcoding reads, taxon x PCR
#                                   replicate, several markers (MV1 used)
#   qPCR_2019A.rds                  [RAW_DIR] the qPCR lab's 2019 analysis A:
#                                   hake wells (Ct, plate, dilution, ...),
#                                   standards, controls; depth_cat per sample
#   data/processed/qpcr_standard_curves.rds  per-plate hake standard curves
#                                   (run R/fit_qpcr_standards.r first)
#   effort.csv, sightings.csv       OBIS-SEAMAP dataset 2147 (SWFSC CCES 2018, cruise 1651)
#   SITEINFO_CCE_2018_Pp_Lo_Mn.csv  [RAW_DIR] per-sighting perpendicular distance (pdist, km)
# Env: RAW_DIR; ETOH ("drop", default, or "keep": the bottles washed with 30%
#   ethanol by mistake, flagged 30EtOH / 30EtOHpaired).
# Bathymetry: NOAA ETOPO 2022 60-arc-second bedrock, downloaded once via marmap
#   and cached in data/processed/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(sf); library(data.table)
})
source("R/functions_visual.R")   # lt_species_params(), HSGP4eDNA functions

TRUNC_QUANTILE <- 0.95           # truncation: drop the largest 5% of distances, per species
UTM_EPSG     <- STUDY_DOMAIN$utm_epsg   # 32610 (HSGP4eDNA R/functions.R)
LAT_MIN      <- STUDY_DOMAIN$lat_min    # 35°N: every analysis is restricted to north of this
BFT_MAX      <- 5                       # visual effort: segments with mean Beaufort above this are
                                        # dropped (and their sightings); the domain is unaffected
DOMAIN_BUFFER_KM <- 10           # margin around the data (km)
DOMAIN_ROUND_KM  <- 10           # origin and extent are multiples of this (km)
RAW_DIR      <- Sys.getenv("RAW_DIR", "data")
raw <- function(f) {
  p <- file.path(RAW_DIR, f)
  if (!file.exists(p)) stop("missing ", p, " (private raw file; set RAW_DIR to its folder)")
  p
}
ETOH         <- match.arg(Sys.getenv("ETOH", "drop"), c("drop", "keep"))
MB_MARKER    <- "MV1"                     # MarVer1; the compiled table also has MFU, C16, DL, ...
# Taxa dropped from the junk column AND the replicate total (contamination;
# exact for the remaining categories). Human and cattle are already absent.
CONTAMINANTS <- c("Gallus", "Gallus gallus", "Homo sapiens", "Bos taurus", "Sus scrofa")
OUT_DIR      <- "data/processed"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- helpers ----------------------------------------------------------------
lonlat_to_utm <- function(lon, lat) {
  pts <- sf::st_transform(sf::st_as_sf(data.frame(lon = lon, lat = lat),
                                       coords = c("lon", "lat"), crs = 4326), UTM_EPSG)
  xy <- sf::st_coordinates(pts)
  cbind(E = xy[, 1], N = xy[, 2])
}

# ---- model domain: the study area north of LAT_MIN ----------------------------
# X = (easting - E0)/1000, Y = (northing - N0)/1000 (km), with (E0, N0) the
# south-west corner of the box, so the domain is [0, X_km_max] x [0, Y_km_max]
# as the downstream code assumes. Only effort segments whose midpoint is at or
# north of LAT_MIN are used (and their sightings); every eDNA bottle is north of
# it. The box covers the kept segments (begin, middle and end points) and every
# bottle, plus DOMAIN_BUFFER_KM to the west, east and north; its southern edge
# is LAT_MIN at the zone's central meridian (123°W), rounded down to the
# DOMAIN_ROUND_KM grid. It must equal HSGP4eDNA's STUDY_DOMAIN, which the
# simulators use, so simulated and real fits share one domain.
.mv   <- sub("^([^-]+-){2}(\\d+)-(\\d+)-.*$", "\\2-\\3",
             grep(sprintf("^[^-]+-%s-", MB_MARKER), names(fread(raw("M3_compiled_taxon_table_wide.csv"), nrows = 0)), value = TRUE))
.q    <- fread(raw("Hake_2019_metadata.csv"))[sampleID %in% .mv]   # the bottles with metabarcoding
.ef   <- read.csv("data/effort.csv", stringsAsFactors = FALSE)
stopifnot(nrow(.q) > 0, all(.q$lat >= LAT_MIN))       # all bottles are north of 35°N
.ef   <- .ef[.ef$mlat >= LAT_MIN, ]
.pts  <- rbind(lonlat_to_utm(.q$lon, .q$lat),
               lonlat_to_utm(.ef$longitude_begin, .ef$latitude_begin),
               lonlat_to_utm(.ef$longitude_end,   .ef$latitude_end),
               lonlat_to_utm(.ef$mlon,            .ef$mlat))
.step <- DOMAIN_ROUND_KM * 1000; .buf <- DOMAIN_BUFFER_KM * 1000
E0 <- floor((min(.pts[, "E"]) - .buf) / .step) * .step
N0 <- floor(min(.pts[, "N"], lonlat_to_utm(-123, LAT_MIN)[, "N"]) / .step) * .step
DOMAIN <- list(X_km_max = ceiling((max(.pts[, "E"]) + .buf - E0) / .step) * DOMAIN_ROUND_KM,
               Y_km_max = ceiling((max(.pts[, "N"]) + .buf - N0) / .step) * DOMAIN_ROUND_KM)
if (!isTRUE(all.equal(c(E0, N0, DOMAIN$X_km_max, DOMAIN$Y_km_max),
                      unlist(STUDY_DOMAIN[c("E0", "N0", "X_km_max", "Y_km_max")], use.names = FALSE))))
  stop("the data's domain (E0 ", E0, ", N0 ", N0, ", ", DOMAIN$X_km_max, " x ", DOMAIN$Y_km_max,
       " km) differs from HSGP4eDNA's STUDY_DOMAIN; update STUDY_DOMAIN there")
rm(.q, .ef, .pts, .mv)

lonlat_to_model <- function(lon, lat) {
  en <- lonlat_to_utm(lon, lat)
  tibble(X_utm = en[, "E"], Y_utm = en[, "N"],
         X = (en[, "E"] - E0) / 1000, Y = (en[, "N"] - N0) / 1000)
}
in_domain <- function(X, Y) X >= 0 & X <= DOMAIN$X_km_max & Y >= 0 & Y <= DOMAIN$Y_km_max

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
# eDNA (2019): one bottle = one sample number. Each station x depth usually has
# two bottles (Niskin biological replicates); they are separate samples at the
# same location. Bottles, qPCR wells and metabarcoding replicates are matched
# on the sample number (sampleID "52193-<n>"; "52193_<n>" in the qPCR file).
# =============================================================================
meta <- fread(raw("Hake_2019_metadata.csv"))[!is.na(sample)]
stopifnot(!anyDuplicated(meta$sampleID))
A  <- readRDS(raw("qPCR_2019A.rds"))
qs <- as.data.table(A$samples)[species == "hake"][, sid := sub("_", "-", sample)]

# ---- metabarcoding columns: <run>-<marker>-<project>-<sample>-d<fold>-<rep>[-sr<k>][_S<n>]
m3 <- fread(raw("M3_compiled_taxon_table_wide.csv"))
RE <- "^([A-Za-z0-9]+)-([A-Za-z0-9]+)-(\\d+)-(\\d+)-d(\\d+)-(\\d+)(-sr\\d+)?(_S\\d+)?$"
stopifnot(all(grepl(RE, names(m3)[-(1:2)])))
mbcol <- data.table(col = names(m3)[-(1:2)])
mbcol[, `:=`(run = sub(RE, "\\1", col), marker = sub(RE, "\\2", col),
             sampleID = paste(sub(RE, "\\3", col), sub(RE, "\\4", col), sep = "-"),
             fold = as.integer(sub(RE, "\\5", col)), rep = as.integer(sub(RE, "\\6", col)))]
mbcol <- mbcol[marker == MB_MARKER]
tax <- m3$BestTaxon
SP_TARGET <- c(hake = "Merluccius productus", humpback = "Megaptera novaeangliae",
               pwsd = "Lagenorhynchus obliquidens")
stopifnot(all(SP_TARGET %in% tax), !anyDuplicated(tax))
cnt <- as.matrix(m3[, mbcol$col, with = FALSE]); cnt[is.na(cnt)] <- 0
rownames(cnt) <- tax
mbcol[, `:=`(total = colSums(cnt), contam = colSums(cnt[tax %in% CONTAMINANTS, , drop = FALSE]))]

# ---- bottles: metabarcoded samples, not controls, north of LAT_MIN ----------
flow <- list(mb_samples = uniqueN(mbcol$sampleID))
bot <- meta[sampleID %in% mbcol$sampleID & control == "no" & is.na(field.negative.type)]
flow$not_control <- nrow(bot)
bot[, etoh := !is.na(drop.sample) & grepl("EtOH", drop.sample)]
flow$etoh <- sum(bot$etoh)
if (ETOH == "drop") bot <- bot[etoh == FALSE]
bot <- bot[lat >= LAT_MIN]
flow$kept <- nrow(bot)

# Target sampling depth (depth_cat, m) from the qPCR file, one per sample
dcat <- unique(qs[, .(sid, depth_cat)])
stopifnot(!anyDuplicated(dcat$sid), all(bot$sampleID %in% dcat$sid))
bot[, depth_cat := dcat$depth_cat[match(sampleID, dcat$sid)]]
DEPTH_LEVELS <- sort(unique(bot$depth_cat))

stn_key <- paste(round(bot$lat, 5), round(bot$lon, 5))
xy <- lonlat_to_model(bot$lon, bot$lat)
samples <- tibble(
  sample_id = seq_len(nrow(bot)), sampleID = bot$sampleID,
  station = match(stn_key, unique(stn_key)), station_name = bot$station, niskin = bot$Niskin,
  lat = bot$lat, lon = bot$lon) |>
  bind_cols(xy) |>
  # bottom depth: the metadata's bottom.depth.consensus (the lab's choice);
  # ETOPO at the bottle kept for comparison (segments use ETOPO)
  mutate(Z_bathy = bot$bottom.depth.consensus,
         Z_bathy_etopo = bottom_depth(lon, lat),
         # target sampling depth: an ordered category (m); the recorded depth is kept
         Z_sample  = as.numeric(bot$depth_cat),
         depth_cat = factor(bot$depth_cat, levels = DEPTH_LEVELS, ordered = TRUE),
         depth_idx = as.integer(depth_cat),
         depth_recorded = bot$depth,
         location = match(paste(station, depth_cat), unique(paste(station, depth_cat))),
         vol_filtered = bot$volume,                        # litres, per bottle
         date = as.Date(bot$date, format = "%m/%d/%Y"),
         etoh_wash = bot$etoh, cruise = "52193")   # project ID: the 2019 hake survey
stopifnot(all(in_domain(samples$X, samples$Y)), !anyNA(samples$vol_filtered), !anyNA(samples$date))
N <- nrow(samples)

# ---- qPCR (hake): every analysis-A well of the kept bottles ------------------
# Analysis A is cleaned upstream (no inhibited wells; one dilution per sample;
# outliers dealt with), so all wells are used as given. Non-detects carry
# Ct = -99 in the file.
q <- qs[sid %in% samples$sampleID]
stopifnot(all(q[, uniqueN(dilution), by = sid]$V1 == 1), all(q$Ct_bin %in% 0:1),
          all(q$Ct[q$Ct_bin == 1] > 0))
curves_f <- file.path(OUT_DIR, "qpcr_standard_curves.rds")
if (!file.exists(curves_f)) stop("run R/fit_qpcr_standards.r first (writes ", curves_f, ")")
sc <- readRDS(curves_f)
stopifnot(all(q$qPCR %in% sc$curves$plate))
q_plate <- match(q$qPCR, sc$curves$plate)
shared  <- setNames(sc$shared$mean, sc$shared$variable)

# ---- metabarcoding: per bottle, the dilution with the most reads ------------
mb <- mbcol[sampleID %in% samples$sampleID & total - contam > 0]
flow$mb_reps_with_reads <- nrow(mb)
# the same sample x dilution x replicate sequenced twice: keep the larger run
mb <- mb[order(-total)][!duplicated(mb[order(-total), .(sampleID, fold, rep)])]
best <- mb[, .(reads = sum(total - contam)), by = .(sampleID, fold)][order(-reads)][!duplicated(sampleID)]
flow$mb_bottles_several_dilutions <- sum(mb[, uniqueN(fold), by = sampleID]$V1 > 1)
mb <- mb[best, on = .(sampleID, fold), nomatch = NULL][order(match(sampleID, samples$sampleID), rep)]
flow$mb_reps_kept <- nrow(mb)
# bottles whose metabarcoding replicates all had no reads keep their qPCR
flow$bottles_without_mb <- sum(!samples$sampleID %in% mb$sampleID)
mb_reads <- t(cnt[SP_TARGET, mb$col, drop = FALSE]); colnames(mb_reads) <- names(SP_TARGET)
mb_reads <- cbind(mb_reads, junk = mb$total - mb$contam - rowSums(mb_reads))
storage.mode(mb_reads) <- "integer"
stopifnot(all(mb_reads >= 0L), all(rowSums(mb_reads) == mb$total - mb$contam))

edna <- list(
  meta = list(
    n_species = 3L,
    sp_names  = c(hake = "Merluccius_productus", humpback = "Megaptera_novaeangliae",
                  pwsd = "Lagenorhynchus_obliquidens"),
    sp_common = c("Pacific hake", "Humpback whale", "Pacific white-sided dolphin"),
    conv_factor  = c(hake = 10, humpback = 200, pwsd = 110, junk = 10),  # PLACEHOLDER (simulator)
    # The eDNA model (HSGP4eDNA v1.5) takes one filtered volume; per-bottle
    # volumes are in design$samples$vol_filtered (2.5 L for almost every bottle).
    vol_filtered = stats::median(samples$vol_filtered),
    vol_aliquot  = 2,      # uL of template per reaction: qPCR (qPCR vignette) and
                           # MarVer1 PCRs (confirmed by the lab, 2026-10-22)
    vol_elution  = 100,    # uL, assumed (HSGP4eDNA's value); a constant, so it only
                           # rescales conv_factor
    N = N, N_qpcr_long = nrow(q), N_mb_long = nrow(mb),
    X_km_max = DOMAIN$X_km_max, Y_km_max = DOMAIN$Y_km_max, lat_min = LAT_MIN,
    utm_origin = c(E0 = E0, N0 = N0), utm_epsg = UTM_EPSG,
    placeholders = c("conv_factor", "elution volume (assumed 100 uL)"),
    ct_available = TRUE,
    qpcr_source = "qPCR_2019A.rds, hake, analysis A",
    mb_source = sprintf("M3_compiled_taxon_table_wide.csv, marker %s; per bottle the dilution with the most reads", MB_MARKER),
    mb_contaminants_dropped = CONTAMINANTS,
    etoh = ETOH, flow = flow,
    bathy_source = "Hake_2019_metadata.csv bottom.depth.consensus (ETOPO in Z_bathy_etopo)",
    sample_depths = DEPTH_LEVELS, depth_levels = DEPTH_LEVELS
  ),
  design = list(samples = samples),
  observed = list(
    qpcr_sample_idx = match(q$sid, samples$sampleID),
    qpcr_detect     = as.integer(q$Ct_bin),
    qpcr_ct         = ifelse(q$Ct_bin == 1, q$Ct, NA_real_),
    mb_sample_idx   = match(mb$sampleID, samples$sampleID),
    mb_reads        = mb_reads,
    mb_total        = as.integer(rowSums(mb_reads)),
    # extras (not used by the v1.5 formatter): per-well and per-replicate detail
    qpcr_plate = q$qPCR, qpcr_plate_idx = q_plate, qpcr_well = q$well,
    qpcr_dilution = q$dilution, qpcr_log_dilution = log(q$dilution),   # fraction: 0.2 = 1:5
    qpcr_inhibition = q$inhibition_rate, qpcr_ipc_ct = q$IPC_Ct, qpcr_copies_ul = q$copies_ul,
    qpcr_alpha_ct = sc$curves$alpha_ct[q_plate], qpcr_beta_ct = sc$curves$beta_ct[q_plate],
    mb_column = mb$col, mb_run = mb$run, mb_rep = mb$rep,
    mb_dilution = mb$fold, mb_log_dilution = -log(mb$fold),             # fold: 5 = 1:5
    mb_contaminant_reads = as.integer(mb$contam)
  ),
  truth = list(
    zsample_effect = matrix(1, N, 3, dimnames = list(NULL, names(SP_TARGET))),
    # One curve for the v1.5 model: the mean over plates of the per-plate
    # curves, with the shared noise and detection terms. Per-plate curves:
    # qpcr_curves and observed$qpcr_alpha_ct / qpcr_beta_ct.
    qpcr_params = list(alpha_ct = mean(sc$curves$alpha_ct), beta_ct = mean(sc$curves$beta_ct),
                       kappa = shared[["kappa"]], gamma_0 = shared[["gamma0"]],
                       gamma_1 = shared[["gamma1"]], sigma_0 = shared[["sigma0"]]),
    qpcr_curves = sc$curves,
    bathy_spline = NULL
  )
)
stopifnot(!anyNA(edna$observed$mb_sample_idx), !anyNA(edna$observed$qpcr_sample_idx))
message(sprintf("eDNA: %d bottles (%d stations, depths %s m); %d qPCR wells on %d plates; %d metabarcoding replicates",
                N, max(samples$station), paste(DEPTH_LEVELS, collapse = "/"), nrow(q),
                length(unique(q$qPCR)), nrow(mb)))

# =============================================================================
# Visual
# =============================================================================
SPECIES_CODE <- c(humpback = 76L, pwsd = 22L)   # SWFSC species codes (Mn, Lo)

ef_all <- read.csv("data/effort.csv", stringsAsFactors = FALSE) |>
  mutate(t_begin = as.POSIXct(datetime_begin, tz = "UTC"),
         t_end   = as.POSIXct(datetime_end,   tz = "UTC"),
         t_mid   = as.POSIXct(mdatetime,      tz = "UTC")) |>
  arrange(t_begin)
stopifnot(all(ef_all$ds_type == "lneff"))   # on-effort line-transect rows only
# Study area: segments whose midpoint is at or north of LAT_MIN, with mean
# Beaufort at most BFT_MAX (the usual line-transect limit; decided 2026-10-22).
# The domain above is built from all effort north of LAT_MIN, so it does not move.
ef <- ef_all |> filter(mlat >= LAT_MIN, avgbft <= BFT_MAX)

xy_seg <- lonlat_to_model(ef$mlon, ef$mlat)
segments <- bind_cols(
  tibble(seg_id = seq_len(nrow(ef)), transect_id = as.integer(factor(ef$section_id, unique(ef$section_id))),
         X = xy_seg$X, Y = xy_seg$Y, seg_l = ef$length_km,
         Z_bathy = bottom_depth(ef$mlon, ef$mlat)),
  tibble(X_utm = xy_seg$X_utm, Y_utm = xy_seg$Y_utm, lat = ef$mlat, lon = ef$mlon,
         datetime_begin = ef$t_begin, datetime_end = ef$t_end, datetime_mid = ef$t_mid,
         cruise = ef$cruise, section_id = ef$section_id, segnum = ef$segnum, row_id = ef$row_id,
         mode = ef$mode, efftype = ef$efftype, eswsides = ef$eswsides,
         speed_kph = ef$speed_kph, bft = ef$avgbft, vis = ef$avgvis, avgbft = ef$avgbft, avgvis = ef$avgvis,
         avgswellhght = ef$avgswellhght)
)
stopifnot(all(in_domain(xy_seg$X, xy_seg$Y)), all(in_domain(samples$X, samples$Y)))

# ---- sightings --------------------------------------------------------------
sg <- read.csv("data/sightings.csv", stringsAsFactors = FALSE) |>
  mutate(sp_code = suppressWarnings(as.integer(species_name)),
         t = as.POSIXct(datetime, tz = "UTC")) |>
  filter(sp_code %in% SPECIES_CODE)
si <- read.csv(raw("SITEINFO_CCE_2018_Pp_Lo_Mn.csv"), stringsAsFactors = FALSE) |>
  transmute(sightno = snum, sp_code = species, pdist)
stopifnot(!anyDuplicated(si[, c("sightno", "sp_code")]))

# On effort := sighting time lies inside an effort segment's [begin, end].
# Matched against ALL effort, then sightings on segments south of LAT_MIN or
# above BFT_MAX are dropped (counted in sight_flow as n_south_of_lat_min and
# n_bft_above_max).
seg_of_all <- vapply(sg$t, function(x) {
  w <- which(ef_all$t_begin <= x & x <= ef_all$t_end); if (length(w)) w[1] else NA_integer_
}, integer(1))
south <- !is.na(seg_of_all) & ef_all$mlat[seg_of_all] < LAT_MIN
rough <- !is.na(seg_of_all) & !south & ef_all$avgbft[seg_of_all] > BFT_MAX
seg_of <- match(ef_all$row_id[seg_of_all], ef$row_id)      # index among kept segments
sg <- sg |> mutate(seg_id = seg_of, south_of_lat_min = south, bft_above_max = rough) |>
  left_join(si, by = c("sightno", "sp_code")) |>
  mutate(size = as.integer(round(gsspbest)))

# Truncation distance per species = 95th percentile of the on-effort
# perpendicular distances that have a group size (5% of detections dropped).
# `obs` holds ALL such sightings; format_stan_data_visual() truncates at
# meta$w[[species]] and rebuilds seg_count.
cand <- sg |> filter(!is.na(seg_id), !is.na(pdist), !is.na(size))
W_SP <- vapply(SPECIES_CODE, function(code)
  as.numeric(quantile(cand$pdist[cand$sp_code == code], TRUNC_QUANTILE)), numeric(1))

sight_flow <- sg |> group_by(sp_code) |>
  summarise(n_all = n(), n_south_of_lat_min = sum(south_of_lat_min),
            n_bft_above_max = sum(bft_above_max),
            n_on_effort = sum(!is.na(seg_id)),
            n_on_effort_with_pdist = sum(!is.na(seg_id) & !is.na(pdist)),
            n_on_effort_pdist_missing = sum(!is.na(seg_id) & is.na(pdist)),
            n_size_missing = sum(!is.na(seg_id) & !is.na(pdist) & is.na(size)),
            .groups = "drop") |>
  mutate(w = W_SP[match(sp_code, SPECIES_CODE)],
         n_used = n_on_effort_with_pdist - n_size_missing,
         n_beyond_w = vapply(seq_len(n()), function(i)
           sum(cand$pdist[cand$sp_code == sp_code[i]] > w[i]), integer(1)))

make_obs <- function(code) {
  d <- cand |> filter(sp_code == code) |> arrange(seg_id, t)
  data.frame(seg_id = as.integer(d$seg_id), distance = d$pdist, size = d$size,
             sightno = d$sightno, datetime = d$t, lat = d$latitude, lon = d$longitude,
             bft = d$bft, vis = d$vis, swellhght = d$swellhght,
             mixed = d$mixed, nsp = d$nsp, row_id = d$row_id)
}
observed_vis <- lapply(SPECIES_CODE, function(code) {
  o <- make_obs(code)
  list(obs = o, seg_count = tabulate(o$seg_id, nbins = nrow(segments)))
})

# Detection-function settings per species (exploratory AIC on this data set,
# forward selection, segment-level covariates; see notes/muri_data_formatting.md):
#   humpback  hazard-rate key; group size only (no Beaufort); use_size_covar must be set here because
#             lt_species_params() has it off for humpback (the simulation truth)
#   pwsd      half-normal key; group size + Beaufort + visibility (the AIC set; visibility
#             added 2026-10-21, ROADMAP_2026-10-21 section 5 item 5)
# Beaufort is the segment-level `bft` column and visibility the `vis` column (segment
# means, binned by the formatter). Visibility is the observers' estimate of the distance
# in nautical miles at which a dolphin could be seen surfacing against the water (SWFSC
# DAS format, "DasAll Events and Associated Data Fields", event W field 5; the swfscDAS
# package's inst/DAS_Format.pdf).
# Segment counts: negative binomial for both species (use_nb; decided 2026-10-21,
# notes/count_dispersion_check.md: harmless on Poisson data, needed if sightings cluster).
sp_params <- lt_species_params()[names(SPECIES_CODE)]
sp_params$humpback$det_key        <- 1L
sp_params$humpback$use_size_covar <- 1L   # group size selected by AIC (the simulator's humpback truth has none)
sp_params$humpback$use_bft        <- FALSE
sp_params$pwsd$use_bft      <- TRUE
sp_params$pwsd$use_vis      <- TRUE
sp_params$humpback$use_nb <- 1L
sp_params$pwsd$use_nb     <- 1L

visual <- list(
  meta = list(seed = NA, species = names(SPECIES_CODE), w = W_SP,
              trunc_quantile = TRUNC_QUANTILE,
              domain = DOMAIN, lat_min = LAT_MIN, bft_max = BFT_MAX,
              utm_origin = c(E0 = E0, N0 = N0), utm_epsg = UTM_EPSG,
              bathy_source = BATHY_SOURCE,
              pdist_source = "data/SITEINFO_CCE_2018_Pp_Lo_Mn.csv (pdist, km), joined on sightno x species",
              sight_flow = sight_flow),
  design = list(segments = as.data.frame(segments)),
  observed = observed_vis,
  truth = list(sp_params = sp_params, spline_setup = NULL)
)

saveRDS(edna,   file.path(OUT_DIR, "muri_edna.rds"))
saveRDS(visual, file.path(OUT_DIR, "muri_visual.rds"))
message("wrote ", OUT_DIR, "/muri_{edna,visual}.rds; domain ", DOMAIN$X_km_max, " x ",
        DOMAIN$Y_km_max, " km, origin (E0, N0) = (", E0, ", ", N0, ")")

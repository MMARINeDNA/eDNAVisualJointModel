# =============================================================================
# make_pred_grid.r
#
# The real-data prediction grid for the joint model's maps: cell centres every
# SPACING km over the 35N domain (the processed objects' UTM origin), kept
# where ETOPO 2022 gives water, with the bottom depth there (bilinear, as
# R/01_load_muri.r) and each cell's distance to the nearest bottle or segment.
#
#   Rscript R/make_pred_grid.r
# Env: SPACING (km, default 20)
# Needs sf, terra, marmap and data/processed/etopo2022_60s_cce.rds (the
# loader's cache). Out: data/processed/pred_grid_real_<SPACING>km.rds
# =============================================================================
source("R/functions_joint.R")
SPACING <- as.numeric(Sys.getenv("SPACING", "20"))
edna <- readRDS("data/processed/muri_edna.rds"); vis <- readRDS("data/processed/muri_visual.rds")
dom <- vis$meta$domain
g <- pred_grid_domain(dom, SPACING)
E0 <- vis$meta$utm_origin[["E0"]]; N0 <- vis$meta$utm_origin[["N0"]]
pts <- sf::st_as_sf(data.frame(E = E0 + 1000 * g$X, N = N0 + 1000 * g$Y), coords = c("E", "N"),
                    crs = vis$meta$utm_epsg)
ll <- sf::st_coordinates(sf::st_transform(pts, 4326))
r <- terra::rast(marmap::as.raster(readRDS("data/processed/etopo2022_60s_cce.rds")))
elev <- terra::extract(r, ll, method = "bilinear")[[1]]
g$lon <- ll[, 1]; g$lat <- ll[, 2]; g$Z_bathy <- ifelse(elev < 0, -elev, NA_real_)
# keep the open ocean: water cells connected (edge-adjacent) to the domain's
# western edge, which drops San Francisco Bay and lakes
ix <- match(g$X, sort(unique(g$X))); iy <- match(g$Y, sort(unique(g$Y)))
wet <- matrix(FALSE, max(ix), max(iy)); wet[cbind(ix, iy)] <- !is.na(g$Z_bathy)
reach <- matrix(FALSE, nrow(wet), ncol(wet)); reach[1, ] <- wet[1, ]
repeat {
  nb <- reach
  nb[-1, ] <- nb[-1, ] | reach[-nrow(reach), ]; nb[-nrow(nb), ] <- nb[-nrow(nb), ] | reach[-1, ]
  nb[, -1] <- nb[, -1] | reach[, -ncol(reach)]; nb[, -ncol(nb)] <- nb[, -ncol(nb)] | reach[, -1]
  nb <- nb & wet
  if (sum(nb) == sum(reach)) break
  reach <- nb
}
cat(sprintf("water cells %d, connected to the open ocean %d\n", sum(wet), sum(reach)))
g <- g[reach[cbind(ix, iy)], ]
# and the Salish Sea (Strait of Juan de Fuca inward, Puget Sound), which joins the ocean
# at this resolution
g <- g[!((g$lon > -124.4 & g$lat > 47.9) | (g$lon > -123.4 & g$lat > 46.9)), ]
xy <- rbind(as.matrix(edna$design$samples[, c("X", "Y")]), as.matrix(vis$design$segments[, c("X", "Y")]))
g$dist_km <- apply(sqrt(outer(g$X, xy[, 1], "-")^2 + outer(g$Y, xy[, 2], "-")^2), 1, min)
rownames(g) <- NULL
attr(g, "spacing") <- SPACING
out <- sprintf("data/processed/pred_grid_real_%gkm.rds", SPACING)
saveRDS(g, out)
cat(sprintf("%d water cells of %g km (of %d); depth %.0f-%.0f m; saved %s\n", nrow(g), SPACING,
            nrow(pred_grid_domain(dom, SPACING)), min(g$Z_bathy), max(g$Z_bathy), out))

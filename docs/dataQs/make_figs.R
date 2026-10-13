# =============================================================================
# docs/dataQs/make_figs.R
#
# Numbers and figures for docs/dataQs.md (qPCR vs metabarcoding data in hand).
#
#   Rscript docs/dataQs/make_figs.R [path/to/qPCR_2019A.rds]
#
# Inputs: data/hake_qPCR_MURI_df.csv, data/MV1_MURI_df.csv, qPCR_2019A.rds
# (default data/qPCR_2019A.rds). Writes docs/dataQs/fig_*.png and prints the
# numbers quoted in the document.
# =============================================================================

suppressMessages({ library(dplyr); library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
f19  <- if (length(args)) args[1] else "data/qPCR_2019A.rds"
out  <- "docs/dataQs"

q  <- read.csv("data/hake_qPCR_MURI_df.csv", stringsAsFactors = FALSE) |>
  mutate(sample = as.character(sample))
mb <- fread("data/MV1_MURI_df.csv")
A  <- readRDS(f19)
xs <- A$samples |> filter(species == "hake") |> mutate(sample = sub("^52193_", "", sample))
st <- A$standards |> filter(species == "hake")

# A MURI well is "in 2019A" if sample, plate and inhibition rate all match.
key <- function(d) paste(d$sample, d$qPCR, round(d$inhibition_rate, 3))
q$in2019A <- key(q) %in% key(xs)
xm <- xs |> filter(sample %in% q$sample)

# ---- numbers -----------------------------------------------------------------
mbr <- unique(mb[, .(pcr_replicate_id, location_id, Dilution, Rep)])
cat("bottles (location_id):", n_distinct(q$location_id), " qPCR sample IDs:", n_distinct(q$sample), "\n")
cat("bottles with 2 qPCR sample IDs:", sum(table(unique(q[, c("location_id", "sample")])$location_id) == 2), "\n")
cat("qPCR wells (MURI csv):", nrow(q), "  MB PCR replicates:", nrow(mbr), "\n")
cat("MB reps per bottle:"); print(table(table(mbr$location_id)))
cat("MB template dilution:"); print(table(mbr$Dilution))
cat("qPCR wells per sample x plate x dilution (MURI):"); print(table(table(paste(q$sample, q$qPCR, q$dilution))))
cat("qPCR wells per sample (2019A, our samples):"); print(table(table(xm$sample)))
cat("MURI wells by dilution x in2019A:\n"); print(table(dilution = q$dilution, in2019A = q$in2019A))
rr <- unique(q$sample[q$dilution < 1 & !q$in2019A])
cat("samples with diluted reruns:", length(rr), "\n")
cat("inhibition, undiluted wells of rerun samples: median",
    median(q$inhibition_rate[q$sample %in% rr & q$dilution == 1]),
    " share > 0.75", mean(q$inhibition_rate[q$sample %in% rr & q$dilution == 1] > 0.75), "\n")
cat("inhibition, 0.1/0.5 reruns: median", median(q$inhibition_rate[q$dilution %in% c(0.1, 0.5)]), "\n")
cat("inhibition, 0.2 wells kept in 2019A: median",
    median(q$inhibition_rate[q$dilution == 0.2 & q$in2019A]), "\n")
mbmax <- mbr[, .(mb_dil = max(Dilution) > 1), by = location_id]
b <- q |> group_by(location_id) |> summarise(qpcr_rerun = any(dilution < 1), .groups = "drop") |>
  left_join(as.data.frame(mbmax), by = "location_id")
cat("qPCR rerun x MB diluted (bottles):\n"); print(table(qpcr_rerun = b$qpcr_rerun, mb_diluted = b$mb_dil))
cat("2019A, our samples: drop.sample:"); print(table(xm$drop.sample, useNA = "ifany"))
cat("  samples flagged:", n_distinct(xm$sample[!is.na(xm$drop.sample)]), "\n")
cat("  dilution:"); print(table(xm$dilution))
cat("  detected wells with Ct > 37.5:", sum(xm$Ct > 37.5), " of detected", sum(xm$Ct > 0),
    "; max Ct", max(xm$Ct), "\n")
cat("hake standards:", nrow(st), "wells on", n_distinct(st$qPCR), "plates; MURI plates without standards:",
    length(setdiff(unique(q$qPCR), unique(st$qPCR))), "; inhibited (> 0.75):", sum(st$inhibition_rate > 0.75),
    "; non-detects:", sum(st$Ct == -99), "\n")
cat("dates:", format(range(xm$date)), "\n")

# ---- Figure 1: what happens to one bottle -----------------------------------
box <- function(x, y, w, h, label, fill) list(
  annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
           fill = fill, colour = "grey30", linewidth = 0.3),
  annotate("text", x = x, y = y, label = label, size = 3, lineheight = 0.95))
arr <- function(x1, y1, x2, y2, lty = "solid") annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
  arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.4, linetype = lty, colour = "grey25")

p1 <- ggplot() + xlim(0, 10.2) + ylim(0, 6.2) + theme_void() +
  box(1.2, 3.1, 2.0, 1.3, "Water bottle\n(2.5 L filtered)\n551 bottles", "grey92") +
  box(3.8, 3.1, 2.0, 1.3, "DNA extract\n(100 uL)", "grey92") +
  box(7.6, 5.2, 4.4, 1.3, "qPCR (hake): 3 wells x 2 uL, undiluted;\nif inhibited (moa IPC): rerun diluted 1:5\nor undiluted on a later plate (kept in 2019A)", "#cfe3f3") +
  box(7.6, 3.3, 4.4, 0.9, "further reruns at 1:2 and 1:10\n(in MURI csv, not in 2019A)", "#f3dccf") +
  box(7.6, 1.2, 4.4, 1.3, "Metabarcoding (MV1): 3 PCR reps\n(Rep 1-3; some 4-6)\ntemplate diluted 1x to 100x", "#d7ecd2") +
  arr(2.2, 3.1, 2.8, 3.1) + arr(4.8, 3.4, 5.35, 4.9) + arr(4.8, 2.8, 5.35, 1.3) +
  arr(7.6, 4.55, 7.6, 3.75, "dashed") +
  annotate("text", x = 5, y = 6.1, size = 3.2, fontface = "bold",
           label = "Per bottle: about 3 qPCR wells and 3 metabarcoding PCRs (~1:1, not 1:3)")
ggsave(file.path(out, "fig_workflow.png"), p1, width = 7, height = 3.6, dpi = 200, bg = "white")

# ---- Figure 2: MURI qPCR wells by plate, dilution and 2019A status ----------
pl <- q |> mutate(plate_num = as.integer(sub("^H", "", qPCR)),
                  status = case_when(in2019A ~ "in 2019A",
                                     dilution == 1 ~ "not in 2019A: undiluted, inhibited",
                                     TRUE ~ "not in 2019A: diluted rerun"),
                  dil = factor(paste0("dilution ", dilution), levels = paste0("dilution ", c(1, 0.5, 0.2, 0.1))))
p2 <- ggplot(pl, aes(reorder(qPCR, plate_num), fill = status)) + geom_bar(width = 0.8) +
  facet_grid(dil ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("in 2019A" = "#4a8cc2", "not in 2019A: diluted rerun" = "#e08a4a",
                               "not in 2019A: undiluted, inhibited" = "#b23a3a")) +
  labs(x = "qPCR plate", y = "hake qPCR wells", fill = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5),
                                  strip.text.y = element_text(angle = 0))
ggsave(file.path(out, "fig_wells_by_plate.png"), p2, width = 7, height = 4.6, dpi = 200, bg = "white")
cat("figures written to", out, "\n")

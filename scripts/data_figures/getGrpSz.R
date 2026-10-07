# Build the empirical group-size pools used by the distance-sampling
# simulator (R/functions_visual.R: lt_species_params()). Run from the project root.
# Outputs are small and tracked: data/grpsz/{humpback,pwsd}.rds

library(dplyr)

dir.create("./data/grpsz", showWarnings = FALSE, recursive = TRUE)

sightings <- read.csv("./data/sightings.csv")

si.humpback <- sightings %>%
  filter(spcode == 76) %>%
  select(group_size)

si.humpback <- na.omit(si.humpback$group_size)

saveRDS(si.humpback, file = "./data/grpsz/humpback.rds")

si.pwsd <- sightings %>%
  filter(spcode == 22) %>%
  select(group_size)

si.pwsd <- na.omit(si.pwsd$group_size)

saveRDS(si.pwsd, file = "./data/grpsz/pwsd.rds")

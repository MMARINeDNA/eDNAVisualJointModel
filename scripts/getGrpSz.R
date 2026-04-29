library(dplyr)

si.humpback <- read.csv("./data/sightings.csv") %>%
  filter(spcode == 76) %>%
  select(group_size)

si.humpback <- na.omit(si.humpback$group_size)

save(si.humpback, file = "./distance/humpback_grpsz.RData")

si.pwsd <- read.csv("./data/sightings.csv") %>%
  filter(spcode == 22) %>%
  select(group_size) 

si.pwsd <- na.omit(si.pwsd$group_size)

save(si.pwsd, file = "./distance/pwsd_grpsz.RData")

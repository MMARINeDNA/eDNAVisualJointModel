library(dplyr)
library(ggplot2)
library(usmap)
library(sf)

load("./ProcessedData/detect_data.Rdata")

pos.data <- filter(detect_data, Detected == 1) %>%
  st_as_sf(coords = c("lon", "lat"))

st_crs(pos.data) <- 4326

pos.data$Type <- "eDNA"

map.data <- us_map(include = c("OR", "WA", "CA"))

map.sf <- st_transform(map.data, 4326)

map.crop <- st_crop(map.sf, st_bbox(pos.data))

##### Harbour porpoise

pp.sf <- filter(pos.data, BestTaxon == "Phocoena phocoena" )

pp.crop <- st_crop(map.sf, st_bbox(pp.sf))

ggplot() +
  geom_sf(data = pp.crop) +
  geom_sf(data = pp.sf, aes(color = primer)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw()

###### Pacific white-sided dolphin

lags.sf <- filter(pos.data, BestTaxon == "Lagenorhynchus obliquidens")

lags.sf <- lags.sf %>% select(BestTaxon, geometry, Type) %>%
  rename(Species = BestTaxon) 

lags.vis <- filter(dat_sf, Species == "Lagenorhynchus obliquidens") %>%
  select(Species, geometry, Type) %>%
  st_crop(st_bbox(pos.data))

lags <- rbind(lags.sf, lags.vis)

ggplot() +
  geom_sf(data = map.crop) +
  geom_sf(data = lags, aes(color = Type)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() +
  facet_wrap(~Type)

#### Beaked whales


beaker.sf <- filter(pos.data, BestTaxon %in% c("Ziphius cavirostris","Berardius bairdii"))

beaker.sf <- beaker.sf %>% select(BestTaxon, geometry, Type) %>%
  rename(Species = BestTaxon) 

beaker.vis <- filter(dat_sf, Species %in% c("Ziphius cavirostris", "ziphiid whale", "Mesoplodon sp.", "Berardius bairdii")) %>%
  select(Species, geometry, Type) %>%
  st_crop(st_bbox(pos.data))

beaker <- rbind(beaker.sf, beaker.vis)

ggplot() +
  geom_sf(data = map.crop) +
  geom_sf(data = beaker, aes(color = Type)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() +
  facet_wrap(~Type)

##### Humpback whales

mn.sf <- filter(pos.data, BestTaxon %in% c("Megaptera novaeangliae"))

mn.sf <- mn.sf %>% select(BestTaxon, geometry, Type) %>%
  rename(Species = BestTaxon) 

mn.vis <- filter(dat_sf, Species %in% c("Megaptera novaeangliae")) %>%
  select(Species, geometry, Type) %>%
  st_crop(st_bbox(pos.data))

mn <- rbind(mn.sf, mn.vis)

ggplot() +
  geom_sf(data = map.crop) +
  geom_sf(data = mn, aes(color = Type)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() +
  facet_wrap(~Type)

##### Dall's Porpoise

pd.sf <- filter(pos.data, BestTaxon %in% c("Phocoenoides dalli"))

pd.sf <- pd.sf %>% select(BestTaxon, geometry, Type) %>%
  rename(Species = BestTaxon) 

pd.vis <- filter(dat_sf, Species %in% c("Phocoenoides dalli")) %>%
  select(Species, geometry, Type) %>%
  st_crop(st_bbox(pos.data))

pd <- rbind(pd.sf, pd.vis)

ggplot() +
  geom_sf(data = map.crop) +
  geom_sf(data = pd, aes(color = Type)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() +
  facet_wrap(~Type)


##### Fin whale

bp.sf <- filter(pos.data, BestTaxon %in% c("Balaenoptera physalus"))

bp.sf <- bp.sf %>% select(BestTaxon, geometry, Type) %>%
  rename(Species = BestTaxon) 

bp.vis <- filter(dat_sf, Species %in% c("Balaenoptera physalus")) %>%
  select(Species, geometry, Type) %>%
  st_crop(st_bbox(pos.data))

bp <- rbind(bp.sf, bp.vis)

ggplot() +
  geom_sf(data = map.crop) +
  geom_sf(data = bp, aes(color = Type)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() +
  facet_wrap(~Type)


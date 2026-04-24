library(dplyr)
library(ggplot2)
library(usmap)
library(sf)

load("Data/detect_data.Rdata")

effort <- read.csv("Data/effort.csv")
sightings <- read.csv("Data/sightings.csv")

bp_sightings <- sightings %>%
  filter(common_name == "Fin Whale") %>%
  select(longitude, latitude, groupsize_count, ) %>%
  st_as_sf(coords = c("longitude", "latitude"))%>%
  st_set_crs(4326) %>%
  mutate(Type = "Visual", What = "Detections") %>%
  st_crop(st_bbox(stations))


lines_sf <- effort %>%
  select(longitude_begin, longitude_end, latitude_begin, latitude_end) %>%
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(matrix(c(longitude_begin, longitude_end, 
                                                  latitude_begin, latitude_end), ncol = 2)))) %>%
  st_as_sf() %>%
  st_set_crs(4326) %>%
  mutate(Type = "Visual", What = "Effort") %>%
  st_crop(st_bbox(stations))


stations <- detect_data %>% 
  select(lon, lat) %>%
  distinct() %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(4326) %>%
  mutate(Type = "Visual", What = "Effort")

map.data <- us_map(include = c("OR", "WA", "CA")) %>%
  st_transform(4326) %>%
  st_crop(st_bbox(stations))

bp_eDNA <- filter(detect_data, Detected == 1 & BestTaxon == "Balaenoptera physalus") %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(4326) %>%
  mutate(Type = "Visual", What = "Detections")

detections <- rbind.data.frame(select(bp_eDNA, geometry, Type, What),
                               select(bp_sightings, geometry, Type, What))


ggplot() +
  geom_sf(data = map.crop) +
  geom_sf(data = stations, color = "darkgrey", size = 0.5) +
#  geom_sf(data = lines_sf) +
  geom_sf(data = bp_eDNA, color = "#00BFC4") +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() 

ggplot() +
  geom_sf(data = map.crop) +
    geom_sf(data = lines_sf, color = "darkgrey") +
  geom_sf(data = bp_sightings, color = "#F8766D") +
  scale_x_continuous(n.breaks = 2) +
  theme_bw() 


# revisions for Biological Invasions submission
# Caitlin C. Mothes
# 03/08/24

library(tidyverse)
library(terra)
library(sf)
library(mapview)
library(spThin)
library(spdep)


# LANCOVER ASSOCIATION --------------------
## read in data -------------

## raw 10m landcover data
lc <- rast("data/landcover_mdc.tif")

## reclassify based on crosswalk doc

reclass <- read_csv("data/CLC_to_FNAINC_Crosswalk_20160701.csv")


land_reclass <- terra::classify(lc, as.matrix(reclass %>% dplyr::select(SITE, reclass)))
# values 2134 and 5300 in raster but not in table...
land_reclass[land_reclass == 5300] <- 16
land_reclass[land_reclass == 2134] <- 10

## read in occurrences

## COMBINE AND REMOVE DUPLICATES
### Agama
ag_all <- read_csv("data/agama_occ_mdc.csv") %>% select(LONG = Longitude, LAT = Latitude)

agama_pts <- read_csv("data/agama_sites.csv") %>% select(LONG, LAT) 

ag_combined <- bind_rows(ag_all, agama_pts) %>% distinct(LONG, LAT) %>% 
  mutate(species = "Agama")


### Basilisk
bas_all <- read_csv("data/basilisk_occ_mdc.csv") %>% select(LONG = Longitude, LAT = Latitude)

bas_pts <- read_csv("data/basilisk_sites.csv") %>% select(LONG, LAT)

bas_combined <- bind_rows(bas_all, bas_pts) %>% distinct(LONG, LAT) %>% 
 mutate(species = "Basilisk")

mapview::mapview(bas_combined)

## LC AGAMA ---------------

### JOIN COUNTS TEST ---------------------

## thin points (increasing thin.par as needed)
ag_thin <- thin(ag_combined, lat.col = "LAT", long.col = "LONG", spec.col = "species",
                thin.par = 0.25,
                reps = 50, write.files = FALSE, locs.thinned.list.return = TRUE)


## test for autocorrelation with join counts test, using last df
ag_thin_sf <- ag_thin[[50]] %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = st_crs(land_reclass))

## create raster for test
r_ag <- rast(extent = ext(ag_thin_sf), resolution = 100, crs = crs(land_reclass))

r_ag2 <- rasterize(ag_thin_sf, r_ag, field = 1)
# make NA values 0
r_ag2[is.na(r_ag2)] <- 0

# check how many points were lost to rasterization at slightly larger res
freq(r_ag2) # none

## set up connectivity matrix, rook default
nb <- cell2nb(nrow = nrow(r_ag2), ncol = ncol(r_ag2))

## convert to w weights list using Binary weights matrix
lwb <- nb2listw(nb, style = "B")

#join counts test and permuation test
ag_jc <- joincount.test(as.factor(values(r_ag2)), lwb, alternative = "greater")

ag_jc_perm <- joincount.mc(as.factor(values(r_ag2)), lwb, nsim = 999, alternative = "greater")

 ### re-calculate LC stats -----------
ag_lc <- ag_thin_sf %>% 
  mutate(landcover = as.numeric(extract(lc, .)[[2]])) %>%
  left_join(reclass, by = c("landcover" = "SITE")) %>% 
  st_drop_geometry() %>% 
  count(NAME_STATE) %>% 
  mutate(ag_percent = n/sum(n)*100,
         ag_point_estimate = (n+2)/(sum(n)+4),
         ag_L95 = ag_point_estimate-1.96*sqrt((ag_point_estimate*(1-ag_point_estimate))/(sum(n)+4)),
         ag_U95 = ag_point_estimate+1.96*sqrt((ag_point_estimate*(1-ag_point_estimate))/(sum(n)+4)))


## BASILISK LC -------------------

## thin points (increasing thin.par as needed)
bas_thin <- thin(bas_combined, lat.col = "LAT", long.col = "LONG", spec.col = "species",
                thin.par = 0.25,
                reps = 50, write.files = FALSE, locs.thinned.list.return = TRUE)


## test for autocorrelation with join counts test, using last df
bas_thin_sf <- bas_thin[[50]] %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = st_crs(land_reclass))

## create raster for test
r_bas <- rast(extent = ext(bas_thin_sf), resolution = 100, crs = crs(land_reclass))

r_bas2 <- rasterize(bas_thin_sf, r_bas, field = 1)
# make NA values 0
r_bas2[is.na(r_bas2)] <- 0

# check how many points were lost to rasterization at slightly larger res
freq(r_bas2) # only 2

## set up connectivity matrix, rook default
nb <- cell2nb(nrow = nrow(r_bas2), ncol = ncol(r_bas2))

## convert to w weights list using Binary weights matrix
lwb <- nb2listw(nb, style = "B")

#join counts test and permuation test
bas_jc <- joincount.test(as.factor(values(r_bas2)), lwb, alternative = "greater")

bas_jc_perm <- joincount.mc(as.factor(values(r_bas2)), lwb, nsim = 999, alternative = "greater")


## re-run landcover stats ---------------

bas_lc <- bas_thin_sf %>% 
  mutate(landcover = as.numeric(extract(lc, .)[[2]])) %>%
  left_join(reclass, by = c("landcover" = "SITE")) %>% 
  st_drop_geometry() %>% 
  count(NAME_STATE) %>% 
  mutate(bas_percent = n/sum(n)*100,
         bas_point_estimate = (n+2)/(sum(n)+4),
         bas_L95 = bas_point_estimate-1.96*sqrt((bas_point_estimate*(1-bas_point_estimate))/(sum(n)+4)),
         bas_U95 = bas_point_estimate+1.96*sqrt((bas_point_estimate*(1-bas_point_estimate))/(sum(n)+4)))


# set up ----------------------

library(raster)
library(rgdal)
library(tidyverse)


# landcover ---------------------

#path to geodatabase
gdb <- "CLC_v3_3_RASTER/CLC_v3_3_RASTER/CLC_v3_3_RASTER.gdb/"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
gdb_list <- ogrListLayers(gdb)
print(gdb_list)
#does not return any layers...

fc <- readOGR(dsn=gdb,layer="CLC_v3_3_SITE")

require(sf)
fc <- sf::st_read("CLC_v3_3_RASTER/CLC_v3_3_RASTER/CLC_v3_3_RASTER.gdb", layer = "CLC_v3_3_SITE")
fcstate <- sf::st_read("CLC_v3_3_RASTER/CLC_v3_3_RASTER/CLC_v3_3_RASTER.gdb", layer = "CLC_v3_3_STATE")


# try reading in SITE file exported from arcmap

landcover <- raster("D:/lizard_genomics/Land_vars/CLC_v3_3_SITE1.tif")

#huge dataset, clip to MDC

# read in mdc polygon

mdc <- readOGR("~/Desktop/NicheProject/MDC_Boundary.shp")

# set to same projection as landcover

mdc.prj <- spTransform(mdc, crs(landcover))

ext <- extent(mdc.prf)

land.crop <- crop(landcover, y = ext, snap="out") 

cropped <- setValues(land.crop, NA) 

mdc.raster <- rasterize(mdc.prj, cropped) 

predictors.masked <- mask(x=land.crop, mask=mdc.raster, 
                          filename = "D:/lizard_genomics/Land_vars/landcover_mdc.tif")

# load in mdc landcover

land.mdc <- raster("D:/lizard_genomics/land_vars/landcover_mdc.tif")

# reclassify based on crosswalk doc

reclass <- read.csv("data/CLC_to_FNAINC_Crosswalk_20160701.csv")


land.reclass <- raster::reclassify(land.mdc, as.matrix(reclass %>% dplyr::select(SITE, reclass)))
# values 2134 and 5300 in raster but not in table...
land.reclass[land.reclass == 5300] <- 16
land.reclass[land.reclass == 2134] <- 10

# make bounding box around points

basilisk_pts <- read.csv("data/basilisk_sites.csv") %>% filter(ID != "BV072")
# REMOVE spatial OUTLIERS BV072 (48), maybe BV056 (40)?


bas_sp <- SpatialPointsDataFrame(basilisk_pts[,c('LONG', 'LAT'),], data = basilisk_pts)

proj4string(bas_sp) <- CRS("+proj=longlat +datum=WGS84")


bas_sp <- spTransform(bas_sp, crs(land.100))


agama_pts <- read.csv("data/agama_sites.csv")

ag_sp <- SpatialPointsDataFrame(agama_pts[,c("LONG", "LAT"),], data = agama_pts)

proj4string(ag_sp) <-  CRS("+proj=longlat +datum=WGS84")

ag_sp <- spTransform(ag_sp, crs(land.reclass))

#make bounding box around points to clip landcover to

## add 3km to points extent

ext <- extent(xmin(bas_sp)+3000, xmax(bas_sp)+3000, ymin(bas_sp)+3000, ymax(bas_sp)+3000)


land.box <- crop(land.reclass, ext)

# get min dist between two points
bas_dist <- as.data.frame(bas_sp) %>% dplyr::select(LONG, LAT) %>% dist(.)
min(bas_dist)
#0.003616684  ~ 361m = .3-.4km

ag.dist <- as.data.frame(ag_sp) %>% dplyr::select(LONG,LAT) %>% dist(.)
min(ag.dist)
#0.004386035

#do 100m res has 720,291 cells
land.100 <- raster::aggregate(land.reclass, fact = 10, fun = modal, na.rm = TRUE,
                              filename = "D:/lizard_genomics/land_vars/land_mdc_100m.tif")

#do 250m res, has 115,520 cells
land.250 <- raster::aggregate(land.reclass, fact = 25, fun = modal, na.rm = TRUE)



# read in roads raster and clip to MDC

roads <- raster("D:/lizard_genomics/land_vars/traffic.tif")

# project to landcover crs

roads.prj <- projectRaster(roads, land.100)

# set non roads (i.e. NA) to 0 and then clip to mdc boundary (to remove ocean)

roads.prj[is.na(roads.prj)] <- 0


cropped <- setValues(roads.prj, NA) 

mdc.raster <- rasterize(mdc.prj, cropped) 

roads.mdc <- mask(x=roads.prj, mask=mdc.raster, 
                          filename = "D:/lizard_genomics/land_vars/roads_mdc.tif")

#canals variable


canals <- raster("D:/lizard_genomics/land_vars/canals.tif")

# project to landcover crs

#NEED TO CHANGE METHOD TO NGB FOR CAT VARIABLES
canals.prj <- projectRaster(canals, land.100, method = "ngb")

# set non canals to 0 and then clip to mdc boundary (to remove ocean)

canals.prj[canals.prj > 23] <- 0 #NA vlues became 28


#set all other values to 1

canals.prj[canals.prj != 0] <- 1

# crop to MDC boundary (re-do with edited canals layer)

cropped <- setValues(canals.prj, NA) 

mdc.raster <- rasterize(mdc_land.prj, cropped) 

canals.mdc <- mask(x=canals.prj, mask=mdc.raster, 
                  filename = "D:/lizard_genomics/land_vars/canals_mdc_update.tif")


# canopy cover

canopy <- raster("D:/WoodTurtleSDM/data/Predictors/current_new/canopy_prj.tif")

ext <- drawExtent()

canopy.crop <- crop(canopy, ext)

canopy.prj <- projectRaster(canopy.crop, land.100)

# remove water cells

canopy.prj[canopy.prj > 100] <- NA

#clip to MDC boundary (this will remove the rest of the ocean cells)
canopy.mdc <- mask(x=canopy.prj, mask=mdc.raster, 
                   filename = "D:/lizard_genomics/land_vars/canopy_mdc.tif")


#stack all vars

#read them in
land.100 <- raster("D:/lizard_genomics/land_vars/land_mdc_100m.tif")
roads.mdc <- raster("D:/lizard_genomics/land_vars/roads_mdc.tif")
canals.mdc <- raster("D:/lizard_genomics/land_vars/canals_mdc.tif")
canopy.mdc <- raster("D:/lizard_genomics/land_vars/canopy_mdc.tif")


liz_vars <- stack(land.100, roads.mdc, canals.mdc, canopy.mdc)


# mask to mdc_land to remove islands

mdc_land <- readOGR("MDC_land_boundary.shp")

mdc_land.prj <- spTransform(mdc_land, crs(land))


#rasterize
cropped <- setValues(land, NA) 

mdc_land.raster <- rasterize(mdc_land.prj, cropped) 

liz_vars.land <- mask(liz_vars, mask = mdc_land.raster)


#crop to extent around points

plot(liz_vars.land[[1]])


#make bounding box around points to clip landcover to

## add 2km to points extent

points(ag_sp, col = "blue")
points(bas_sp, col = "red") # removed outlier basilisk BV072 pt in everglades

#combine points
all_pts <- rbind(ag_sp, bas_sp)

ext <- extent(xmin(all_pts)-2000, xmax(all_pts)+2000, ymin(all_pts)-2000, ymax(all_pts)+2000)

liz.vars.crop <- crop(liz_vars.land, ext)

#edit land cover vars (remove or combine)

## combine canal/ditch (16) with artificial pond (14) = artificial waterways

liz.vars.crop[[1]][liz.vars.crop[[1]] == 16] <- 14

## combine 10 and 11 = wetlands in general

liz.vars.crop[[1]][liz.vars.crop[[1]] == 11] <- 10

#save variables

for (i in 1:nlayers(liz.vars.crop)){
  writeRaster(liz.vars.crop[[i]], paste0("data/", names(liz.vars.crop[[i]]), "_final.tif"), overwrite = TRUE)
}

#save as ascii
names(liz.vars.crop) <- c("land", "roads", "canals", "canopy")
for (i in 1:nlayers(liz.vars.crop)){
  writeRaster(liz.vars.crop[[i]], paste0("data/predictors_ascii/", names(liz.vars.crop[[i]]), ".asc"), overwrite = TRUE)
}


# save as ascii for resistanceGA

canals <- raster("data/predictors/canals_mdc_final.tif")

writeRaster(canals, "data/predictors_ascii/canals.asc")

canopy <- raster("data/predictors/canopy_mdc_final.tif")
writeRaster(canopy, "data/predictors_ascii/canopy.asc")

roads <- raster("data/predictors/roads_mdc_final.tif")
writeRaster(roads, "data/predictors_ascii/roads.asc")

land <- raster("data/predictors/land_mdc_100m_final.tif")
raster::writeRaster(land, "data/predictors_ascii/land.asc", format = "ascii", prj = TRUE, overwrite = TRUE)


# clip updated canals layer to extent of others

canals.update <- crop(canals.mdc, extent(canopy))
writeRaster(canals.update, "data/predictors/canals_mdc_finalUpdate.tif")
writeRaster(canals.update, "data/predictors/canals_mdc_finalUpdate.asc",
            format = "ascii", overwrite = TRUE) 

#create distance to canals layer and re-run SS_optim

canals <- raster("data/predictors/canals_mdc_finalUpdate.tif") 

canals[canals != 1] <- NA

dist2canals <- raster::distance(canals)
# mask to other layers to get NA in ocean
canals <- raster("data/predictors/canals_mdc_finalUpdate.tif") 

dist2canals[is.na(canals)] <- NA

writeRaster(dist2canals, "data/predictors/dist2canals.asc")
writeRaster(dist2canals, "data/predictors/dist2canals.tif")


#create roads as feature

roads <- raster("data/predictors/roads_mdc_final.tif")
roads[roads > 0] <- 1

writeRaster(roads, "data/predictors/roadsFeature.asc")


#create distance to roads layer
roadsFeature <- raster("data/predictors/roadsFeature.asc")

roadsFeature[roadsFeature != 1] <- NA

dist2road <- raster::distance(roadsFeature)
dist2road[is.na(roads)] <- NA

# check variation in point values
crs(dist2road) <- crs(roads)

raster::extract(dist2road, bas_sp) #decent variation, lots of 0's

#now try focal statistics

## sum # road cells within 1sq km (10x10 moving window)
roadsFeature <- raster("data/predictors/roadsFeature.asc") #use 1/0 layer

focalRoad <- raster::focal(roadsFeature, w = matrix(1,9,9), fun = sum)
#check variation among points 
crs(focalRoad) <- crs(roads)

raster::extract(focalRoad, bas_sp) #good variation, probably better than dist2road, but 3 NA sites..
#turn all NA to 0 then water to NA
focalRoad[is.na(focalRoad)] <- 0
focalRoad[is.na(roads)] <- NA

raster::extract(focalRoad, bas_sp) #fixed, keep this one

writeRaster(focalRoad, "data/predictors/focalRoad.asc")



## percent impervious surface
developed <- raster("D:/WoodTurtleSDM/data/Predictors/current_new/imperv_prj.tif")

ext <- drawExtent()

dev.crop <- crop(developed, ext)

dev.prj <- projectRaster(dev.crop, land)

# remove water cells

dev.prj[dev.prj > 100] <- NA

#clip to MDC boundary (this will remove the rest of the ocean cells)
dev.mdc <- mask(x=dev.prj, mask=mdc_land.raster, 
                   filename = "data/predictors/percentImpervious.tif")
writeRaster(dev.mdc, "data/predictors/percentImpervious.asc")



# landcover analysis ----------------------------------------------------------------


# read in 10m, reclassified landcover layer

land.mdc <- raster("D:/lizard_genomics/land_vars/landcover_mdc.tif")

#reclassify
reclass <- read.csv("data/CLC_to_FNAINC_Crosswalk_20160701.csv") %>% dplyr::select(SITE, NAME_SITE, NAME_STATE, reclass)

#land.reclass <- raster::reclassify(land.mdc, as.matrix(reclass %>% dplyr::select(SITE, reclass)))

# read in points

agama_pts <- read.csv("data/agama_sites.csv")

ag_sp <- SpatialPointsDataFrame(agama_pts[,c("LONG", "LAT"),], data = agama_pts)

proj4string(ag_sp) <-  CRS("+proj=longlat +datum=WGS84")

ag_sp <- spTransform(ag_sp, crs(land.mdc))

agama_land <- data.frame(SITE = raster::extract(land.mdc, ag_sp))
agama_land_summary <- left_join(agama_land, reclass, by = "SITE") %>% group_by(NAME_STATE) %>% count(NAME_STATE)
barplot(agama_land_summary$n, names = agama_land_summary$NAME_STATE, las = 2)


# all known localities
ag_all <- read.csv("data/agama_occ_mdc.csv")

ag_all_sp <- SpatialPointsDataFrame(ag_all[,c("Longitude", "Latitude"),], data = ag_all)

proj4string(ag_all_sp) <-  CRS("+proj=longlat +datum=WGS84")

ag_all_sp <- spTransform(ag_all_sp, crs(land.mdc))

agama_all_land <- data.frame(SITE = raster::extract(land.mdc, ag_all_sp))
agama_all_land_sum <- left_join(agama_all_land, reclass, by = "SITE") %>% group_by(NAME_STATE) %>% count(NAME_STATE)


write.csv(agama_land_summary, "data/agama_land_association.csv")
write.csv(agama_all_land_sum, "data/agama_allocc_land_association.csv")


## COMBINE AND REMOVE DUPLICATES
ag_all <- read.csv("data/agama_occ_mdc.csv") %>% select(LONG = Longitude, LAT = Latitude)

agama_pts <- read.csv("data/agama_sites.csv") %>% select(LONG, LAT) 

combined <- rbind(ag_all, agama_pts) %>% distinct(LONG, LAT)

combined_sp <- SpatialPointsDataFrame(combined[,c("LONG", "LAT"),], data = combined)

proj4string(combined_sp) <-  CRS("+proj=longlat +datum=WGS84")

combined_sp <- spTransform(combined_sp, crs(land.mdc))

combined_land <- data.frame(SITE = raster::extract(land.mdc, combined_sp))
combined_land_sum <- left_join(combined_land, reclass, by = "SITE") %>% group_by(NAME_STATE) %>% count(NAME_STATE)

write.csv(combined_land_sum, "data/agama_combined_land_association.csv")

# BASILISKS

bas_pts <- read.csv("data/basilisk_sites.csv")

bas_sp <- SpatialPointsDataFrame(bas_pts[,c("LONG", "LAT"),], data = bas_pts)

proj4string(bas_sp) <-  CRS("+proj=longlat +datum=WGS84")


bas_sp <- spTransform(bas_sp, crs(land.mdc))

bas_land <- data.frame(SITE = raster::extract(land.mdc, bas_sp))
bas_land_sum <- left_join(bas_land, reclass, by = "SITE") %>% group_by(NAME_STATE) %>% count(NAME_STATE)

# canals do not seem to be appropriately recognized in landcover dataset (maybe due to res)

# analyzing dist2canals
raster::extract(dist2canal, bas_sp) %>% .[. < 100]
# 19 found at canals
raster::extract(dist2canal, bas_sp) %>% .[. < 200]
# 34 found within 1 pixel (100m)
raster::extract(dist2canal, bas_sp) %>% .[. < 300]



#COMBINED BASILISK

## COMBINE AND REMOVE DUPLICATES
bv_all <- read.csv("data/basilisk_occ_mdc.csv") %>% select(LONG = Longitude, LAT = Latitude)

bv_pts <- read.csv("data/basilisk_sites.csv") %>% select(LONG, LAT) 

combined <- rbind(bv_all, bv_pts) %>% distinct(LONG, LAT)

combined_sp <- SpatialPointsDataFrame(combined[,c("LONG", "LAT"),], data = combined)

proj4string(combined_sp) <-  CRS("+proj=longlat +datum=WGS84")

combined_sp <- spTransform(combined_sp, crs(land.mdc))

combined_land <- data.frame(SITE = raster::extract(land.mdc, combined_sp))
combined_land_sum <- left_join(combined_land, reclass, by = "SITE") %>% group_by(NAME_STATE) %>% count(NAME_STATE)
#only 26% urban, 26% transportation


write.csv(combined_land_sum, "data/basilisk_combined_land_association.csv")


#extract canopy cover

canopy <- raster("data/predictors/canopy_mdc_final.tif")

bas_pts <- read.csv("data/basilisk_sites.csv")

bas_sp <- SpatialPointsDataFrame(bas_pts[,c("LONG", "LAT"),], data = bas_pts)

proj4string(bas_sp) <-  CRS("+proj=longlat +datum=WGS84")


bas_sp <- spTransform(bas_sp, crs(canopy))

combined_sp <- spTransform(combined_sp, crs(canopy))


bas_canopy <- data.frame(PercentCanopy = raster::extract(canopy, bas_sp))
#average 32% canopy cover, median is 27.4%

#all localities
combined_canopy <- data.frame(PercentCanopy = raster::extract(canopy, combined_sp))
#average 28%, median 23%

#try original 30m canopy cover

canopy <- raster("D:/WoodTurtleSDM/data/Predictors/current_new/canopy_prj.tif")

bas_canopy30 <- raster::extract(canopy, bas_sp)


## contingency analysis ---------------------------------------------

basilisk_land <- read_csv("data/basilisk_combined_land_association.csv") %>% dplyr::select(-X1) %>% 
  rename(land_cover = NAME_STATE, bas_occ = n)

agama_land <- read_csv("data/agama_combined_land_association.csv") %>% dplyr::select(-X1) %>% 
  rename(land_cover = NAME_STATE, agama_occ = n)

land_ass <- full_join(basilisk_land, agama_land, by = "land_cover") %>% 
  mutate(bas_occ = if_else(is.na(bas_occ), 0, bas_occ)) %>% 
  mutate(agama_occ = if_else(is.na(agama_occ), 0, agama_occ)) %>%
  #remove cover types with < 10 occ for both species
  filter(bas_occ > 10 & agama_occ > 10) %>% 
  column_to_rownames(var = "land_cover") 
  

x <- chisq.test(land_ass)
# X-squared = 176.6, df = 5, p-value < 2.2e-16


## confident interval analysis


ci <- read_csv("conf_int.csv") %>% janitor::clean_names()


ap <- ci %>% filter(land_cover_type == "Artificial Pond") %>% 
  ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("red", "blue"))+
  theme(axis.title.x = element_blank())+
  ylab("Point Estimate")+
  labs(title = "Artificial Pond")+
  theme_classic()

cd <- ci %>% filter(land_cover_type == "Canal/Ditch") %>% 
  ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("red", "blue"))+
  theme(axis.title.x = element_blank())+
  ylab("Point Estimate")+
  labs(title = "Canal/Ditch")+
  theme_classic()


crop <- ci %>% filter(land_cover_type == "Cropland/Pasture") %>% 
  ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("red", "blue"))+
  theme(axis.title.x = element_blank())+
  ylab("Point Estimate")+
  labs(title = "Cropland/Pasture")+
  theme_classic()


hu <- ci %>% filter(land_cover_type == "High Intensity Urban") %>% 
  ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("red", "blue"))+
  theme(axis.title.x = element_blank())+
  ylab("Point Estimate")+
  labs(title = "High Intensity Urban")+
  theme_classic()

lu <- ci %>% filter(land_cover_type == "Low Intensity Urban") %>% 
  ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("red", "blue"))+
  theme(axis.title.x = element_blank())+
  ylab("Point Estimate")+
  labs(title = "Low Intensity Urban")+
  theme_classic()

trans <- ci %>% filter(land_cover_type == "Transportation") %>% 
  ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("red", "blue"))+
  theme(axis.title.x = element_blank())+
  ylab("Point Estimate")+
  labs(title = "Transportation")+
  theme_classic()


ggpubr::ggarrange(ap, cd, crop, hu, lu, trans, nrow = 2, ncol = 3)

#try faceting
ci %>% 
ggplot(aes(x = species, y = point_estimate, color = species))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))+
  scale_color_manual(name = "Species", values = c("blue", "red"))+
  ylab("Probability of Occurence")+
  theme_bw()+
  facet_wrap(~ land_cover_type, nrow = 1, strip.position = "bottom")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")



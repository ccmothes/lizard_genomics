#introduction history

library(tidyverse)
library(CoordinateCleaner)
library(ggpubr)

# AGAMA --------------------------------------------------------------
# agama GBIF              
ag_raw <- read_delim("data/gbif_agama/occurrence.txt", delim = "\t")


#select columns of interest
ag <- ag_raw %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)

# remove records without coordinates and dupilcate coordinates
ag <- ag %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude)) %>% 
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>% 
  rename(Longitude = decimalLongitude, Latitude = decimalLatitude) %>% 
  select(year, Longitude, Latitude)

#extract points within mdc

mdc <- readOGR("~/Desktop/NicheProject/MDC_Boundary.shp")

coordinates(ag) <- c("Longitude", "Latitude")
proj4string(ag) <- CRS("+proj=longlat +datum=WGS84")

ag_mdc <- ag[complete.cases(over(ag, mdc)),]
ag_mdc <- ag_mdc[!is.na(ag_mdc$year),]




#agama Eddmaps

library(janitor)

ag_edd <- read.csv("data/eddMaps_agama/mappings_edit.csv")

ag_eddLoc <- ag_edd %>% select(ObsYear, Location, Latitude, Longitude, Datum, Method, CoordAcc) %>% 
  filter(!is.na(ObsYear)) %>% rename(year = ObsYear) %>% 
  select(year, Latitude, Longitude)

ag_eddLoc <- ag_eddLoc %>% remove_empty(which = "rows") %>% drop_na(Latitude, Longitude) %>% filter(Latitude != "NULL") %>% 
  filter(Longitude != "NULL") %>% mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>%
  drop_na(Latitude, Longitude) %>% distinct(Longitude, Latitude, .keep_all = TRUE)

coordinates(ag_eddLoc) <- c("Longitude", "Latitude")
proj4string(ag_eddLoc) <- CRS("+proj=longlat +datum=WGS84")


ag_eddmdc <- ag_eddLoc[complete.cases(over(ag_eddLoc, mdc)),]

ag_all <- rbind(ag_mdc, ag_eddmdc)

write.csv(ag_all, "data/agama_occ_mdc.csv")

#plot points by year

mdc@data$seq_id <- seq(1:nrow(mdc@data))
#fortify
mdc@data$id <- rownames(mdc@data)
#make df
mdcdata <- fortify(mdc, region = "id")
#merge with spatial object
mdcdf <- merge(mdcdata, mdc@data, by = "id")

ag_alldf <- as.data.frame(ag_all)


ggplot()+
  geom_polygon(data = mdcdf, aes(x = long, y = lat, group = group), fill = "lightgrey")+
  theme(panel.background = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
  geom_point(data=as.data.frame(ag_all), aes(Longitude,Latitude), inherit.aes = FALSE, colour = "red")+
  coord_equal()+
  facet_wrap(vars(year))

#plot specific range of years
years <- unique(ag_alldf$year) %>% sort()

p <- vector("list", length = length(years))
for (i in 1:length(years)) {
  p[[i]] <- ggplot() +
    geom_polygon(data = mdcdf,
                 aes(x = long, y = lat, group = group),
                 fill = "lightgrey") +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0,01,0.1,0), "lines")
    ) +
    geom_point(
      data = ag_alldf[ag_alldf$year <= years[i], ],
      aes(Longitude, Latitude),
      inherit.aes = FALSE,
      shape = 15,
      colour = "blue"
    ) +
    coord_equal()
}


#arrange all together
ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]],
          p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]],
          p[[15]], p[[16]], p[[17]], p[[18]], p[[19]], p[[20]],
          nrow = 5, ncol = 4, labels = years[1:20])


# BASILISK ------------------------------------------------------

# agama GBIF              
bv_raw <- read_delim("data/gbif_basilisk/occurrence.txt", delim = "\t")


#select columns of interest
bv <- bv_raw %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)

# remove records without coordinates and dupilcate coordinates
bv <- bv %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude)) %>% 
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>% 
  rename(Longitude = decimalLongitude, Latitude = decimalLatitude) %>% 
  dplyr::select(year, Longitude, Latitude)

#extract points within mdc

mdc <- readOGR("~/Desktop/NicheProject/MDC_Boundary.shp")

coordinates(bv) <- c("Longitude", "Latitude")
proj4string(bv) <- CRS("+proj=longlat +datum=WGS84")

bv_mdc <- bv[complete.cases(over(bv, mdc)),]
bv_mdc <- bv_mdc[!is.na(bv_mdc$year),]

#basilisk Eddmaps

library(janitor)

bv_edd <- read.csv("data/eddMaps_basilisk/mappings_edit.csv")

bv_eddLoc <- bv_edd %>% dplyr::select(Year, Location, Latitude, Longitude, Datum, Method, CoordAcc) %>% 
  filter(!is.na(Year)) %>% rename(year = Year) %>% 
  dplyr::select(year, Latitude, Longitude)

bv_eddLoc <- bv_eddLoc %>% remove_empty(which = "rows") %>% drop_na(Latitude, Longitude) %>% filter(Latitude != "NULL") %>% 
  filter(Longitude != "NULL") %>% mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>%
  drop_na(Latitude, Longitude) %>% distinct(Longitude, Latitude, .keep_all = TRUE)

coordinates(bv_eddLoc) <- c("Longitude", "Latitude")
proj4string(bv_eddLoc) <- CRS("+proj=longlat +datum=WGS84")


bv_eddmdc <- bv_eddLoc[complete.cases(over(bv_eddLoc, mdc)),]

bv_all <- rbind(bv_mdc, bv_eddmdc)

write.csv(bv_all, "data/basilisk_occ_mdc.csv")

#plot points by year

mdc@data$seq_id <- seq(1:nrow(mdc@data))
#fortify
mdc@data$id <- rownames(mdc@data)
#make df
mdcdata <- fortify(mdc, region = "id")
#merge with spatial object
mdcdf <- merge(mdcdata, mdc@data, by = "id")

bv_alldf <- as.data.frame(bv_all)


ggplot()+
  geom_polygon(data = mdcdf, aes(x = long, y = lat, group = group), fill = "lightgrey")+
  theme(panel.background = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
  geom_point(data=as.data.frame(ag_all), aes(Longitude,Latitude), inherit.aes = FALSE, colour = "red")+
  coord_equal()+
  facet_wrap(vars(year))

#plot specific range of years
years <- unique(bv_alldf$year) %>% sort()

p <- vector("list", length = length(years))
for (i in 1:length(years)) {
  p[[i]] <- ggplot() +
    geom_polygon(data = mdcdf,
                 aes(x = long, y = lat, group = group),
                 fill = "lightgrey") +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0,01,0.1,0), "lines")
    ) +
    geom_point(
      data = bv_alldf[bv_alldf$year <= years[i], ],
      aes(Longitude, Latitude),
      inherit.aes = FALSE,
      shape = 17,
      colour = "red"
    ) +
    coord_equal()
}


#arrange all together, removed 1997 because same as 1996 to fit in same dimenstions
# as agama
ggarrange(p[[1]],  p[[3]], p[[4]], p[[5]], p[[6]], p[[7]],
          p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]],
          p[[15]], p[[16]], p[[17]], p[[18]], p[[19]], p[[20]], p[[21]],
          nrow = 5, ncol = 4, labels = years[c(1, 3:21)])


   
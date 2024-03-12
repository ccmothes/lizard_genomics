# SETUP -------------------------------------------------------------------
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("vcfR")
library(tidyverse)
library(janitor)
library(sp)
library(rgdal)
library(raster)
library(conStruct)

#in /data folder

setwd("data/")

basilisk <- structure2conStruct("basilisk_85percent", onerowperind = TRUE,
                             start.loci = 3, start.samples = 3, missing.datum = 0,
                             outfile = "basilisk_conStruct") #samples start at 3 because of 2 headers


# edit rownames (make smaller)
rownames(basilisk) <- substr(rownames(basilisk), 1, 5)

# load coords file

bas_pts <- read.csv("data/basilisk_sites.csv") %>% dplyr::select(LONG, LAT) %>% 
  as.matrix()


#distance matrix

bas_dist <- dist(bas_pts) %>% as.matrix(.)

# run validation run to choose best k value

#set wd to basilisk folder for output files

setwd("data/basilisk/construct_run2/")

choose.k <- x.validation(n.reps = 8, K = 1:5, freqs = basilisk, geoDist = bas_dist,
                               coords = bas_pts, prefix = "bv_choose.k.10fold", save.files = TRUE, 
                               make.figs = TRUE, n.iter = 5000
)

# do again for k 6:10
setwd("data/basilisk/construct_run3/")


choose.k <- x.validation(n.reps = 8, K = 6:10, freqs = basilisk, geoDist = bas_dist,
                         coords = bas_pts, prefix = "bv_choose.k.10fold", save.files = TRUE, 
                         make.figs = TRUE, n.iter = 5000
)


library(raster)
library(conStruct)
library(tidyverse)

# structure to conStruc 

load("data/agama_conStruct.RData")
agama <- freqs #matrix loaded in with name 'freqs'

# edit rownames (make smaller)
rownames(agama) <- substr(rownames(agama), 1, 5)

# load coords file

agama_pts <- read.csv("data/agama_sites.csv") %>% dplyr::select(LONG, LAT) %>% 
  as.matrix()

## need to switch lat and long column order

#distance matrix

agama_dist <- dist(agama_pts) %>% as.matrix(.)

# run validation run to choose best k value

setwd("data/agama/agama_constructRun2")

choose.k.again <- x.validation(n.reps = 10, K = 9:10, freqs = agama, geoDist = agama_dist,
                               coords = agama_pts, prefix = "ag_choose.k.10fold", save.files = TRUE, 
                               make.figs = TRUE, n.iter = 5000)
                               
                               
                               
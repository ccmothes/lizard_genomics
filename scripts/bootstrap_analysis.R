# Bootstrap analysis for Basilisk

library(ResistanceGA)
library(dplyr)


#read in files for each output

## make string of layered names paired woth canopy

layers <- c("canals", "land", "roadsFeature", "traffic", "urban")

## make a list of resistance distance matrices

mat.list <- vector("list", length=length(layers))

for (i in 1:length(layers)){
  
  file <- list.files(path =  paste0("data/basilisk_ms_canopy_", layers[i], "/"),
                     pattern = "jlResistMat", full.names = TRUE)
  
  mat.list[[i]] <- read.csv(file, header = FALSE) %>% as.matrix() 

  
}

# add canopy ss to list

mat.list[[6]] <- read.csv("data/basilisk_ssoptim_update/Results/canopy_jlResistMat.csv",
                          header = FALSE) %>% as.matrix()

layers <- c(layers, "canopy_solo")

## Make list of k (in order of layers)

k <- c(6, 19, 6, 7, 7)

k <- c(k, 4)

## read in genetic distance matrix
bv_dps <- read.csv("data/basilisk_dps.csv")
bv_dps <- bv_dps[,-1]
bv_dps <- bv_dps[-48,-48]

# Run Bootstrap

AIC.boot <- Resist.boot(mod.names = layers, 
                        dist.mat = mat.list,
                        n.parameters = k,
                        iters = 1000,
                        obs = 49,
                        genetic.mat = bv_dps)

# Re-do, forgot to add canopy solo layer

AIC.boot2 <- Resist.boot(mod.names = layers, 
                        dist.mat = mat.list,
                        n.parameters = k,
                        iters = 1000,
                        obs = 49,
                        genetic.mat = bv_dps)




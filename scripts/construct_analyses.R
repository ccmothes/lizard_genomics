# conStruct Analyses

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
library("conStruct")

#load Inauguration palette
#devtools::install_github("ciannabp/inauguration")
library(inauguration)

# color palettes

agama_cols <- c("#e96318", "#2d3437", "#f0b253", "#969a9f")

myidols <- c("#749dae", "#5c1a33", "#f3c483", "#5445b1", "#cd3341","#f7dc6a",
             "#4DA896", "#E07882", "#465952", "#50506D")

# conStruct -----------------------------




# load coords file

agama_pts <- read.csv("data/agama_sites.csv") %>% dplyr::select(LONG, LAT) %>% 
  as.matrix()

ag_id <- read.csv("data/agama_sites.csv")

bas_pts <- read.csv("data/basilisk_sites.csv") %>% dplyr::select(LONG, LAT) %>% 
  as.matrix()

bas_id <- read.csv("data/basilisk_sites.csv")


# * Assess results -----------------------------------

## ** BASILISK -----------------------------------------

# read in results of both runs and combine
sp.results1_5 <- as.matrix(
  read.table("data/basilisk/construct_run2/bv_choose.k.10fold_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

sp.results6_10 <-  as.matrix(
  read.table("data/basilisk/construct_run3/bv_choose.k.10fold_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

sp.results <- rbind(sp.results1_5, sp.results6_10)
rownames(sp.results) <- 1:10

nsp.results1_5 <- as.matrix(
  read.table("data/basilisk/construct_run2/bv_choose.k.10fold_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

nsp.results6_10 <-  as.matrix(
  read.table("data/basilisk/construct_run3/bv_choose.k.10fold_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

nsp.results <- rbind(nsp.results1_5, nsp.results6_10)
rownames(nsp.results) <- 1:10

# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:5 with 8 replicates

plot(rowMeans(sp.results),
     pch=19,col="#cd3341", cex = 2,
     ylab="Predictive Accuracy",xlab="Values of K",
     ylim=range(sp.results,nsp.results))
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "#cd3341",lwd=3)
points(rowMeans(nsp.results),col="#f7dc6a",pch=17, cex = 1.8)
segments(x0 = 1:nrow(nsp.results),
         y0 = nsp.CIs[1,],
         x1 = 1:nrow(nsp.results),
         y1 = nsp.CIs[2,],
         col = "#f7dc6a",lwd=3)

#zoom in on k 4 - 6 spatial
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results[4:6,]),
     main="cross-validation results")


# t-test for sig difference between models
t.test(rowMeans(sp.results), rowMeans(nsp.results), paired = TRUE)
## spatial model sig. better than non spatial, p = 0.003

# sig difference in accuracy between dif values of k?
t.test(sp.results[4,], sp.results[5,], paired = TRUE)
#OoOoOo they are not sig different! close, p = 0.061

# try 3 and 5
t.test(sp.results[3,], sp.results[5,], paired = TRUE)
# 5 is sig better

# 3 and 4
t.test(sp.results[3,], sp.results[4,], paired = TRUE)
# 4 sig better: p = 0.0181


## * * layer contributions -------------------------------

#calculate for each run, then average layer contributions


# SPATIAL MODEL 

layer.contributions <- vector("list", length = 8)

for (i in 1:8) {
  layer.contributions[[i]] <- matrix(NA, nrow = 10, ncol = 10)
  
  #do k1 first
  load(
    paste0(
      "data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep",
      i,
      "K1_conStruct.results.Robj"
    )
  )
  load(
    paste0(
      "data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep",
      i,
      "K1_data.block.Robj"
    )
  )
  
  # calculate layer contributions
  layer.contributions[[i]][, 1] <-
    c(calculate.layer.contribution(conStruct.results[[1]], data.block),
      rep(0, 9))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions
  
  #now add the rest
  for (j in 2:10) {
    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(sprintf(
      paste0(
        "data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep",
        i,
        "K%s_conStruct.results.Robj"
      ),
      j
    ))
    load(sprintf(
      paste0(
        "data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep",
        i,
        "K%s_data.block.Robj"
      ),
      j
    ))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <-
      match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
    
    # calculate layer contributions
    layer.contributions[[i]][, j] <-
      c(
        calculate.layer.contribution(
          conStruct.results = conStruct.results[[1]],
          data.block =
            data.block,
          layer.order =
            tmp.order
        ),
        rep(0, 10 - j)
      )
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order]
  }
}

# average across matrices
layer.contributions.all <- purrr::reduce(layer.contributions, `+`) / length(layer.contributions)

#plot layer contributions spatial
barplot(layer.contributions.all,
        col=myidols,
        xlab="",
        ylab="Layer Contributions",
        names.arg=paste0("K=",1:10))


# NON SPATIAL

layer.contributions.nsp <- vector("list", length = 8)

for (i in 1:8) {
  layer.contributions.nsp[[i]] <- matrix(NA, nrow = 10, ncol = 10)
  
  #do k1 first
  load(
    paste0(
      "data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep",
      i,
      "K1_conStruct.results.Robj"
    )
  )
  load(
    paste0(
      "data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep",
      i,
      "K1_data.block.Robj"
    )
  )
  
  # calculate layer contributions
  layer.contributions.nsp[[i]][, 1] <-
    c(calculate.layer.contribution(conStruct.results[[1]], data.block),
      rep(0, 9))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions
  
  #now add the rest
  for (j in 2:10) {
    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(sprintf(
      paste0(
        "data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep",
        i,
        "K%s_conStruct.results.Robj"
      ),
      j
    ))
    load(sprintf(
      paste0(
        "data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep",
        i,
        "K%s_data.block.Robj"
      ),
      j
    ))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <-
      match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
    
    # calculate layer contributions
    layer.contributions.nsp[[i]][, j] <-
      c(
        calculate.layer.contribution(
          conStruct.results = conStruct.results[[1]],
          data.block =
            data.block,
          layer.order =
            tmp.order
        ),
        rep(0, 10 - j)
      )
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order]
  }
}

# average across matrices
layer.contributions.allnsp <- purrr::reduce(layer.contributions.nsp, `+`) / length(layer.contributions.nsp)

#plot layer contributions spatial
barplot(layer.contributions.allnsp,
        col=myidols,
        xlab="",
        ylab="Nonspatial Layer Contributions",
        names.arg=paste0("K=",1:10))



## * * STRUCTURE PLOTS ----------------------------------------------

#cannot average across runs because layer permutation may differ (e.g. layer 1 in run 1
# may be layer 3 in run 2)

## * * * spatial -----------------------------------------

###K = 2
load("data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep1K2_conStruct.results.Robj")
admix_sp <- conStruct.results$chain_1$MAP$admix.proportions


#order by latitude
make.structure.plot(admix.proportions = admix_sp,
                    sample.order = order(-bas_pts[,2]), layer.colors = myidols,
                    sample.names = bas_id$ID) #order by lat or long


###K = 3
load("data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep1K3_conStruct.results.Robj")
admix_sp3 <- conStruct.results$chain_1$MAP$admix.proportions
# 2= 1, 1 = 2, 3 = 3
admix_sp3 <- admix_sp3[,c(2,1,3)]

#order by latitude
make.structure.plot(admix.proportions = admix_sp3, layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long
###K = 4
load("data/basilisk/construct_run2/bv_choose.k.10fold_sp_rep1K4_conStruct.results.Robj")
admix_sp4 <- conStruct.results$chain_1$MAP$admix.proportions
# 2 = 1, 3 = 2, 1 = 3, 4 = 4
admix_sp4 <- admix_sp4[,c(2, 3, 1, 4)]
#1 = 1, 4 = 2, 2 = 3, 3 = 4
admix_sp4 <- admix_sp4[,c(1,4,2,3)]

#order by latitude
make.structure.plot(admix.proportions = admix_sp4, layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long


## * * * non-spaital ---------------------------------------------------

###K = 2
load("data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep1K2_conStruct.results.Robj")
admix_nsp <- conStruct.results$chain_1$MAP$admix.proportions
#switch to match sp
admix_nsp <- admix_nsp[,c(2,1)]

#order by latitude
make.structure.plot(admix.proportions = admix_nsp, layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long


###K = 3
load("data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep1K3_conStruct.results.Robj")
admix_nsp3 <- conStruct.results$chain_1$MAP$admix.proportions
# 1 = 1, 3 = 2, 2 = 3
admix_nsp3 <- admix_nsp3[,c(1,3,2)]
#switch 1 and 2
admix_nsp3 <- admix_nsp3[,c(2,1,3)]

make.structure.plot(admix.proportions = admix_nsp3, layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]),
                    sample.names = bas_id$ID)


###K = 4
load("data/basilisk/construct_run2/bv_choose.k.10fold_nsp_rep1K4_conStruct.results.Robj")
admix_nsp4 <- conStruct.results$chain_1$MAP$admix.proportions
# 3 = 1, 4 = 2, 1 = 3, 2 = 4
admix_nsp4 <- admix_nsp4[,c(3,4,1,2)]
#switch 1 and 2
admix_nsp4 <- admix_nsp4[,c(2,1,3,4)]

#order by latitude
make.structure.plot(admix.proportions = admix_nsp4, layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long



# * * * structure -----------------------------------------------------------------

# K = 2
struct_k2 <- read.csv("basilisk_gen/basilisk_gen_update/K2_run1_admix.csv", header = FALSE)
# switch 1 and 2
struct_k2 <- struct_k2[,c(2,1)]

#order by latitude
make.structure.plot(admix.proportions = as.matrix(struct_k2), layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long

# K = 3

struct_k3 <- read.csv("basilisk_gen/basilisk_gen_update/K3_run1_admix.csv", header = FALSE)
#rearrange cluster (column) order to visually match K2, 3 = 1; 1 = 2, 2 = 3
struct_k3 <- struct_k3[,c(3,1,2)]
#switch 1 and 2
struct_k3 <- struct_k3[,c(2,1,3)]

#order by latitude
make.structure.plot(admix.proportions = as.matrix(struct_k3), layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long

# K = 4

struct_k4 <- read.csv("basilisk_gen/basilisk_gen_update/K4_run1_admix.csv", header = FALSE)
#rearrange again, 4 = 1, 1 = 2, 2 = 3, 3 = 4
struct_k4 <- struct_k4[,c(4,1,2,3)]
#switch 1 and 2
struct_k4 <- struct_k4[,c(2,1,3,4)]


#order by latitude
make.structure.plot(admix.proportions = as.matrix(struct_k4), layer.colors = myidols,
                    sample.order = order(-bas_pts[,2]), sample.names = bas_id$ID) #order by lat or long


# * * PIE CHART -----------------------------------------------
# add county boundary
library(rgdal)
mdc <- readOGR("~/Desktop/NicheProject/MDC_Boundary.shp")

mdc.crop <- crop(mdc, extent(agama_pts)+0.05)

#remove extra margin space
par(oma=c(0,0,0,0)) 
par(mar=c(0,0,0,0))

plot(mdc.crop, col="#f2f2f2", bg="lightgray", lwd=0.25, border=0 )
# # make background map map
# maps::map(database = "county", xlim = range(bas_pts[,1]) + c(-0.01, 0.01), 
#           ylim = range(bas_pts[,2]) + c(-0.01,0.01), col="gray")


# * * * spatial ----------------------------------------------

# K = 2
make.admix.pie.plot(admix.proportions = admix_sp,
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# K = 3

make.admix.pie.plot(admix.proportions = admix_sp3,
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# K = 4
make.admix.pie.plot(admix.proportions = admix_sp4,
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# * * * non- spatial ----------------------------------------------

# K = 2
make.admix.pie.plot(admix.proportions = admix_nsp,
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# K = 3

make.admix.pie.plot(admix.proportions = admix_nsp3,
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# K = 4
make.admix.pie.plot(admix.proportions = admix_nsp4,
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# * * * structure ----------------------------------------------

# K = 2
make.admix.pie.plot(admix.proportions = as.matrix(struct_k2),
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# K = 3

make.admix.pie.plot(admix.proportions = as.matrix(struct_k3),
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# K = 4
make.admix.pie.plot(admix.proportions = as.matrix(struct_k4),
                    coords = bas_pts,
                    layer.colors = myidols, radii = 2.8,
                    add = TRUE)

# * AGAMA ------------------------------------------------------

# read in results of both runs and combine
sp.results1_5 <- as.matrix(
  read.table("data/agama/choose.k.10fold_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

sp.results6_8 <-  as.matrix(
  read.table("data/agama/ag_choose.k.10fold_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

sp.results9_10 <- as.matrix(
  read.table("data/agama/agama_constructRun2/ag_choose.k.10fold_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

sp.results <- rbind(sp.results1_5, sp.results6_8, sp.results9_10)
rownames(sp.results) <- 1:10

nsp.results1_5 <- as.matrix(
  read.table("data/agama/choose.k.10fold_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

nsp.results6_8 <-  as.matrix(
  read.table("data/agama/ag_choose.k.10fold_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

nsp.results9_10 <-  as.matrix(
  read.table("data/agama/agama_constructRun2/ag_choose.k.10fold_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)


nsp.results <- rbind(nsp.results1_5, nsp.results6_8, nsp.results9_10)
rownames(nsp.results) <- 1:10

# first, get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:5 with 8 replicates

plot(rowMeans(sp.results),
     pch=19,col="#cd3341", cex = 2,
     ylab="Predictive Accuracy",xlab="Values of K",
     ylim=range(sp.results,nsp.results))
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "#cd3341",lwd=3)
points(rowMeans(nsp.results),col="#f7dc6a",pch=17, cex = 1.8)
segments(x0 = 1:nrow(nsp.results),
         y0 = nsp.CIs[1,],
         x1 = 1:nrow(nsp.results),
         y1 = nsp.CIs[2,],
         col = "#f7dc6a",lwd=3)

#zoom in on k 2- 4 spatial
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results[2:4,]),
     xlim=range(2:4),
     main="cross-validation results")


## * * layer contributions -------------------------------

#calculate for each run, then average layer contributions

#SPATIAL    

layer.contributions.ag <- vector("list", length = 10)

for (i in 1:10) {
  layer.contributions.ag[[i]] <- matrix(NA, nrow = 10, ncol = 10)
  
  #do k1 first
  load(
    paste0(
      "data/agama/choose.k.10fold_sp_rep",
      i,
      "K1_conStruct.results.Robj"
    )
  )
  load(
    paste0(
      "data/agama/choose.k.10fold_sp_rep",
      i,
      "K1_data.block.Robj"
    )
  )
  
  # calculate layer contributions
  layer.contributions.ag[[i]][, 1] <-
    c(calculate.layer.contribution(conStruct.results[[1]], data.block),
      rep(0, 9))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions
  
  #now add the rest
  for (j in 2:5) {
    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(sprintf(
      paste0(
        "data/agama/choose.k.10fold_sp_rep",
        i,
        "K%s_conStruct.results.Robj"
      ),
      j
    ))
    load(sprintf(
      paste0(
        "data/agama/choose.k.10fold_sp_rep",
        i,
        "K%s_data.block.Robj"
      ),
      j
    ))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <-
      match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
    
    # calculate layer contributions
    layer.contributions.ag[[i]][, j] <-
      c(
        calculate.layer.contribution(
          conStruct.results = conStruct.results[[1]],
          data.block =
            data.block,
          layer.order =
            tmp.order
        ),
        rep(0, 10 - j)
      )
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order]
  }
  
  for (k in 6:10) {
    load(sprintf(
      paste0(
        "data/agama/ag_choose.k.10fold_sp_rep",
        i,
        "K%s_conStruct.results.Robj"
      ),
      k
    ))
    load(sprintf(
      paste0(
        "data/agama/ag_choose.k.10fold_sp_rep",
        i,
        "K%s_data.block.Robj"
      ),
      k
    ))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <-
      match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
    
    # calculate layer contributions
    layer.contributions.ag[[i]][, k] <-
      c(
        calculate.layer.contribution(
          conStruct.results = conStruct.results[[1]],
          data.block =
            data.block,
          layer.order =
            tmp.order
        ),
        rep(0, 10 - k)
      )
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order] 
    
  }
}

# average across matrices
layer.contributions.agall <- purrr::reduce(layer.contributions.ag, `+`) / length(layer.contributions.ag)

#plot layer contributions
barplot(layer.contributions.agall,
        col=myidols,
        xlab="",
        ylab="Layer Contributions",
        names.arg=paste0("K=",1:10))



#############    NON SPATIAL 

layer.contributions.agnsp <- vector("list", length = 10)

for (i in 1:10) {
  layer.contributions.agnsp[[i]] <- matrix(NA, nrow = 10, ncol = 10)
  
  #do k1 first
  load(
    paste0(
      "data/agama/choose.k.10fold_nsp_rep",
      i,
      "K1_conStruct.results.Robj"
    )
  )
  load(
    paste0(
      "data/agama/choose.k.10fold_nsp_rep",
      i,
      "K1_data.block.Robj"
    )
  )
  
  # calculate layer contributions
  layer.contributions.agnsp[[i]][, 1] <-
    c(calculate.layer.contribution(conStruct.results[[1]], data.block),
      rep(0, 9))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions
  
  #now add the rest
  for (j in 2:5) {
    # load the conStruct.results.Robj and data.block.Robj
    #   files saved at the end of a conStruct run
    load(sprintf(
      paste0(
        "data/agama/choose.k.10fold_nsp_rep",
        i,
        "K%s_conStruct.results.Robj"
      ),
      j
    ))
    load(sprintf(
      paste0(
        "data/agama/choose.k.10fold_nsp_rep",
        i,
        "K%s_data.block.Robj"
      ),
      j
    ))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <-
      match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
    
    # calculate layer contributions
    layer.contributions.agnsp[[i]][, j] <-
      c(
        calculate.layer.contribution(
          conStruct.results = conStruct.results[[1]],
          data.block =
            data.block,
          layer.order =
            tmp.order
        ),
        rep(0, 10 - j)
      )
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order]
  }
  
  for (k in 6:10) {
    load(sprintf(
      paste0(
        "data/agama/ag_choose.k.10fold_nsp_rep",
        i,
        "K%s_conStruct.results.Robj"
      ),
      k
    ))
    load(sprintf(
      paste0(
        "data/agama/ag_choose.k.10fold_nsp_rep",
        i,
        "K%s_data.block.Robj"
      ),
      k
    ))
    
    # match layers up across runs to keep plotting colors consistent
    #   for the same layers in different runs
    tmp.order <-
      match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
    
    # calculate layer contributions
    layer.contributions.agnsp[[i]][, k] <-
      c(
        calculate.layer.contribution(
          conStruct.results = conStruct.results[[1]],
          data.block =
            data.block,
          layer.order =
            tmp.order
        ),
        rep(0, 10 - k)
      )
    tmp <- conStruct.results[[1]]$MAP$admix.proportions[, tmp.order] 
    
  }
}

# average across matrices
layer.contributions.agallnsp <- purrr::reduce(layer.contributions.agnsp, `+`) / length(layer.contributions.agnsp)

#plot layer contributions
barplot(layer.contributions.agallnsp,
        col=myidols,
        xlab="",
        ylab="Layer Contributions",
        names.arg=paste0("K=",1:10))



## * * STRUCTURE PLOTS ----------------------------------------------

#cannot average across runs because layer permutation may differ (e.g. layer 1 in run 1
# may be layer 3 in run 2)

## * * * spatial -----------------------------------------

###K = 2
load("data/agama/choose.k.10fold_sp_rep1K2_conStruct.results.Robj")
admix_sp <- conStruct.results$chain_1$MAP$admix.proportions


#order by latitude
make.structure.plot(admix.proportions = admix_sp,
                    sample.order = order(-agama_pts[,2]), layer.colors = myidols,
                    sample.names = ag_id$ID) #order by lat or long


###K = 3
load("data/agama/choose.k.10fold_sp_rep1K3_conStruct.results.Robj")
admix_sp3 <- conStruct.results$chain_1$MAP$admix.proportions
# 3= 1, 1 = 2, 2 = 3
admix_sp3 <- admix_sp3[,c(3,1,2)]

#order by latitude
make.structure.plot(admix.proportions = admix_sp3, layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long
###K = 4
load("data/agama/choose.k.10fold_sp_rep1K4_conStruct.results.Robj")
admix_sp4 <- conStruct.results$chain_1$MAP$admix.proportions
# 1 = 1, 3 = 2, 2 = 3, 4 = 4
admix_sp4 <- admix_sp4[,c(1, 3, 2, 4)]


#order by latitude
make.structure.plot(admix.proportions = admix_sp4, layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long


## * * * non-spaital ---------------------------------------------------

###K = 2
load("data/agama/choose.k.10fold_nsp_rep1K2_conStruct.results.Robj")
admix_nsp <- conStruct.results$chain_1$MAP$admix.proportions
#switch to match sp
admix_nsp <- admix_nsp[,c(2,1)]

#order by latitude
make.structure.plot(admix.proportions = admix_nsp, layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long


###K = 3
load("data/agama/choose.k.10fold_nsp_rep1K3_conStruct.results.Robj")
admix_nsp3 <- conStruct.results$chain_1$MAP$admix.proportions
# 1 = 1, 3 = 2, 2 = 3
admix_nsp3 <- admix_nsp3[,c(1,3,2)]


make.structure.plot(admix.proportions = admix_nsp3, layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]),
                    sample.names = ag_id$ID)


###K = 4
load("data/agama/choose.k.10fold_nsp_rep1K4_conStruct.results.Robj")
admix_nsp4 <- conStruct.results$chain_1$MAP$admix.proportions
# 1 = 1, 3 = 2, 2 = 3, 4= 4
admix_nsp4 <- admix_nsp4[,c(1,3,2,4)]


#order by latitude
make.structure.plot(admix.proportions = admix_nsp4, layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long



# * * * structure -----------------------------------------------------------------

# K = 2
struct_k2 <- read.csv("data/agama_structurek2.csv", header = FALSE)
# switch 1 and 2
struct_k2 <- struct_k2[,c(2,1)]

#order by latitude
make.structure.plot(admix.proportions = as.matrix(struct_k2), layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long

# K = 3

struct_k3 <- read.csv("data/agama_structurek3.csv", header = FALSE)
#rearrange cluster (column) order to visually match K2, 3 = 1; 1 = 2, 2 = 3
struct_k3 <- struct_k3[,c(2,3,1)]


#order by latitude
make.structure.plot(admix.proportions = as.matrix(struct_k3), layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long

# K = 4

struct_k4 <- read.csv("data/agama_structurek4.csv", header = FALSE)
#rearrange again, 4 = 1, 1 = 2, 2 = 3, 3 = 4
struct_k4 <- struct_k4[,c(4,2,3,1)]



#order by latitude
make.structure.plot(admix.proportions = as.matrix(struct_k4), layer.colors = myidols,
                    sample.order = order(-agama_pts[,2]), sample.names = ag_id$ID) #order by lat or long


# * * PIE CHART -----------------------------------------------
# add county boundary
library(rgdal)
mdc <- readOGR("~/Desktop/NicheProject/MDC_Boundary.shp")

mdc.crop <- crop(mdc, extent(agama_pts)+0.05)

#remove extra margin space
par(oma=c(0,0,0,0)) 
par(mar=c(0,0,0,0))

plot(mdc.crop, col="#f2f2f2", bg="lightgray", lwd=0.25, border=0 )

# # make background map map
# maps::map(database = "county", xlim = range(bas_pts[,1]) + c(-0.01, 0.01), 
#           ylim = range(bas_pts[,2]) + c(-0.01,0.01), col="gray")


# * * * spatial ----------------------------------------------

# K = 2
make.admix.pie.plot(admix.proportions = admix_sp,
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# K = 3

make.admix.pie.plot(admix.proportions = admix_sp3,
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# K = 4
make.admix.pie.plot(admix.proportions = admix_sp4,
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# * * * non- spatial ----------------------------------------------

# K = 2
make.admix.pie.plot(admix.proportions = admix_nsp,
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# K = 3

make.admix.pie.plot(admix.proportions = admix_nsp3,
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# K = 4
make.admix.pie.plot(admix.proportions = admix_nsp4,
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# * * * structure ----------------------------------------------

# K = 2
make.admix.pie.plot(admix.proportions = as.matrix(struct_k2),
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# K = 3

make.admix.pie.plot(admix.proportions = as.matrix(struct_k3),
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

# K = 4
make.admix.pie.plot(admix.proportions = as.matrix(struct_k4),
                    coords = agama_pts,
                    layer.colors = myidols, radii = 3.6,
                    add = TRUE)

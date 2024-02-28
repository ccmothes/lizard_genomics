# Population genomic analyses

# Popgen Metrics --------------------------------------------------

library(hierfstat)
library(adegenet)
library(vcfR)
library(pegas)
library(ade4)
library(raster)
library(tidyverse)

# * Ho and He ------------------------------------------------

agama_vcf <- read.vcfR("D:/lizard_genomics/UNEAK/hapMap/agama_update.vcf")

agama_genind <- vcfR2genind(agama_vcf)
agama_genind@pop <- as.factor(rep(1:nInd(agama_genind))) #assign each to own pop

ag_sum <- summary(agama_genind)
# 8% missing data
# Ho = 0.3126457 +- 0.193SD, Hexp = 0.3797112
# Fis = 1 - ho/he = 
1-0.3126/0.3797
ag_stat <- basic.stats(agama_genind)
# Fis = 0.1885
plot(ag_sum$Hexp, ag_sum$Hobs, main="Observed vs expected heterozygosity")
abline(0,1,col="red")

t.test(ag_sum$Hexp,ag_sum$Hobs,pair=T,var.equal=TRUE,alter="greater")
#Ho significantly smaller than He, p < 2.2e-16, mean diff = 0.06706549

bas_vcf <- read.vcfR("D:/lizard_genomics/UNEAK/hapMap/basilisk_update.vcf")

bas_genind <- vcfR2genind(bas_vcf)
bas_genind@pop <- as.factor(rep(1, nInd(bas_genind)))

bv_sum <- summary(bas_genind)
#9.5% missing data
# Ho = 0.1951562
# He = 0.3143765
1-0.19516/0.31437
bv_stat <- basic.stats(bas_genind)
#Fis = 0.3888

t.test(bv_sum$Hexp,bv_sum$Hobs,pair=T,var.equal=TRUE,alter="greater")
# ho sig lower than He, p < 2.2e-16, mean diff = 0.119



# Kinship -------------------------------------------------------------

BiocManager::install("SeqArray")

library("SNPRelate")
library("SeqArray")

browseVignettes("SNPRelate")

## * AGAMAS ------------------------------------------------------------
#convert to gds file
seqVCF2GDS("D:/lizard_genomics/UNEAK/hapMap/agama_update.vcf", "data/agama.gds")

#extract gds file

genofile <- seqOpen("data/agama.gds")

#LD SNP pruning

set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.4)
## leaving only a few snps....keep original dataset for now, other studies don't
## mention filter out LD snps

# calculate IDB relatedness using plink method of moments MOM

ibd.kin <- snpgdsIBDMoM(genofile,
                    maf=0.05, missing.rate=0.05, num.thread=2, kinship = TRUE)


# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd.kin)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

#save pairwise kinship file
write.csv(ibd.kin$kinship, "data/agama_kinshipPlink.csv")

max(ibd.kin$kinship[lower.tri(ibd.kin$kinship)]) #0.1936073, all relatively unrelated
#this most related pair is ag018 (metro path by campus) and ag046 (homestead library)
quantile(ibd.kin$kinship[lower.tri(ibd.kin$kinship)])
#        0%        25%        50%        75%       100% 
#0.00000000 0.06586653 0.09127849 0.11283983 0.19360726 

#kinship classes
lower(ag_kinship) %>% as.numeric %>% .[. > 0.375]
# 0 first order
lower(ag_kinship) %>% as.numeric %>% .[. > 0.1875]
# 2 pairs second order


## * BASILISKS -----------------------------------------------
#convert to gds file
seqVCF2GDS("D:/lizard_genomics/UNEAK/hapMap/basilisk_update.vcf", "data/basilisk.gds")

#extract gds file

genofile <- seqOpen("data/basilisk.gds")



# calculate IDB relatedness using plink method of moments MOM
## made stricter filter b/c 5% missing was returning uninformative results

ibd.kin.bas <- snpgdsIBDMoM(genofile,
                        maf=0.01, missing.rate=0.01, num.thread=2, kinship = TRUE)

ibd.coeff.bas <- snpgdsIBDSelection(ibd.kin.bas)
head(ibd.coeff.bas)

plot(ibd.coeff.bas$k0, ibd.coeff.bas$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

#save pairwise kinship file
write.csv(ibd.kin.bas$kinship, "data/basilisk_kinshipPlink.csv")

max(ibd.kin.bas$kinship[lower.tri(ibd.kin.bas$kinship)]) #0.3549202
#the most related ind are two from northern miami area
quantile(ibd.kin.bas$kinship[lower.tri(ibd.kin.bas$kinship)])
#0%       25%       50%       75%      100% 
#0.0000000 0.1535936 0.1998862 0.2409288 0.3549202 

#kinship classes
lower(bas_kinship) %>% as.numeric() %>% .[. > 0.375]
# 0 first order
lower(bas_kinship) %>% as.numeric() %>% .[. > 0.1875]
# 731 pairs second order
lower(bas_kinship) %>% as.numeric() %>% .[. < 0.1875]
# 494 distantly related


# IBD --------------------------------------------------------------------------
library(tidyverse)


## * AGAMAS ----------------------------------------------------------

# load coords file

agama_pts <- read.csv("data/agama_sites.csv")

#load site info to add sex
ag_info <- read_csv("collection_sites_updated.csv") %>% dplyr::select(ID, SEX)
agama_pts <- left_join(agama_pts, ag_info, by = "ID")


ag_sp <- SpatialPointsDataFrame(agama_pts[,c("LONG", "LAT"),], data = agama_pts)

proj4string(ag_sp) <-  CRS("+proj=longlat +datum=WGS84")

land <- raster("data/predictors/land_mdc_100m_final.tif")

ag_sp <- spTransform(ag_sp, crs(land)) %>% as.data.frame() %>% 
  dplyr::select(LONG.1, LAT.1) %>% as.matrix()


ag_gendist <- read.csv("data/agama_chord_update.csv")

agama_dps <- propShared(agama_genind)
agama_dps <- read.csv("data/agama_dps.csv")
agama_dps <- agama_dps[,-1] %>% as.dist()


ag_gendist <- ag_gendist[,-1] %>% as.dist(.)

ag_geodist  <- dist(ag_sp) #max dist = ~44km

library(adegenet)

## Mantel test

ibd <- mantel.randtest(ag_gendist, ag_geodist)
# p = 0.761

ibd2 <- mantel.randtest(agama_dps, ag_geodist)
# p = 0.359

ibd_imput <- mantel.randtest(agama_chord_imputed, dist(ag_sp))

library(MASS)
dens <- kde2d(ag_geodist,agama_dps, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(ag_geodist, agama_dps, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(agama_dps~ag_geodist))

# try with euclidean dist
ibd_eucl <- mantel.randtest(dist(agama_genind@tab), ag_geodist)
plot(ibd_eucl, nclass = 30) #p = 0.09 with imputed data


#mantel test in ECODIST
library(ecodist)

mantel(agama_dps ~ ag_geodist, mrank = TRUE)
#  mantelr        pval1        pval2        pval3    llim.2.5%   ulim.97.5% 
#0.009746913  0.405000000  0.596000000  0.829000000 -0.021940004  0.045849736

#Mantel Correlogram

#test data 
Dgen <- dist(spcaIllus$dat2A$tab)
Dgeo <- dist(spcaIllus$dat2A$other$xy)
mantCor <- mgram(Dgen, Dgeo, nperm = 499, nclass = 8)
plot(mantCor)


#Dgen <- dist(agama_genind$tab)
#Dgeo <- dist(ag_sp)

mantCor <- mgram(agama_dps, ag_geodist, nboot = 99999, 
                 stepsize = 5000, 
                 mrank = TRUE)
plot(mantCor)
segments(x0 = mantCor$mgram[,1],
         y0 = mantCor$mgram[,5],
         x1 = mantCor$mgram[,1],
         y1 = mantCor$mgram[,6],
         col = "lightgrey",lwd=3)

#plain dots mean a significant correlation



## * BASILISKS --------------------------------------------------


bas_pts <- read.csv("data/basilisk_sites.csv")

bas_sp <- SpatialPointsDataFrame(bas_pts[,c("LONG", "LAT"),], data = bas_pts)

proj4string(bas_sp) <-  CRS("+proj=longlat +datum=WGS84")

land <- raster("data/predictors/land_mdc_100m_final.tif")

bas_sp <- spTransform(bas_sp, crs(land)) %>% as.data.frame() %>% 
  dplyr::select(LONG.1, LAT.1) %>% as.matrix()


bas_gendist <- read.csv("data/basilisk_chord_update.csv")
bas_gendist <- bas_gendist[,-1] %>% as.dist(.)

bas_dps <- read.csv("data/basilisk_dps.csv")
bas_dps <- bas_dps[,-1] %>% as.dist()


#compare with old chord dist w/ 6,000 SNPs
bas_gendistold <- read.csv("data/basilisk_chord.csv")
bas_gendistold <- bas_gendistold[,-1] %>% as.dist(.)
# no difference, these files are almost perfectly correlated (cor=0.996)

#distance matrix

bas_geodist <- dist(bas_sp) #max dist = 78km

ibd_bas <- mantel.randtest(bas_gendist, bas_geodist)
# p = 0.733
ibd_basold <- mantel.randtest(bas_gendistold, bas_geodist)
# p = 0.823
ibd_bas_new <- mantel.randtest(bas_dps, bas_geodist)
# p = 0.02


dens <- kde2d(bas_geodist,bas_dps, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(bas_geodist, bas_dps, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(bas_dps~bas_geodist))


#mantel in ecodist
mantel(bas_dps ~ bas_geodist, mrank = TRUE)
# p = 0.005

#mantel correlgram
mantCor2 <- mgram(bas_dps, bas_geodist, nboot = 99999, stepsize = 1000,
                  mrank = TRUE)
plot(mantCor2)


#combine agama and basilisk plots
par(mfrow=c(2,1))
plot(mantCor)
plot(mantCor2)

## PCA -----------------------------------------------------------

## * AGAMAS -----------------------------------------------------
agama_vcf <- read.vcfR("data/agama_update.vcf")

agama_genind <- vcfR2genind(agama_vcf)
#fix individual names
rownames(agama_genind@tab) <-  substr(rownames(agama_genind@tab), 1, 5)

#impute missing data
ag_x <- tab(agama_genind, freq = TRUE, NA.method = "mean")


ag.pca <- dudi.pca(ag_x,  scale = TRUE) #keep 40 (almost all) PCs
# plot eigenvalues
barplot(ag.pca$eig, main = "Agama dataset - PCA eigenvalues",
        col = heat.colors(length(ag.pca$eig)))

#Check for outliers

loadingplot(ag.pca$c1^2) #nothing stands out too much, but check TP159371.0
X <- tab(agama_genind)
TP159371.0 <- X[, "TP159371.0"]
table(TP159371.0)
#which individuals are outliers?
rownames(X)[TP159371.0 < 2] # ooo UM agamas 2 of the 4


#
ag_pts <- agama_pts %>% select(x = LONG, y = LAT)
s.value(cbind(1:11,rep(1,11)), -5:5, cleg=0)
text(1:11,rep(1,11), -5:5, col="red",cex=1.5)

s.value(ag_pts, ag.pca$li[,1], add.p=TRUE, cleg=0.5)
#plot

s.label(ag.pca$li, boxes = FALSE, addaxes = TRUE)

colorplot(ag.pca$li, ag.pca$li, transp = TRUE, cex = 3, xlab = "PC1",
          ylab = "PC 2")
abline(v=0,h=0,col="grey", lty=2)
#higher difference in color signifies more differentiation

#try with third pc axix
colorplot(ag.pca$li[c(1,3)], ag.pca$li, transp = TRUE, cex = 3, xlab = "PC1",
          ylab = "PC 3")
abline(v=0,h=0,col="grey", lty=2)

# DAPC -----------------------------------------------------------------------------
#example
data("dapcIllus")
x <- dapcIllus$a
grp <- find.clusters(x, max.n.clust = 40)
dapc1 <- dapc(x, grp$grp)
scatter(dapc1)


## * AGAMAS -----------------------------------------------------
#
agama_vcf <- read.vcfR("data/agama_update.vcf")

agama_genind <- vcfR2genind(agama_vcf)
#fix individual names
rownames(agama_genind@tab) <-  substr(rownames(agama_genind@tab), 1, 5)
#agama_genind@pop <- as.factor(1:nInd(agama_genind)) #assign pop info

#impute missing data
ag_x <- tab(agama_genind, freq = TRUE, NA.method = "mean")

# find clusters
grp2 <- find.clusters(ag_x, max.n.clust = 40)
#keep all PCs
# 1 cluster had lowest BIC, but has to keep >= 2
#keep k = 2 which has a deltaBIC ~2 from K = 1
grp2$grp
#UM agamas are separate group

dapc2 <- dapc(ag_x, grp2$grp, n.da = 100, n.pca = 50)

temp2 <- optim.a.score(dapc2) #optim PC = 1...

dapc2_optima <- dapc(ag_x, grp2$grp, n.da = 100, n.pca = 1)

#kept all PCs and 1 discriminant functions (only 1 available)
scatter(dapc2_optima) 

compoplot(dapc2_optima, posi="bottomright",
          txt.leg=paste("Cluster", 1:2), lab="",
           xlab="individuals")

#try k = 3
grp3 <- find.clusters(ag_x, max.n.clust = 45)
dapc3 <- dapc(ag_x, grp3$grp, n.da = 100, n.pca = 50)
#40PCs and 2 discriminant functions (only 2 available)
scatter(dapc3)

#structure like plot
compoplot(dapc3, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), lab="",
           xlab="individuals")

# a-score
temp3 <- optim.a.score(dapc3) #optim around PC = 6

dapc_optima <- dapc(ag_x, grp3$grp, n.da = 100, n.pca = 6)

scatter(dapc_optima)
compoplot(dapc_optim, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), lab="",
          xlab="individuals")
which(apply(dapc_optim$posterior,1, function(e) all(e<0.9)))

#optim # PCs with cross-validation
xval3 <- xvalDapc(ag_x, grp3$grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval3[2:6] #this says 14 is optimal # PCs
scatter(xval$DAPC)
compoplot(dapc_optim2, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), lab="",
          xlab="individuals")
which(apply(dapc_optim2$posterior,1, function(e) all(e<0.9)))

# k = 4
grp4 <- find.clusters(ag_x, max.n.clust = 40)
dapc4 <- dapc(ag_x, grp4$grp, n.da = 100, n.pca = 50)

scatter(dapc4)

#structure like plot
compoplot(dapc4, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), lab="",
          xlab="individuals")

# a-score
temp4 <- optim.a.score(dapc4) #optim around PC = 8

dapc4_optima <- dapc(ag_x, grp4$grp, n.da = 100, n.pca = 8)

scatter(dapc4_optima)
compoplot(dapc4_optima, posi="bottomright",
          txt.leg=paste("Cluster", 1:4), lab="",
          xlab="individuals")

#optim # PCs with cross-validation
xval4 <- xvalDapc(ag_x, grp4$grp, n.pca.max = 300, training.set = 0.9,
                  result = "groupMean", center = TRUE, scale = FALSE,
                  n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval4[2:6] #this says 8 is optimal # PCs, same as a score
scatter(xval4$DAPC)
compoplot(xval4$DAPC, posi="bottomright",
          txt.leg=paste("Cluster", 1:3), lab="",
          xlab="individuals")




## * BASILISKS --------------------------------------------------------
bas_vcf <- read.vcfR("data/basilisk_update.vcf")

bas_genind <- vcfR2genind(bas_vcf)
bas_genind@pop <- as.factor(rep(1, nInd(bas_genind)))
#impute missing data to use in dapc
x <- tab(bas_genind, freq= TRUE, NA.method = "mean")
bas_genind@tab <- tab(bas_genind, freq= TRUE, NA.method = "mean")
sum(is.na(bas_genind$tab))

grp <- find.clusters(x, max.n.clust = 40)
dapc2 <- dapc(x, grp$grp)
scatter(dapc2)

#basically, DAPC analysis says that both species best described by 1 cluster/group

#try web server (shiny app)
#save genind objects
save(agama_genind, file = "data/agama_genind.RData")
save(bas_genind, file = "data/basilisk_genind.RData")

adegenetServer(what=c("DAPC"))

# sPCA -----------------------------------------------------------------------

# * AGAMAS -------------------------------------------------------------

# use points and genind created above
##note- this function automatically imputes NAs

mySpca <- spca(agama_genind, ag_sp, type=2,scannf=FALSE,plot.nb=FALSE,nfposi=2,nfnega=0)
save(mySpca, file = "data/agama_spca2eig.RData")

#view eigenvectors
barplot(mySpca$eig, main="A variant of the plotnn of sPCA eigenvalues",
        col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
#split eigen values into spatial and variance components
screeplot(mySpca)

#test for presence of global and local structures, use imputed data from before
myGtest <- global.rtest(ag_x,mySpca$lw,nperm=999)
# Monte-Carlo test
# Call: global.rtest(X = ag_x, listw = mySpca$lw, nperm = 999)
# 
# Observation: 0.02465926 
# 
# Based on 999 replicates
# Simulated p-value: 0.072 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 1.562649e+00 2.401281e-02 1.711390e-07 

myLtest <- local.rtest(ag_x,mySpca$lw,nperm=999)
#not sig, p = 0.988



#plot first global structure
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- ag_sp[,1]
y <- ag_sp[,2]
temp <- interp(x, y, mySpca$li[,1])
image(temp, col=azur(100))
points(x,y)

#make smoother
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY)
image(temp, col=azur(100))
points(x,y)

myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
title("sPCA - interpolated map of individual scores")
points(x,y)
}
filled.contour(temp, color.pal=myPal, nlev=50,
               key.title=title("lagged nnscore 1"), plot.title=annot())


# * BASILISKS

bas_spca <- spca()

# MRM EcoDist ------------------------------------------------------------------------

library(ecodist)

# * AGAMAS --------------------------

ag_dps <- read.csv("data/agama_dps.csv")
ag_dps <- ag_dps[,-1]
ag_dps_dist <- as.dist(ag_dps)

geodist <- read.csv("data/agama_ssoptim_update/Results/Distance_jlResistMat.csv", header = FALSE)
canaldist <- read.csv("data/agama_ssoptim_update/Results/canals_jlResistMat.csv", header = FALSE)
canopydist <- read.csv("data/agama_ssoptim_update/Results/canopy_jlResistMat.csv", header = FALSE)
landdist <- read.csv("data/agama_ssoptim_update/Results/land_jlResistMat.csv", header = FALSE)
trafficdist <- read.csv("data/agama_ssoptim_update/Results/roads_jlResistMat.csv", header = FALSE)
impervdist <- read.csv("data/agama_ssoptim_impervious/Results/percentImpervious_jlResistMat.csv", header = FALSE)
roadsdist <- read.csv("data/agama_ssoptim_roadsFeature/Results/roadsFeature_jlResistMat.csv", header = FALSE)


mod <- MRM(ag_dps_dist ~ log(dist(geodist)) + dist(canaldist) + dist(canopydist) +
             dist(landdist) + dist(trafficdist) + dist(impervdist) + dist(roadsdist), 
           mrank = TRUE, nperm = 10000)

#none are significant

# * BASILISKS

#NEW Dps instead of chord
bv_dps <- read.csv("data/basilisk_dps.csv")
bv_dps <- bv_dps[,-1]
bv_dps <- bv_dps[-48,-48]

geodist <- read.csv("data/basilisk_ssoptim_update/Results/Distance_jlResistMat.csv", header = FALSE)
canaldist <- read.csv("data/basilisk_ssoptim_update/Results/canals_jlResistMat.csv", header = FALSE)
canopydist <- read.csv("data/basilisk_ssoptim_update/Results/canopy_jlResistMat.csv", header = FALSE)
landdist <- read.csv("data/basilisk_ssoptim_update/Results/land_jlResistMat.csv", header = FALSE)
trafficdist <- read.csv("data/basilisk_ssoptim_update/Results/roads_jlResistMat.csv", header = FALSE)

mod2 <- MRM(as.dist(bv_dps) ~ dist(geodist) + dist(canaldist) + dist(canopydist) +
              dist(landdist) + dist(trafficdist), 
            mrank = TRUE, nperm = 10000)

##----------------------------------------------------------------------------##
## Author: Christopher P. Catano                                              ## 
## Creation: August 08, 2016                                                  ## 
##                                                                            ## 
## Discription: Analysis of how functional beta diversity changes across      ##
##              productivty levels in the TRC Forest Dynamics Plot, MO        ##
##                                                                            ##
## Input:  Herbaceous plant census wide data frame (HerbaceousCensus2015.csv) ##
##         Trait data                                                         ##
##                                                                            ##
## Output: TBD                                                                ##
##----------------------------------------------------------------------------##


rm(list = ls())

library(reshape2)
library(plyr)
library(vegan)
library(FD)
library(ggplot2)

setwd("C:/Users/Chris/Desktop/Beta-productivity Project/Data/TraitAnalysis_datasets")


#-------------------------------------------------------------------------------
# 1. Import data and prepare for analyses. Species names must be in same order 
#    in both files.
#-------------------------------------------------------------------------------

# Load herbaceous plant census data
species_low <- read.csv("DataSummer_Low.csv", header = TRUE)
species_med <- read.csv("DataSummer_Med.csv", header = TRUE)
species_high <- read.csv("DataSummer_High.csv", header = TRUE)
str(species_low)
str(species_med)
str(species_high)

# Load population mean trait values
traits_low <- read.csv("Species_traits_Low.csv", row.names = 1, header = TRUE)
traits_med <- read.csv("Species_traits_Med.csv", row.names = 1, header = TRUE)
traits_high <- read.csv("Species_traits_High.csv", row.names = 1, header = TRUE)
str(traits_low)
str(traits_med)
str(traits_high)

newlow <- colSums(species_low)
newlow2 <- rbind(newlow, species_low)
newlow3 <- newlow2[-c(2:31), ]

newmed <- colSums(species_med)
newmed2 <- rbind(newmed, species_med)
newmed3 <- newmed2[-c(2:31), ]

newhigh <- colSums(species_high)
newhigh2 <- rbind(newhigh, species_high)
newhigh3 <- newhigh2[-c(2:31), ]

# Check to make sure same number of species in both files (CRITICAL!)
if(dim(species_low)[2]!=dim(traits_low)[1])
  stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")
if(dim(species_med)[2]!=dim(traits_med)[1])
  stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")
if(dim(species_high)[2]!=dim(traits_high)[1])
  stop("error:differentnumberofspeciesin'traits'and'abundances'matrices")

# Remove ITV columns and traits that include petiole
traits_low_nopet <- traits_low[, -c(2:4,6:8,10:12,14,16)]
traits_med_nopet <- traits_med[, -c(2:4,6:8,10:12,14,16)]
traits_high_nopet <- traits_high[, -c(2:4,6:8,10:12,14,16)]

# Remove ITV columns and keep traits that include petiole
traits_low_pet <- traits_low[, -c(1:2,4:6,8:10,12,14,16)]
traits_med_pet <- traits_med[, -c(1:2,4:6,8:10,12,14,16)]
traits_high_pet <- traits_high[, -c(1:2,4:6,8:10,12,14,16)]

# Scale traits to unit variance 
traits1 <- apply(traits_low_nopet, 2, scale)
traits2 <- apply(traits_med_nopet, 2, scale)
traits3 <- apply(traits_high_nopet, 2, scale)


#-------------------------------------------------------------------------------
# 2. Run functional diversity code to calculate community weighted means and
#    other metrics.
#-------------------------------------------------------------------------------

# Functional diversity with un-scaled traits
low_FD_nopet <- dbFD(traits_low_nopet, species_low, corr = c("lingoes"),w.abun=TRUE,stand.x=TRUE)
med_FD_nopet <- dbFD(traits_med_nopet, species_med, corr = c("lingoes"),w.abun=TRUE,stand.x=TRUE)
high_FD_nopet <- dbFD(traits_high_nopet, species_high, corr = c("lingoes"), w.abun=TRUE,stand.x=TRUE)


# Functional diversity with scaled traits so units are comparable 
traits1_low <- read.csv("FD_low_scaled.csv", row.names = 1, header = TRUE)
str(traits1_low)
traits2_med <- read.csv("FD_med_scaled.csv", row.names = 1, header = TRUE)
str(traits2_med)
traits3_high <- read.csv("FD_high_scaled.csv", row.names = 1, header = TRUE)
str(traits3_high)

low_FD_scaled <- dbFD(traits1_low, species_low, corr = c("lingoes"),w.abun=TRUE,stand.x=TRUE)
med_FD_scaled <- dbFD(traits2_med, species_med, corr = c("lingoes"),w.abun=TRUE,stand.x=TRUE)
high_FD_scaled <- dbFD(traits3_high, species_high, corr = c("lingoes"), w.abun=TRUE,stand.x=TRUE)

#write.csv(low_FD_scaled, file = "low_FD_scaled1.csv")
#write.csv(med_FD_scaled, file = "med_FD_scaled1.csv")
#write.csv(high_FD_scaled, file = "high_FD_scaled1.csv")


# aaa <- calc.FD(species_low, traits_low_nopet)


# plot species trait composition in multivariate space  ------------------------

traits1<-apply(traits_low_nopet,2,scale)
traits2<-apply(traits_med_nopet,2,scale)
traits3<-apply(traits_high_nopet,2,scale)

D1 <- dist(traits1)
mds1 <- metaMDS(D1, k=2)
mds.points1 <- mds1$points
low <- as.matrix(colSums(species_low))
newlow <- cbind(mds.points1,low)
newlow[,3] <- newlow[,3]/sum(newlow[,3])

D2 <- dist(traits2)
mds2 <- metaMDS(D2, k=2)
mds.points2 <- mds2$points
med <- as.matrix(colSums(species_med))
newmed <- cbind(mds.points2,med)
newmed[,3] <- newmed[,3]/sum(newmed[,3])

D3 <- dist(traits3)
mds3 <- metaMDS(D3, k=2)
mds.points3 <- mds3$points
high <- as.matrix(colSums(species_high))
newhigh <- cbind(mds.points3,high)
newhigh[,3] <- newhigh[,3]/sum(newhigh[,3])

par(mfrow = c(1,3))
plot(mds.points1[,2] ~ mds.points1[,1], pch=20, col = "black", cex=(20*newlow[,3]),
     ylim = c(-3,6), xlim = c(-3, 5)) 

plot(mds.points2[,2]~mds.points2[,1], pch=20, col = "red", cex=(20*newmed[,3]),
     ylim = c(-3,6), xlim = c(-3, 5))

plot(mds.points3[,2]~mds.points3[,1], pch=20, col = "blue", cex=(20*newhigh[,3]),
     ylim = c(-3,6), xlim = c(-3, 5))

#par(mfrow = c(1,1))
#plot(mds.points1[,2]~mds.points1[,1], pch=20, col = "black", cex=(20*newlow[,3]),
#ylim = c(-3,4), xlim = c(-4, 5))
#points(mds.points2[,2]~mds.points2[,1], pch=19, col = "red", cex=(20*newmed[,3]))
#points(mds.points3[,2]~mds.points3[,1], pch=18, col = "blue", cex=(20*newhigh[,3]))

#---------------------------------------------------------------------------------

##write table with FD metrics
#write.csv(low_FD_nopet, file = "FD_low_nopet_st.csv")  
#write.csv(med_FD_nopet, file = "FD_med_nopet_st.csv")  
#write.csv(high_FD_nopet, file = "FD_high_nopet_st.csv")  

#write.csv(traits1, file = "FD_low_scaled.csv")  
#write.csv(traits2, file = "FD_med_scaled.csv")  
#write.csv(traits3, file = "FD_high_scaled.csv")  

#write.csv(traits_low_nopet, file = "low_names.csv")  
#write.csv(traits_med_nopet, file = "med_names.csv")  
#write.csv(traits_high_nopet, file = "high_names.csv")  

# for spfd
low_SPFD <- dbFD(traits1_low, newlow3, corr = c("lingoes"))
med_SPFD <- dbFD(traits2_med, newmed3, corr = c("lingoes"))
high_SPFD <- dbFD(traits3_high, newhigh3, corr = c("lingoes"))

##write table with FD metrics
#write.csv(low_SPFD, file = "SPFD_low_scaled.csv")  
#write.csv(med_SPFD, file = "SPFD_med_scaled.csv")  
#write.csv(high_SPFD, file = "SPFD_high_scaled.csv")  

FunDiv <- read.csv("FunctionalDiversity_nopet.csv", header = TRUE)
FunDiv_scaled <- read.csv("FunctionalDiversity_Scaled_nopet.csv", header = TRUE)

# functional beta diversity based on trait composition using community weighted
# means of all traits
DIS <- vegdist(FunDiv[, 10:15], method = "euclidean")
DIS_scaled <- vegdist(FunDiv_scaled[, 12:17], method = "euclidean")

mod_CWM <- betadisper(DIS, FunDiv$productivity, type = "median", 
                      bias.adjust=FALSE)
mod_CWM_scaled <- betadisper(DIS_scaled, FunDiv_scaled$productivity, type = "median", 
                             bias.adjust=FALSE)

# get different results using scaled or unscaled functional diversity metrics.
# unscaled shows functional beta signficantly lower than either high or medium
# productivity habitats.
# scaled shows functional beta significantly higher in medium compared to high, 
# but not low productivity habitats. 

par(mfrow = c(1,1))
plot(mod_CWM, axes = c(1,2), hull = FALSE)
boxplot(mod_CWM)
permutest(mod_CWM, pairwise = TRUE, permutations = 9999)

mod_CWM$distances <- log(mod_CWM$distances) # transform distances, highly right skewed
plot(mod_CWM, axes = c(1,2), hull = FALSE)
boxplot(mod_CWM)
permutest(mod_CWM, pairwise = TRUE, permutations = 9999)
# P = 0.0536 (low:med = 0.0625, low:high = 0.0283) (low has sig less beta)

plot(mod_CWM_scaled, axes = c(1,2), hull = FALSE)
boxplot(mod_CWM_scaled)
permutest(mod_CWM_scaled, pairwise = TRUE, permutations = 9999)

mod_CWM_scaled$distances <- log(mod_CWM_scaled$distances) # transform distances, highly right skewed
plot(mod_CWM_scaled, axes = c(1,2), hull = FALSE)
boxplot(mod_CWM_scaled)
permutest(mod_CWM_scaled, pairwise = TRUE, permutations = 9999)
# P = 0.0683 (high:med = 0.015) (med has sig high beta)



# functional beta diversity based on trait composition using community weighted
# means of LA, SLA, Hmax, and seed mass (Lth and LDMC highly correlated with SLA)
# RESULT IS SAME AS ABOVE SO USE ALL TRAITS
cor(FunDiv[, 10:15])

DIS1 <- vegdist(FunDiv[, c(10:11,14:15)], method = "euclidean")
mod_CWM1 <- betadisper(DIS1, FunDiv$productivity, type = "median", 
                       bias.adjust=FALSE)  

mod_CWM1$distances <- log(mod_CWM1$distances) # transform distances, highly right skewed
boxplot(mod_CWM1) 
permutest(mod_CWM1, pairwise = TRUE, permutations = 9999) 
# P = 0.0561 low is lowest, same as first group of traits


DIS1_scaled <- vegdist(FunDiv_scaled[, c(12:13,16:17)], method = "euclidean")
mod_CWM1_scaled <- betadisper(DIS1_scaled, FunDiv_scaled$productivity, type = "median", 
                              bias.adjust=FALSE)  

mod_CWM1_scaled$distances <- log(mod_CWM1_scaled$distances) # transform distances, highly right skewed
boxplot(mod_CWM1_scaled) 
permutest(mod_CWM1_scaled, pairwise = TRUE, permutations = 9999) 
# P = 0.0353 (med sig higher than low or high)




# functional beta diversity based on trait composition using community weighted
# means of LA, Hmax, and seed mass 

DIS3 <- vegdist(FunDiv[, c(10,14:15)], method = "euclidean")
mod_CWM3 <- betadisper(DIS3, FunDiv$productivity, type = "median", 
                       bias.adjust=FALSE)  

mod_CWM3$distances <- log(mod_CWM3$distances)
boxplot(mod_CWM3) 
permutest(mod_CWM3, pairwise = TRUE, permutations = 9999) 
# P = 0.0.0016 (pairwise permutations signficant, med:low = 0.0118, med:high = 0.0005 )

DIS3_scaled <- vegdist(FunDiv_scaled[, c(12,16:17)], method = "euclidean")
mod_CWM3_scaled <- betadisper(DIS3_scaled, FunDiv_scaled$productivity, type = "median", 
                              bias.adjust=FALSE)  

mod_CWM3_scaled$distances <- log(mod_CWM3_scaled$distances)
boxplot(mod_CWM3_scaled) 
permutest(mod_CWM3_scaled, pairwise = TRUE, permutations = 9999) 
# P = 0.0.0025 (pairwise permutations signficant, med:low = 0.0096, med:high = 0.0004 )






par(mfrow = c(1, 2))
boxplot(mod_CWM, main = "Functional beta-diversity (LA, SLA, LDMC, Lth, Hmax, Seed mass)")
boxplot(mod_CWM3, main = "Functional beta-diversity (LA, Hmax, Seed mass)") 

#Calculate functional beta diversity as CWM of individual traits
LA <- FunDiv[ ,c(10)]
SLA <- FunDiv[ ,c(11)]
LDMC <- FunDiv[ ,c(12)]
Lth <- FunDiv[ ,c(13)]
Hmax <- FunDiv[ ,c(14)]
seedMass <- FunDiv[ ,c(15)]
Fdis <- FunDiv[ ,c(8)]

Beta_LA <- vegdist(LA, method = "euclidean")
Beta_SLA <- vegdist(SLA, method = "euclidean")
Beta_LDMC <- vegdist(LDMC, method = "euclidean")
Beta_Lth <- vegdist(Lth, method = "euclidean")
Beta_Hmax <- vegdist(Hmax, method = "euclidean")
Beta_seedMass <- vegdist(seedMass, method = "euclidean")
Beta_LA <- vegdist(LA, method = "euclidean")
Beta_Fdis <- vegdist(Fdis, method = "euclidean")

mod_LA <- betadisper(Beta_LA, FunDiv$productivity, type = "median", bias.adjust=FALSE)
mod_SLA <- betadisper(Beta_SLA, FunDiv$productivity, type = "median", bias.adjust=FALSE)
mod_LDMC <- betadisper(Beta_LDMC, FunDiv$productivity, type = "median", bias.adjust=FALSE)
mod_Lth <- betadisper(Beta_Lth, FunDiv$productivity, type = "median", bias.adjust=FALSE)
mod_Hmax <- betadisper(Beta_Hmax, FunDiv$productivity, type = "median", bias.adjust=FALSE)
mod_seedMass <- betadisper(Beta_seedMass, FunDiv$productivity, type = "median", bias.adjust=FALSE)
mod_Fdis <- betadisper(Beta_Fdis, FunDiv$productivity, type = "median", bias.adjust=FALSE)

#Leaf area
mod_LA$distances <- log(mod_LA$distances+1)
boxplot(mod_LA) 
permutest(mod_LA, pairwise = TRUE, permutations = 9999)
# P <0.001 (med sign dif than low or high, high no diff than low)

#Specific leaf area
mod_SLA$distances <- log(mod_SLA$distances+1)
boxplot(mod_SLA) 
permutest(mod_SLA, pairwise = TRUE, permutations = 9999)
# P 0.0411 (low:high = 0.0113, med:high = 0.054)

#leaf dry matter content
mod_LDMC$distances <- log(mod_LDMC$distances)
boxplot(mod_LDMC) 
permutest(mod_LDMC, pairwise = TRUE, permutations = 9999)
# P 0.7262 

# Leaf thickness
#mod_Lth$distances <- log(mod_Lth$distances)
boxplot(mod_Lth) 
permutest(mod_Lth, pairwise = TRUE, permutations = 9999)
#P < 0.0001 (low:high = P < 0.001, low:med = P < 0.001)

# maximum plant height
boxplot(mod_Hmax) 
permutest(mod_Hmax, pairwise = TRUE, permutations = 9999)
# P < 0.0001 (low:high = P = 0.0057, low:med = P = 0.0002)

# seed mass
mod_seedMass$distances <- log(mod_seedMass$distances)
boxplot(mod_seedMass)
permutest(mod_seedMass, pairwise = TRUE, permutations = 9999)
# P = 0.0011 (med:low P = 0.0005, med:high P = 0.031, low:high P = 0.0973)

# functional dispersion
mod_Fdis$distances <- log(mod_Fdis$distances)
boxplot(mod_Fdis)
permutest(mod_Fdis, pairwise = TRUE, permutations = 9999)
# P = 0.0011 (med:low P = 0.0005, med:high P = 0.031, low:high P = 0.0973)


par(mfrow = c(2,3))
boxplot(mod_LA, main = "LA", ylab = "ln(Distance to centroid + 1)")
boxplot(mod_SLA, main = "SLA", ylab = "ln(Distance to centroid + 1)")
boxplot(mod_LDMC, main = "LDMC", ylab = "ln(Distance to centroid)")
boxplot(mod_Lth, main = "Lth")
boxplot(mod_Hmax, main = "Hmax")
boxplot(mod_seedMass, main = "seed mass", "ylab = ln(Distance to centroid)")


















































FunDiv_low <- read.csv("FD_low.csv", header = TRUE)
FunDiv_med <- read.csv("FD_med.csv", header = TRUE)
FunDiv_high <- read.csv("FD_high.csv", header = TRUE)


low_DIS = vegdist(FunDiv_low[,10:15], method="euclidean")
med_DIS = vegdist(FunDiv_med[,10:15], method="euclidean")
high_DIS = vegdist(FunDiv_high[,10:15], method="euclidean")

betadisp.Low <- betadisper(low_DIS, rep(1, length(FunDiv_low$X)),
                           type = "median", bias.adjust=FALSE) 
hist(betadisp.Low$distances)
betadisp.Med <- betadisper(med_DIS, rep(1, length(FunDiv_med$X)),
                           type = "median", bias.adjust=FALSE)
betadisp.High <- betadisper(high_DIS, rep(1, length(FunDiv_high$X)),
                            type = "median", bias.adjust=FALSE)

boxplot(low_DIS, med_DIS, high_DIS,  	
        ylab = "Dissimilarity", names = c("Low", "Medium", "High"), 
        xlab = "Productivity level", col = c("grey90", "grey90", "grey90"))

boxplot(betadisp.Low$distance, betadisp.Med$distance, 
        betadisp.High$distance, main = "Functional beta-diversity", 
        ylab = "Distance-to-median",  
        xlab = "Productivity", names = c("Low", "Medium", "High"),
        col = c("grey90", "grey90", "grey90")) 

mean(low_DIS)
seLOW=sd(low_DIS)/sqrt(640)

TY_DIS = vegdist(TY[,2:7], method="euclidean")
OBSB=mean(TY_DIS)
seTY=sd(TY_DIS)/sqrt(506)

SC_DIS = vegdist(SC, method="euclidean")
mean(SC_DIS)
seTY=sd(SC_DIS)/sqrt(379)

par(mfrow = c(2,3))
boxplot(FunDiv$CWM.LA_MEAN ~ FunDiv$productivity, main = "CWM LA")
boxplot(FunDiv$CWM.SLA_MEAN ~ FunDiv$productivity, main = "CWM SLA")
boxplot(FunDiv$CWM.LDMC_MEAN ~ FunDiv$productivity, main = "CWM LDMC")
boxplot(FunDiv$CWM.Lth_MEAN ~ FunDiv$productivity, main = "CWM Lth")
boxplot(FunDiv$CWM.Hmax_MEAN ~ FunDiv$productivity, main = "CWM Hmax")
boxplot(FunDiv$CWM.SeedMass ~ FunDiv$productivity, main = "CWM seed mass")

par(mfrow = c(2,2))
boxplot(FunDiv$FEve ~ FunDiv$productivity, main = "Functional Evenness (m-2 plot)")

boxplot(FunDiv$FDiv ~ FunDiv$productivity, main = "Functional Divergence (m-2 plot)")

boxplot(log(FunDiv$FRic) ~ FunDiv$productivity, main = "Functional Richness (m-2 plot)") 
# volume of niche space occupied 
# by the species in the communities (Mouchet et al. 2010). FRic should increase
# when niche complementarity enhances probabilities of species occurence.
boxplot(log(FunDiv$FDis+1) ~ FunDiv$productivity, main = "Functional Dispersion (m-2 plot)") 
# measures degree to which the 
# abundance of a community is distributed toward the extremities of occupied 
# trait space. Should increase when niche complementarity enhances species'
# relative abundances. SPecifically FDisp should increase when niche complementarity 
# enhances either, or both, species' occurrence probabilities and abundances

















low_FD_pet <- dbFD(traits_low_pet, species_low, corr = c("lingoes"))
med_FD_pet <- dbFD(traits_med_pet, species_med, corr = c("lingoes"))
high_FD_pet <- dbFD(traits_high_pet, species_high, corr = c("lingoes"))

write.csv(low_FD_pet, file = "FD_low_pet.csv")  
write.csv(med_FD_pet, file = "FD_med_pet.csv")  
write.csv(high_FD_pet, file = "FD_high_pet.csv")  

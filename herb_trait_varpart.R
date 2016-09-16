##----------------------------------------------------------------------------##
## Author: Christopher P. Catano                                              ## 
## Creation: September 7, 2016                                                ## 
##                                                                            ## 
## Discription: Create spatial variables and run variance partitioning to     ##  
##              patterns of beta diversity (taxanomic and function) across    ##
##              environmental productivity levels                             ##
##                                                                            ##
## Input:  Herbaceous plant census wide data frame (HerbaceousCensus2015.csv) ##
##                                                                            ##
##                                                                            ##
## Output: TBD                                                                ##
##----------------------------------------------------------------------------##


rm(list = ls())

library(reshape2)
library(plyr)
library(vegan)
library(ggplot2)
library(packfor)
library(spacemakeR)
library(ade4)
library(spdep)
library(adespatial)


par(mfrow=c(2,3))

#-------------------------------------------------------------------------------
# 1. Import data and prepare for analyses
#-------------------------------------------------------------------------------

# Load plant community data, convert plot to character, and order by plot.
setwd("C:/Users/Chris/Desktop/Beta-productivity Project/Data/TraitAnalysis_datasets")
traits <- read.csv("FunctionalDiversity_Scaled_nopet.csv", header = TRUE)
traits$plot <- as.character(traits$plot)
traits <- traits[order(traits$plot),]

# Load environmental data, convert plot to character, and order by plot.
setwd("C:/Users/Chris/Desktop/Beta-productivity Project/Data/ReadyForR")
data <- read.csv("Environmental_Data_All_08102016.csv")
data$plot <- as.character(data$plot)
data <- data[order(data$plot),]

# extract data for specific to each soil productivity level
traits.low <- traits[traits$productivity == "low", ]
traits.med <- traits[traits$productivity == "med", ]
traits.high <- traits[traits$productivity == "high", ]

# extract data for specific to each soil productivity level
data.low <- data[data$productivity == "low", ]
data.med <- data[data$productivity == "med", ]
data.high <- data[data$productivity == "high", ]

# Create coordinate matrices for each productivity level
low.xy <- data.low[, c(6,7)]
med.xy <- data.med[, c(6,7)]
high.xy <- data.high[, c(6,7)]

plot(low.xy, col="red", pch = 15)
plot(med.xy, col="orange", pch = 15)
plot(high.xy, col="green", pch = 15)

# Low Productivity Site --------------------------------------------------------
# pca on all abiotic environmental variables
env.low <- data.low[, -c(4:12,32,36:37,39,40,43,46)] 
low.pca <- prcomp(env.low[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(low.pca, type = "l")
summary(low.pca)
low.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(low.pca, scale = 0)
pc.scores.low <- as.data.frame(low.pca$x[ ,1:7])
pc.min <- abs(min(pc.scores.low[ ,1:6])) 
env.het.data.low <- pc.scores.low + pc.min
env.het.low <- vegdist(env.het.data.low[ ,1:6], method = "euclidean")
mod_env.het.low <- betadisper(env.het.low, rep(1, length(data.low$plot)),
                              type = "median", bias.adjust=FALSE)

# calculate fractions of variance explained using distance-based redundancy 
# analysis (dbRDA). dbRDA takes distance matrix, then makes an ordination 
# (i.e., principal coordinates analysis), then a constrained ordination in which 
# ordination axes are a linear function of environmental and spatial variables.
# function "capscale" (using a formula with only a constant on the left hand side) 
# takes the distance and creates a multivariate space that preserves those distances 
# (i.e., performs principal coordinates analysis); the main message: the result 
# of principal coordinate analysis (the object that below is named "traits.model") 
# is the response variable used for variance partitioning using function 'varpart'.
traits.model.low <- capscale(traits.low[, 12:17] ~ 1)
traits.low.sc = scores(traits.model.low, display = c("wa"), scaling =1)
traits.eu.low <- vegdist(traits.low[, 12:17], method = "euclidean")


# There are multiple ways to model the spatial component. Linear distance of 
# coordinates and Moran's eigenvector maps (PCNM= dbMEM). Because the study 
# was designed and sampled in a way to minimize and quilize space-env. correlations
# across the three producitivity levels, linear distance is most appropriate (there
# is no reason to expect complex spatial patterns at the scales we are addressing).
# That being said, the code below using both linear distance and dbMEM (PCNM) with
# forward model selection just to check. As expected, both approaches give similar
# almost identical results. Also, after linear detrending the data, no dbMEM axes
# are signficant. These results confirm the simpler linera distance is most appropriate.

# 1. distance-based Moran's eigenvector maps (equal to PCNM)  
#threshold for the minimum spanning tree set to 87 (max distance among plots
#within a local site. Setting this value to null fails)
low.dbmem1 <- dbmem(low.xy, thresh = 87.290463, MEM.autocor = "positive")
# Print the (n-1) non-null eigenvalues
attributes(low.dbmem1)$values # 3 positive eigenvalues
# Plot the associated spatial weighting matrix
s.label(low.xy)
# Plot maps of the first 3 dbMEM eigenfunctions
s.value(low.xy, low.dbmem1[,1:3])
# Compute and test associated Moran's I values
# Eigenvalues are proportional to Moran's I
test <- moran.randtest(low.dbmem1, nrepet = 99)
plot(test$obs, attr(low.dbmem1, "values"), xlab = "Moran's I", ylab = "Eigenvalues")
# Decreasing values of Moran's I for the successive MEM.
# The red line is the expected value of Moran's I under H0.
plot(test$obs, xlab="MEM rank", ylab="Moran's I")
abline(h=-1/(nrow(low.xy) - 1), col="red")
#convert dbMEM result to matrix for RDA
low.dbMEM <- as.matrix(low.dbmem1)


# 2. Run the global dbMEM analysis on the *detrended* data
# Is there a linear trend? If yes, detrend the data before assessing more
# complicated spatial structures.
anova(rda(traits.low.sc, low.xy), step=1000)	# Result: no significant trend
(low.PCNM.rda <- rda(traits.low.sc, low.dbMEM))
anova(low.PCNM.rda, step=1000) # global model not significant

# this code will give you the adjusted r2 for a+b in the variance partitioning
cap.mod <-capscale(traits.eu.low ~ PC1+PC2+PC3+PC4+PC5+PC6, pc.scores.low,
                   dist="euclidean", add =TRUE)
anova(cap.mod, by = "margin")
RsquareAdj(cap.mod)$adj.r.squared
plot(cap.mod)

# comute RDA and variance partitioning using all 7 environmetal PCs and linear coordinates 
#low.RDA <- rda(traits.low[,12:17], pc.scores.low[,1:6], low.xy)
#anova(low.RDA)
low.varpart.linearxy <- varpart(traits.low.sc, ~ PC1+PC2+PC3+PC4+PC5+PC6, low.xy, data = pc.scores.low)
low.varpart.linearxy
low.varpart.linearxy.eu <- varpart(traits.eu.low, ~ PC1+PC2+PC3+PC4+PC5+PC6, low.xy, data = pc.scores.low)
low.varpart.linearxy.eu



# Med Productivity Site --------------------------------------------------------
# pca on all abiotic environmental variables
env.med <- data.med[, -c(4:12,32,36:37,39,40,43,46)] 
med.pca <- prcomp(env.med[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(med.pca, type = "l")
summary(med.pca)
med.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(med.pca, scale = 0)
pc.scores.med <- as.data.frame(med.pca$x[ ,1:7])
pc.min.med <- abs(min(pc.scores.med[ ,1:6])) 
env.het.data.med <- pc.scores.med + pc.min.med
env.het.med <- vegdist(env.het.data.med[ ,1:6], method = "euclidean")
mod_env.het.med <- betadisper(env.het.med, rep(1, length(data.med$plot)),
                              type = "median", bias.adjust=FALSE)

# calculate fractions of variance explained using distance-based redundancy 
# analysis (dbRDA). dbRDA takes distance matrix, then makes an ordination 
# (i.e., principal coordinates analysis), then a constrained ordination in which 
# ordination axes are a linear function of environmental and spatial variables.
# function "capscale" (using a formula with only a constant on the left hand side) 
# takes the distance and creates a multivariate space that preserves those distances 
# (i.e., performs principal coordinates analysis); the main message: the result 
# of principal coordinate analysis (the object that below is named "traits.model") 
# is the response variable used for variance partitioning using function 'varpart'.
traits.model.med <- capscale(traits.med[, 12:17] ~ 1)
traits.med.sc = scores(traits.model.med, display = c("wa"), scaling =1)
traits.eu.med <- vegdist(traits.med[, 12:17], method = "euclidean")

# There are multiple ways to model the spatial component. Linear distance of 
# coordinates and Moran's eigenvector maps (PCNM= dbMEM). Because the study 
# was designed and sampled in a way to minimize and quilize space-env. correlations
# across the three producitivity levels, linear distance is most appropriate (there
# is no reason to expect complex spatial patterns at the scales we are addressing).
# That being said, the code bemed using both linear distance and dbMEM (PCNM) with
# forward model selection just to check. As expected, both approaches give similar
# almost identical results. Also, after linear detrending the data, no dbMEM axes
# are signficant. These results confirm the simpler linera distance is most appropriate.

# 1. distance-based Moran's eigenvector maps (equal to PCNM)  
#threshold for the minimum spanning tree set to 87 (max distance among plots
#within a local site. Setting this value to null fails)
med.dbmem1 <- dbmem(med.xy, thresh = 95.685706, MEM.autocor = "positive")
# Print the (n-1) non-null eigenvalues
attributes(med.dbmem1)$values # 3 positive eigenvalues
# Plot the associated spatial weighting matrix
s.label(med.xy)
# Plot maps of the first 3 dbMEM eigenfunctions
s.value(med.xy, med.dbmem1[,1:3])
# Compute and test associated Moran's I values
# Eigenvalues are proportional to Moran's I
test <- moran.randtest(med.dbmem1, nrepet = 99)
plot(test$obs, attr(med.dbmem1, "values"), xlab = "Moran's I", ylab = "Eigenvalues")
# Decreasing values of Moran's I for the successive MEM.
# The red line is the expected value of Moran's I under H0.
plot(test$obs, xlab="MEM rank", ylab="Moran's I")
abline(h=-1/(nrow(med.xy) - 1), col="red")
#convert dbMEM result to matrix for RDA
med.dbMEM <- as.matrix(med.dbmem1)


# 2. Run the global dbMEM analysis on the *detrended* data
# Is there a linear trend? If yes, detrend the data before assessing more
# complicated spatial structures.
anova(rda(traits.med.sc, med.xy), step=1000)	# Result: no significant trend (0.092)


cap.mod.med <-capscale(traits.eu.med ~ PC1+PC2+PC3+PC4+PC5+PC6, pc.scores.med,
                       dist="euclidean", add =TRUE)
anova(cap.mod.med, by = "margin")
RsquareAdj(cap.mod.med)$adj.r.squared
plot(cap.mod.med)

# comute RDA and variance partitioning using all 7 environmetal PCs and linear coordinates 
med.varpart.linearxy <- varpart(traits.med.sc, ~ PC1+PC2+PC3+PC4+PC5+PC6, med.xy, data = pc.scores.med)
med.varpart.linearxy
med.varpart.linearxy.eu <- varpart(traits.eu.med, ~ PC1+PC2+PC3+PC4+PC5+PC6, med.xy, data = pc.scores.med)
med.varpart.linearxy.eu
med.varpart.linearMEM.eu <- varpart(traits.eu.med, ~ PC1+PC2+PC3+PC4+PC5+PC6, med.dbMEM, data = pc.scores.med)
med.varpart.linearMEM.eu




# High Productivity Site --------------------------------------------------------
# pca on all abiotic environmental variables
env.high <- data.high[, -c(4:12,32,36:37,39,40,43,46)] 
high.pca <- prcomp(env.high[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(high.pca, type = "l")
summary(high.pca)
high.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(high.pca, scale = 0)
pc.scores.high <- as.data.frame(high.pca$x[ ,1:7])
pc.min.high <- abs(min(pc.scores.high[ ,1:6])) 
env.het.data.high <- pc.scores.high + pc.min.high
env.het.high <- vegdist(env.het.data.high[ ,1:6], method = "euclidean")
mod_env.het.high <- betadisper(env.het.high, rep(1, length(data.high$plot)),
                               type = "median", bias.adjust=FALSE)

# calculate fractions of variance explained using distance-based redundancy 
# analysis (dbRDA). dbRDA takes distance matrix, then makes an ordination 
# (i.e., principal coordinates analysis), then a constrained ordination in which 
# ordination axes are a linear function of environmental and spatial variables.
# function "capscale" (using a formula with only a constant on the left hand side) 
# takes the distance and creates a multivariate space that preserves those distances 
# (i.e., performs principal coordinates analysis); the main message: the result 
# of principal coordinate analysis (the object that below is named "traits.model") 
# is the response variable used for variance partitioning using function 'varpart'.
traits.model.high <- capscale(traits.high[, 12:17] ~ 1)
traits.high.sc = scores(traits.model.high, display = c("wa"), scaling =1)
traits.eu.high <- vegdist(traits.high[, 12:17], method = "euclidean")

# There are multiple ways to model the spatial component. Linear distance of 
# coordinates and Moran's eigenvector maps (PCNM= dbMEM). Because the study 
# was designed and sampled in a way to minimize and quilize space-env. correlations
# across the three producitivity levels, linear distance is most appropriate (there
# is no reason to expect complex spatial patterns at the scales we are addressing).
# That being said, the code behigh using both linear distance and dbMEM (PCNM) with
# forward model selection just to check. As expected, both approaches give similar
# almost identical results. Also, after linear detrending the data, no dbMEM axes
# are signficant. These results confirm the simpler linera distance is most appropriate.

# 1. distance-based Moran's eigenvector maps (equal to PCNM)  
#threshold for the minimum spanning tree set to 87 (max distance among plots
#within a local site. Setting this value to null fails)
high.dbmem1 <- dbmem(high.xy, thresh = 94.225987, MEM.autocor = "positive")
# Print the (n-1) non-null eigenvalues
attributes(high.dbmem1)$values # 3 positive eigenvalues
# Plot the associated spatial weighting matrix
s.label(high.xy)
# Plot maps of the first 3 dbMEM eigenfunctions
s.value(high.xy, high.dbmem1[,1:3])
# Compute and test associated Moran's I values
# Eigenvalues are proportional to Moran's I
test <- moran.randtest(high.dbmem1, nrepet = 99)
plot(test$obs, attr(high.dbmem1, "values"), xlab = "Moran's I", ylab = "Eigenvalues")
# Decreasing values of Moran's I for the successive MEM.
# The red line is the expected value of Moran's I under H0.
plot(test$obs, xlab="MEM rank", ylab="Moran's I")
abline(h=-1/(nrow(high.xy) - 1), col="red")
#convert dbMEM result to matrix for RDA
high.dbMEM <- as.matrix(high.dbmem1)


# 2. Run the global dbMEM analysis on the *detrended* data
# Is there a linear trend? If yes, detrend the data before assessing more
# complicated spatial structures.
anova(rda(traits.high.sc, high.xy), step=1000)	# Result: no significant trend (0.001)


cap.mod.high <-capscale(traits.eu.high ~ PC1+PC2+PC3+PC4+PC5+PC6, pc.scores.high,
                        dist="euclidean", add =TRUE)
anova(cap.mod.high, by = "margin")
RsquareAdj(cap.mod.high)$adj.r.squared
plot(cap.mod.high)

# comute RDA and variance partitioning using all 7 environmetal PCs and linear coordinates 
high.varpart.linearxy <- varpart(traits.high.sc, ~ PC1+PC2+PC3+PC4+PC5+PC6, high.xy, data = pc.scores.high)
high.varpart.linearxy
high.varpart.linearxy <- varpart(traits.eu.high, ~ PC1+PC2+PC3+PC4+PC5+PC6, high.xy, data = pc.scores.high)
high.varpart.linearxy
high.varpart.linearxy.eu <- varpart(traits.eu.high, ~ PC1+PC2+PC3+PC4+PC5+PC6, high.dbMEM, data = pc.scores.high)
high.varpart.linearxy.eu


# box plots of differences in environmental heterogeneity (7 PCs)
par(mfrow = c(1,1))
boxplot(mod_env.het.low$distances, mod_env.het.med$distances, 
        mod_env.het.high$distances, main = "Environment", 
        ylab = "Distance-to-median", ylim = c(0, 1), 
        xlab = "Productivity", names = c("Low", "Medium", "High"),
        col = c("grey90", "grey90", "grey90")) 


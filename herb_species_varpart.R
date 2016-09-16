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

setwd("C:/Users/Chris/Desktop/Beta-productivity Project/Data/ReadyForR")


#-------------------------------------------------------------------------------
# 1. Import data and prepare for analyses
#-------------------------------------------------------------------------------

# Load plant community data, convert plot to character, and order by plot.
species <- read.csv("HerbData2015_Summer_wide_07012016.csv")
species$plot <- as.character(species$plot)
species <- species[order(species$plot),]

# Load environmental data, convert plot to character, and order by plot.
data <- read.csv("Environmental_Data_All_08102016.csv")
data$plot <- as.character(data$plot)
data <- data[order(data$plot),]

# extract data for specific to each soil productivity level
species.low <- species[species$productivity == "low", ]
species.med <- species[species$productivity == "med", ]
species.high <- species[species$productivity == "high", ]

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

# Transform species composition data using hellinger distance. THis is necessary
# for the RDA (equivilent to dbRDA with hellinger distance). Note: I also checked
# dbRDA using Bray0Curtis but this led to poorer models (less explained var)
low.h <- decostand(species.low[, 4:68], "hellinger")

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
anova(rda(low.h, low.xy), step=1000)	# Result: significant trend
# Computation of linearly detrended data
low.h.det <- resid(lm(as.matrix(low.h) ~ ., data=low.xy))
(low.PCNM.rda <- rda(low.h.det, low.dbMEM))
anova(low.PCNM.rda, step=1000) # after linear detrending, gloal dbMEM not signficant,
# no need for model selection of further analysis with db MEM, using linear xy coords

# comute RDA and variance partitioning using all 7 environmetal PCs and linear coordinates 
low.RDA <- rda(low.h, pc.scores.low[,1:6], low.xy)
anova(low.RDA)
low.varpart.linearxy <- varpart(low.h, ~ PC1+PC2+PC3+PC4+PC5+PC6, low.xy, data = pc.scores.low)
low.varpart.linearxy
#adjusted R-squared
#a   = 0.17204
#a+b = 0.17612 (0.37499 r2)
#b+c = 0.04518
#c   = 0.04111
#total explained = 21.723 (0.22840 with PC1-6)
#total unexplained = 0.78277

# instead of RDA with hellinger transform, consider dbRDA with bray-curtis

low.bc <- vegdist(species.low[, 4:68], method = "bray")
low.dbRDA.BC <- capscale(low.bc ~ PC1+PC2+PC3+PC4+PC5+PC6 + Condition(scores(low.xy)), data = pc.scores.low)
anova(low.dbRDA.BC)
low.varpart.bc <- varpart(low.bc, ~ PC1+PC2+PC3+PC4+PC5+PC6, low.xy, data = pc.scores.low)
low.varpart.bc
#unexplained 0.84

low.mh <- vegdist(species.low[, 4:68], method = "horn")
low.dbRDA.MH <- capscale(low.mh ~ PC1+PC2+PC3+PC4+PC5+PC6 + Condition(scores(low.xy)), data = pc.scores.low)
anova(low.dbRDA.MH)
low.varpart.MH <- varpart(low.mh, ~ PC1+PC2+PC3+PC4+PC5+PC6, low.xy, data = pc.scores.low)
low.varpart.MH
#unexplained 0.84


# calculate spatial distance matrix, then test for space-environment correlation
# at the whole plot scale
low.spatialdists <- dist(low.xy)
low.envdists <- dist(pc.scores.low[,1:6])
mantel(low.spatialdists, low.envdists, method="pearson", permutations=99)
# Mantel statistic r: 0.4484, Significance: 0.01 


# Med Productivity Site --------------------------------------------------------
# pca on all abiotic environmental variables
env.med <- data.med[, -c(4:12,32,36:37,39,40,43,46)] 
med.pca <- prcomp(env.med[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(med.pca, type = "l")
summary(med.pca)
med.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(med.pca, scale = 0)
pc.scores.med <- as.data.frame(med.pca$x[ ,1:7])

# Transform species composition data using hellinger distance. THis is necessary
# for the RDA (equivilent to dbRDA with hellinger distance). Note: I also checked
# dbRDA using Bray0Curtis but this led to poorer models (less explained var)
med.h <- decostand(species.med[, 4:68], "hellinger")

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
anova(rda(med.h, med.xy), step=1000)	# Result: close significant trend (0.092)
# Computation of linearly detrended data
med.h.det <- resid(lm(as.matrix(med.h) ~ ., data=med.xy))
(med.PCNM.rda <- rda(med.h.det, med.dbMEM))
anova(med.PCNM.rda, step=1000) # after linear detrending, gloal dbMEM not signficant,
# no need for model selection of further analysis with db MEM, using linear xy coords

# comute RDA and variance partitioning using all 7 environmetal PCs and linear coordinates 
med.RDA <- rda(med.h, pc.scores.med[1:6], med.xy)
anova(med.RDA)
med.varpart.linearxy <- varpart(med.h, ~ PC1+PC2+PC3+PC4+PC5+PC6, med.xy, data = pc.scores.med)
med.varpart.linearxy
#adjusted R-squared
#a   = -0.02456
#a+b = 0.01103 (0.24975 r2)
#b+c = 0.03218
#c   = -0.00341
#total explained = 0.00762 (0.22840 with PC1-6)
#total unexplained = 0.99238

# instead of RDA with hellinger transform, consider dbRDA with bray-curtis

med.bc <- vegdist(species.med[, 4:68], method = "bray")
med.dbRDA.BC <- capscale(med.bc ~ PC1+PC2+PC3+PC4+PC5+PC6 + Condition(scores(med.xy)), data = pc.scores.med)
anova(med.dbRDA.BC)
med.varpart.bc <- varpart(med.bc, ~ PC1+PC2+PC3+PC4+PC5+PC6, med.xy, data = pc.scores.med)
med.varpart.bc

med.mh <- vegdist(species.med[, 4:68], method = "horn")
med.dbRDA.MH <- capscale(med.mh ~ PC1+PC2+PC3+PC4+PC5+PC6 + Condition(scores(med.xy)), data = pc.scores.med)
anova(med.dbRDA.MH)
med.varpart.MH <- varpart(med.mh, ~ PC1+PC2+PC3+PC4+PC5+PC6, med.xy, data = pc.scores.med)
med.varpart.MH
#unexplained 0.95 (performed better,but still low)

# calculate spatial distance matrix, then test for space-environment correlation
# at the whole plot scale
med.spatialdists <- dist(med.xy)
med.envdists <- dist(pc.scores.med[,1:6])
mantel(med.spatialdists, med.envdists, method="pearson", permutations=99)
# Mantel statistic r: 0.5157, Significance: 0.01 


# High Productivity Site --------------------------------------------------------
# pca on all abiotic environmental variables
env.high <- data.high[, -c(4:12,32,36:37,39,40,43,46)] 
high.pca <- prcomp(env.high[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(high.pca, type = "l")
summary(high.pca)
high.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(high.pca, scale = 0)
pc.scores.high <- as.data.frame(high.pca$x[ ,1:7])

# Transform species composition data using hellinger distance. THis is necessary
# for the RDA (equivilent to dbRDA with hellinger distance). Note: I also checked
# dbRDA using Bray0Curtis but this led to poorer models (less explained var)
high.h <- decostand(species.high[, 4:68], "hellinger")

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
anova(rda(high.h, high.xy), step=1000)	# Result: significant trend (0.001)
# Computation of linearly detrended data
high.h.det <- resid(lm(as.matrix(high.h) ~ ., data=high.xy))
(high.PCNM.rda <- rda(high.h.det, high.dbMEM))
anova(high.PCNM.rda, step=1000) # after linear detrending, gloal dbMEM not signficant,
# no need for model selection of further analysis with db MEM, using linear xy coords

# comute RDA and variance partitioning using all 7 environmetal PCs and linear coordinates 
high.RDA <- rda(high.h, pc.scores.high, high.xy)
anova(high.RDA)
high.varpart.linearxy <- varpart(high.h, ~ PC1+PC2+PC3+PC4+PC5+PC6, high.xy, data = pc.scores.high)
high.varpart.linearxy
#adjusted R-squared
#a   = 0.08363
#a+b = 0.17307 (0.37268 r2)
#b+c = 0.08953
#c   = 0.00009
#total explained = 0.17316 (0.22840 with PC1-6)
#total unexplained = 0.82684

# instead of RDA with hellinger transform, consider dbRDA with bray-curtis

high.bc <- vegdist(species.high[, 4:68], method = "bray")
high.dbRDA.BC <- capscale(high.bc ~ PC1+PC2+PC3+PC4+PC5+PC6 + Condition(scores(high.xy)), data = pc.scores.high)
anova(high.dbRDA.BC)
high.varpart.bc <- varpart(high.bc, ~ PC1+PC2+PC3+PC4+PC5+PC6, high.xy, data = pc.scores.high)
high.varpart.bc

high.mh <- vegdist(species.high[, 4:68], method = "horn")
high.dbRDA.MH <- capscale(high.mh ~ PC1+PC2+PC3+PC4+PC5+PC6 + Condition(scores(high.xy)), data = pc.scores.high)
anova(high.dbRDA.MH)
high.varpart.MH <- varpart(high.mh, ~ PC1+PC2+PC3+PC4+PC5+PC6, high.xy, data = pc.scores.high)
high.varpart.MH


# calculate spatial distance matrix, then test for space-environment correlation
# at the whole plot scale
high.spatialdists <- dist(high.xy)
high.envdists <- dist(pc.scores.high[,1:6])
mantel(high.spatialdists, high.envdists, method="pearson", permutations=99)
# Mantel statistic r:  0.43,  Significance: 0.01  



#Old working code for PCNMs and spatial model selection  
#################################################################################
#############################################################################
############################################################################

# PCNM
xy.d1 <- dist(low.xy)
spanning <- spantree(xy.d1)
dmin <- max(spanning$dist)
dmin <- 87.290463
# Truncate the distance matrix
xy.d1[xy.d1 > dmin] <- 4*dmin
# PCoA of truncated distance matrix
xy.PCoA <- cmdscale(xy.d1, k=nrow(low.xy)-1, eig=TRUE)
# Count the positive eigenvalues (PCNM with positive AND negative spatial
# correlation)
(nb.ev <- length(which(xy.PCoA$eig > 0.000000001)))
# Construct a data frame containing the PCNM variables
low.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(low.xy), 1:nb.ev])

## 2. Run the global PCNM analysis on the *detrended* mite data
# Is there a linear trend in the mite data?
anova(rda(low.h, low.xy), step=1000)	# Result: significant trend
# Computation of linearly detrended mite data
low.h.det <- resid(lm(as.matrix(low.h) ~ ., data=low.xy))
(low.PCNM.rda <- rda(low.h.det, low.PCNM))
#(low.PCNM.rda <- rda(low.h, low.PCNM)) # if not detrended
anova(low.PCNM.rda, step=1000)

(low.R2a <- RsquareAdj(low.PCNM.rda)$adj.r.squared) 
(low.PCNM.fwd <- forward.sel(low.h.det, as.matrix(low.PCNM), 
                             adjR2thresh=low.R2a))
# According to the R2a criterion, if we retain PCNM 5 we get a model with 
# a R2adj slightly higher than that of the complete model. This slight 
# excess is not too serious, however.
(nb.sig.PCNM <- nrow(low.PCNM.fwd))	# Number of signif. PCNM
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(low.PCNM.fwd[,2]))
# Write the significant PCNMs to a new object
PCNM.red <- low.PCNM[,c(PCNM.sign)]

## 4. New PCNM analysis with 10 significant PCNM variables
##    Adjusted R-square after forward selection: R2adj=0.2713
(low.PCNM.rda2 <- rda(low.h.det ~ PCNM.red)) # if detrended
(low.fwd.R2a <- RsquareAdj(low.PCNM.rda2)$adj.r.squared)
anova(low.PCNM.rda2, step=1000)
(axes.test <- anova(low.PCNM.rda2, by="axis", step=1000))
# Number of significant axes
(nb.ax <- length(which(axes.test[,4] <= 0.05)))

low.PCNM.axes <- scores(low.PCNM.rda2, choices=c(1:nb.ax), display="lc", 
                        scaling=1)
windows(title="PCNM analysis of low data", 8, 6)
par(mfrow=c(1,1))
s.value(low.xy, low.PCNM.axes[,1])	# ade4 function: s.value

# Interpreting the spatial variation: regression of the significant
# canonical axes on the environmental variables
shapiro.test(resid(lm(low.PCNM.axes[,1] ~ ., data=pc.scores.low))) # Normality test
low.PCNM.axis1.env <- lm(low.PCNM.axes[,1]~., data=pc.scores.low)
summary(low.PCNM.axis1.env)

# Maps of the 10 significant PCNM variables
windows(title="10 PCNM variables - low")
par(mfrow=c(1,1))
for(i in 1:ncol(PCNM.red))
{
  s.value(low.xy, PCNM.red, sub=PCNM.sign, csub=2)
}

# variance partitioning with linear spatial coordinates

low.varpart <- varpart(low.h, ~ ., low.xy, data = pc.scores.low)
low.varpart

# variance partitioning with PCNM spatial variables
low.varpart.pcnm <- varpart(low.h, ~ ., PCNM.red, data = pc.scores.low)
low.varpart.pcnm

# follow this code for forward selection on pcnm with regular (not untrended) data


## 3. Test and forward selection of PCNM variables
# Run the global PCNM analysis on the *undetrended* mite data
low.undet.PCNM.rda <- rda(low.h, low.PCNM)
anova(low.undet.PCNM.rda, step=1000)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the PCNM variables
(low.undet.PCNM.R2a <- RsquareAdj(low.undet.PCNM.rda)$adj.r.squared)
(low.undet.PCNM.fwd <- forward.sel(low.h, as.matrix(low.PCNM), 
                                   adjR2thresh=0.1))
# Number of significant PCNMs
(nb.sig.PCNM <- nrow(low.undet.PCNM.fwd))
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(low.undet.PCNM.fwd$order))
# Write the significant PCNMs to a new object
PCNM.red <- low.PCNM[,c(PCNM.sign)]

# variance partitioning with PCNM spatial variables
low.varpart.pcnm.sel <- varpart(low.h, ~ ., PCNM.red, data = env.red)
low.varpart.pcnm.sel

windows(title="Species - environment - space variation partitioning", 12, 6)
par(mfrow=c(1,2))
showvarparts(2)
plot(low.varpart, digits=2)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(low.h, pc.scores.low, low.xy))











env.med <- data.med[, -c(4:12,32,36:37,39,40,43,46)]
med.pca <- prcomp(env.med[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(med.pca, type = "l")
summary(med.pca)
med.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(med.pca, scale = 0)
pc.scores.med <- as.data.frame(med.pca$x[ ,1:7])

med.h <- decostand (species.med[, 4:68], "hellinger")

# PCNM
xy.d1 <- dist(med.xy)
spanning <- spantree(xy.d1)
dmin <- max(spanning$dist)
# Truncate the distance matrix
xy.d1[xy.d1 > dmin] <- 4*dmin
# PCoA of truncated distance matrix
xy.PCoA <- cmdscale(xy.d1, k=nrow(med.xy)-1, eig=TRUE)
# Count the positive eigenvalues (PCNM with positive AND negative spatial
# correlation)
(nb.ev <- length(which(xy.PCoA$eig > 0.0000001)))
# Construct a data frame containing the PCNM variables
med.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(med.xy), 1:nb.ev])

## 2. Run the global PCNM analysis on the *detrended* mite data
# Is there a linear trend in the mite data?
anova(rda(med.h, med.xy), step=1000)	# Result: significant trend
# Computation of linearly detrended mite data
med.h.det <- resid(lm(as.matrix(med.h) ~ ., data=med.xy))
(med.PCNM.rda <- rda(med.h.det, med.PCNM))
anova(med.PCNM.rda, step=1000)

(med.R2a <- RsquareAdj(med.PCNM.rda)$adj.r.squared)
(med.PCNM.fwd <- forward.sel(med.h.det, as.matrix(med.PCNM), 
                             adjR2thresh=0.1))
# According to the R2a criterion, if we retain PCNM 5 we get a model with 
# a R2adj slightly higher than that of the complete model. This slight 
# excess is not too serious, however.
(nb.sig.PCNM <- nrow(med.PCNM.fwd))	# Number of signif. PCNM
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(med.PCNM.fwd[,2]))
# Write the significant PCNMs to a new object
PCNM.red <- med.PCNM[,c(PCNM.sign)]

## 4. New PCNM analysis with 10 significant PCNM variables
##    Adjusted R-square after forward selection: R2adj=0.2713
(med.PCNM.rda2 <- rda(med.h.det ~ PCNM.red))
(med.fwd.R2a <- RsquareAdj(med.PCNM.rda2)$adj.r.squared)
anova(med.PCNM.rda2, step=1000)
(axes.test <- anova(med.PCNM.rda2, by="axis", step=1000))
# Number of significant axes
(nb.ax <- length(which(axes.test[,4] <= 0.05)))

med.PCNM.axes <- scores(med.PCNM.rda2, choices=c(1:nb.ax), display="lc", 
                        scaling=1)
windows(title="PCNM analysis of med data", 8, 6)
par(mfrow=c(1,1))
s.value(med.xy, med.PCNM.axes[,1])	# ade4 function: s.value

# Interpreting the spatial variation: regression of the significant
# canonical axes on the environmental variables
shapiro.test(resid(lm(med.PCNM.axes[,1] ~ ., data=pc.scores.med))) # Normality test
med.PCNM.axis1.env <- lm(med.PCNM.axes[,1]~., data=pc.scores.med)
summary(med.PCNM.axis1.env)

# Maps of the 10 significant PCNM variables
windows(title="10 PCNM variables - med")
par(mfrow=c(1,1))
for(i in 1:ncol(PCNM.red))
{
  s.value(med.xy, PCNM.red, sub=PCNM.sign, csub=2)
}

# variance partitioning with linear spatial coordinates
med.h <- decostand (species.med[, 4:68], "hellinger")
med.varpart <- varpart(med.h, ~ ., med.xy, data = pc.scores.med)
med.varpart

# variance partitioning with PCNM spatial variables
med.varpart.pcnm <- varpart(med.h, ~ ., PCNM.red, data = pc.scores.med)
med.varpart.pcnm


windows(title="Species - environment - space variation partitioning", 12, 6)
par(mfrow=c(1,2))
showvarparts(2)
plot(med.varpart, digits=2)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(med.h, pc.scores.med, med.xy))










med.varpart <- varpart(med.h, ~ ., med.xy, data = pc.scores.med)
med.varpart

windows(title="Species - environment - space variation partitioning", 12, 6)
par(mfrow=c(1,2))
showvarparts(2)
plot(med.varpart, digits=2)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(med.h, pc.scores.med, med.xy))


env.high <- data.high[, -c(4:12,32,36:37,39,40,43,46)]
high.pca <- prcomp(env.high[ , c(4:30)], center = TRUE, scale. = TRUE) 
plot(high.pca, type = "l")
summary(high.pca)
high.pca$rotation[, 1:7] #examine the contribution of each variable to each principal component
biplot(high.pca, scale = 0)
pc.scores.high <- as.data.frame(high.pca$x[ ,1:7])

high.h <- decostand (species.high[, 4:68], "hellinger")

## 1a. Construct the matrix of PCNM variables step by step...
xy.d1 <- dist(high.xy)
spanning <- spantree(xy.d1)
dmin <- max(spanning$dist)
# Truncate the distance matrix
xy.d1[xy.d1 > dmin] <- 4*dmin
# PCoA of truncated distance matrix
xy.PCoA <- cmdscale(xy.d1, k=nrow(high.xy)-1, eig=TRUE)
# Count the positive eigenvalues (PCNM with positive AND negative spatial
# correlation)
(nb.ev <- length(which(xy.PCoA$eig > 0.0000001)))
# Construct a data frame containing the PCNM variables
high.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(high.xy), 1:nb.ev])

## 3. Test and forward selection of PCNM variables
# Run the global PCNM analysis on the *undetrended* high data
high.undet.PCNM.rda <- rda(high.h, high.PCNM)
anova(high.undet.PCNM.rda, step=1000)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the PCNM variables
(high.undet.PCNM.R2a <- RsquareAdj(high.undet.PCNM.rda)$adj.r.squared)
(high.undet.PCNM.fwd <- forward.sel(high.h, as.matrix(high.PCNM), 
                                    adjR2thresh=high.undet.PCNM.R2a))
# Number of significant PCNMs
(nb.sig.PCNM <- nrow(high.undet.PCNM.fwd))
# Identity of significant PCNMs in increasing order
(PCNM.sign <- sort(high.undet.PCNM.fwd$order))
# Write the significant PCNMs to a new object
PCNM.red <- high.PCNM[,c(PCNM.sign)]


high.PCNM.vegan <- pcnm(dist(high.xy))
(dminv <- high.PCNM.vegan$threshold)
(nb.evv <- length(which(high.PCNM.vegan$values > 0.0000001)))
high.PCNMv <- as.data.frame(high.PCNM.vegan$vectors)


high.varpart <- varpart(high.h, ~ ., high.xy, data = pc.scores.high)
high.varpart

high.varpart <- varpart(high.h, ~ ., high.PCNMv, data = pc.scores.high)
high.varpart

high.varpart <- varpart(high.h, ~ ., PCNM.red, data = pc.scores.high)
high.varpart


windows(title="Species - environment - space variation partitioning", 12, 6)
par(mfrow=c(1,2))
showvarparts(2)
plot(high.varpart, digits=2)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(high.h, pc.scores.high, low.PCNMv))





##----------------------------------------------------------------------------##
## Author: Christopher P. Catano                                              ## 
## Created: July 06, 2016                                                     ## 
##                                                                            ## 
## Discription: Null model analysis to determine how the processes driving    ##
##              beta diversity patterns change along productivty levels in    ##
##              the TRC Forest Dynamics Plot, MO                              ##
##                                                                            ##
## Input:  Herbaceous plant data (HerbData2015_Summer_wide_07012016.csv)      ##
##                                                                            ##
##                                                                            ##
## Output: TBD                                                                ##
##----------------------------------------------------------------------------##


#-------------------------------------------------------------------------------
# Details:
# This script runs a null model that simulates data by randomly shuffling 
# individuals among plots while preserving species relative abundance and the 
# number of individuals per plot; thus creating a matrix of simulated species 
# abundances and simulated Bray-Curtis distances (in species composition) b/w 
# each pair of plots. Then, for each pair of plots, the script calculates an 
# effect size that measures the departure of observed from simulated Bray-curtis 
# distances. The scripts uses the effect size # for all pairs of plots as the 
# response in a distance-based redundancy analysis (db-RDA), respective forward 
# mode selection following Blanchet et al. (2008), and variance partitioning 
# between environmental and spatial matrix. Finally, the script calculates 
# bootstrap distributions that will be used (in another script) to test 
# predictions about differences between variance fractions.
#-------------------------------------------------------------------------------


  rm(list = ls())

  library(reshape2)
  library(plyr)
  library(vegan)
  library(ggplot2)

  setwd("C:/Users/Chris/Desktop/Beta-productivity Project/Data/ReadyForR")


#-------------------------------------------------------------------------------
# 1. Import data and prepare for analyses
#-------------------------------------------------------------------------------

# Load herbaceous plant census data.
  dataSummer_wide <- read.csv("HerbData2015_Summer_wide_07012016.csv")
  data <- dataSummer_wide[, -c(1:3)]
  env <- dataSummer_wide[, c(1:3)]

# extract data for specific to each soil productivity level
  dataSummer_Low <- dataSummer_wide[dataSummer_wide$productivity == "low", ]
  dataSummer_Med <- dataSummer_wide[dataSummer_wide$productivity == "med", ]
  dataSummer_High <- dataSummer_wide[dataSummer_wide$productivity == "high", ]


#-------------------------------------------------------------------------------
# 2. create ordination using non-metric multidimensional scaling (NMDS) to 
#    visualize differences in species composition among the three productivity
#    levels (low = south-west facing slopes, med = north to east facing slopes,
#            high = valleys)
#    Reference: http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
#-------------------------------------------------------------------------------  

# NMDS 
  sol <- metaMDS(data, distance = "bray", k = 2, trymax = 100)

# Make a new data frame, add produtivity and site information (used for coloring
# and shape of points  
  NMDS <- data.frame(x = sol$point[, 1], y = sol$point[, 2], 
                     productivity = as.factor(env[, 2]), 
                     site = as.factor(env[, 3]))

# Get spread of points based on productivity level 
  plot(NMDS)
  ord <- ordiellipse(sol, as.factor(env[, 2]), display = "sites", kind = "sd", 
                     conf = 0.95, label = T)

# Data frame df_ell contains values to show ellipses. It is calculated with 
# function veganCovEllipse which is hidden in vegan package. This function is 
# applied to each level of NMDS (group) and it uses also function cov.wt to 
# calculate covariance matrix. 
  veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

# Generate ellipse points (not that ellipse created for med productivity inflates
# the percieved overlap with high productivity due to two outlying plots) 
  df_ell <- data.frame()
  for(g in levels(NMDS$productivity)){
    if(g != "" && (g %in% names(ord))){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$productivity == g, ],
                                                       veganCovEllipse(ord[[g]]$cov, 
                                                                       ord[[g]]$center,
                                                                       ord[[g]]$scale)))
                                    , productivity = g))
    }
  }

  head(df_ell)

# Generate mean values from nMDS plot grouped by productivity leel 
  NMDS.mean <- aggregate(NMDS[, 1:2], list(group = NMDS$productivity), mean)

# Plot nMDS  
  shape_values <- seq(1, 9)
  p <- ggplot(data = NMDS, aes(x, y, colour = productivity)) 

  #p <- p + geom_path(data=df_ell, aes(x = NMDS1, y = NMDS2), size = 1, linetype = 2)
  p <- p + geom_point(aes(size = 2)) + 
       scale_shape_manual(values = shape_values)  
  p <- p + annotate("text", x = NMDS.mean$x, y = NMDS.mean$y, 
                     label = NMDS.mean$group, size = 4)+
       theme_bw() 
  p


# Calculate mean community size, mean richness, species pool richness and 
# relative abundance for each productivity level
  gamma <- specnumber(data, groups = env$productivity) 
  richness <- specnumber(data, groups = env$plot) 
  mean_rich <- cbind(richness, env$productivity)
  boxplot(mean_rich[,1]~ mean_rich[,2], main = "Mean plot richness")
  text(2, 18, "1=high, 2=low, 3=med")

# Calculate mean plot richness for each treatment
  alpha_mean <- tapply(specnumber(data), env$productivity, mean)

# Calculate mean community size (# of individuals/plot) by treatment
  comsize <- as.matrix(rowSums(data))
  mean_size <- cbind(comsize, env$productivity)
  boxplot(mean_size[,1]~ mean_size[,2], main = "Mean community size")
  text(2, 130, "1=high, 2=low, 3=med")

# Calculate regional abundance distribution (evenness for each treatment)
  newData <- cbind(data, env)
  regAbund <- aggregate(. ~ productivity, data = newData, FUN = sum)
  H <- diversity(regAbund[, -1])
  J <- H/log(gamma)
  SP_shannon <- exp(H)


#-------------------------------------------------------------------------------
# 3. Calculate observed beta-diversity and test for differences among three
#    productivity (soil resource availability) levels: Low, Medium, High
#-------------------------------------------------------------------------------

# Calculate community dissimilarities using Jaccard (presence-absence), 
# Bray-Curtis, and Morisita-Horn indices.
  Beta_jc <- vegdist(data, method = "jaccard", binary = TRUE)
  Beta_bc <- vegdist(data, method = "bray")
  Beta_mh <- vegdist(data, method = "horn") 

# Calculate distance to centroids for each prouctivity level using the bootstrapped
# mean distance to spatial medians (conserves appropriatet Type I error rates 
# compared to distance to cetroid)
  mod_jc <- betadisper(Beta_jc, env$productivity, type = "median", bias.adjust=FALSE)
  mod_bc <- betadisper(Beta_bc, env$productivity, type = "median", bias.adjust=FALSE)
  mod_mh <- betadisper(Beta_mh, env$productivity, type = "median", bias.adjust=FALSE)

# Permutation test for F (more conservative than ANOVA by relaxing assumption 
# that data are derived from the normal distribution)
  permutest(mod_jc, pairwise = TRUE, permutations = 9999)
  permutest(mod_bc, pairwise = TRUE, permutations = 9999)
  permutest(mod_mh, pairwise = TRUE, permutations = 9999)

# Plot distance to centroids (spatial medians) for Bray-Curtis dissimilarity
  par(mfrow = c(1,1))
  plot(mod_bc, axes = c(1,2), hull = TRUE, main = "Distance to centroid")
  text(-0.4, 0.0, "high", cex = 1)
  text(-0.2, 0.4, "med", cex = 1)
  text(0.2, -0.1, "low", cex = 1)

  
#-------------------------------------------------------------------------------
# 4. Plots of observed dissimilarity and distance-to-centroid as calculated from 
#    spatial medians.
#------------------------------------------------------------------------------- 

# Low Productivity Plots
  jcLow <- vegdist(dataSummer_Low[ , 4:68], method = "jaccard", binary = TRUE)
  bcLow <- vegdist(dataSummer_Low[ , 4:68], method = "bray")
  mhLow <- vegdist(dataSummer_Low[ , 4:68], method = "horn")

# Medium Productivity Plots
  jcMed <- vegdist(dataSummer_Med[ , 4:68], method = "jaccard", binary = TRUE)
  bcMed <- vegdist(dataSummer_Med[ , 4:68], method = "bray")
  mhMed <- vegdist(dataSummer_Med[ , 4:68], method = "horn")

# High Productivity Plots
  jcHigh <- vegdist(dataSummer_High[ , 4:68], method = "jaccard", binary = TRUE)
  bcHigh <- vegdist(dataSummer_High[ , 4:68], method = "bray")
  mhHigh <- vegdist(dataSummer_High[ , 4:68], method = "horn")

  par(mfrow = c(1,3))
# Jaccard   
  boxplot(jcLow, jcMed, jcHigh, ylim=c(0,1), main = "Jaccard", 	
          ylab = "Dissimilarity", names = c("Low", "Medium", "High"), 
          xlab = "Productivity level", col = c("grey90", "grey90", "grey90"))

# Bray-Curtis  
  boxplot(bcLow, bcMed, bcHigh, ylim=c(0,1), main = "Bray-Curtis", 	
          ylab = "Dissimilarity", names = c("Low", "Medium", "High"), 
          xlab = "Productivity level", col = c("red", "yellow", "green"))  

# Morisita-Horn  
  boxplot(mhLow, mhMed, mhHigh, ylim=c(0,1), main = "Morisita-Horn", 	
          ylab = "Dissimilarity", names = c("Low", "Medium", "High"), 
          xlab = "Productivity level", col = c("red", "yellow", "green"))  

# Jaccard beta-diversity  
  betadisp.jcLow <- betadisper(jcLow, rep(1, length(dataSummer_Low$plot)),
                               type = "median", bias.adjust=FALSE) 
  betadisp.jcMed <- betadisper(jcMed, rep(1, length(dataSummer_Med$plot)),
                               type = "median", bias.adjust=FALSE)
  betadisp.jcHigh <- betadisper(jcHigh, rep(1, length(dataSummer_High$plot)),
                                type = "median", bias.adjust=FALSE)

# Bray-Curtis beta-diversity  
  betadisp.bcLow <- betadisper(bcLow, rep(1, length(dataSummer_Low$plot)), 
                               type = "median", bias.adjust=FALSE)
  betadisp.bcMed <- betadisper(bcMed, rep(1, length(dataSummer_Med$plot)), 
                               type = "median", bias.adjust=FALSE)
  betadisp.bcHigh <- betadisper(bcHigh, rep(1, length(dataSummer_High$plot)), 
                                type = "median", bias.adjust=FALSE)

# Morisita-Horn beta-diversity  
  betadisp.mhLow <- betadisper(mhLow, rep(1, length(dataSummer_Low$plot)),
                               type = "median", bias.adjust=FALSE) 
  betadisp.mhMed <- betadisper(mhMed, rep(1, length(dataSummer_Med$plot)),
                               type = "median", bias.adjust=FALSE) 
  betadisp.mhHigh <- betadisper(mhHigh, rep(1, length(dataSummer_High$plot)),
                                type = "median", bias.adjust=FALSE)

  par(mfrow = c(1,3))
  boxplot(betadisp.jcLow$distance, betadisp.jcMed$distance, 
          betadisp.jcHigh$distance, main = "Jaccard beta-diversity", 
          ylab = "Distance-to-median", ylim = c(0, 1), 
          xlab = "Productivity", names = c("Low", "Medium", "High"),
          col = c("grey90", "grey90", "grey90")) 

  boxplot(betadisp.bcLow$distance, betadisp.bcMed$distance, 
          betadisp.bcHigh$distance, main = "Bray-Curtis beta-diversity", 
          ylab = "Distance-to-median", ylim = c(0, 1), 
          xlab = "Productivity", names = c("Low", "Medium", "High"),
          col = c("grey90", "grey90", "grey90")) 

  boxplot(betadisp.mhLow$distance, betadisp.mhMed$distance, 
          betadisp.mhHigh$distance, main = "Morisita-Horn beta-diversity", 
          ylab = "Distance-to-median", ylim = c(0, 1), 
          xlab = "Productivity", names = c("Low", "Medium", "High"),
          col = c("grey90", "grey90", "grey90"))


#-------------------------------------------------------------------------------
# 5. Create new data sets for null model. The new data sets will contain 
#    individuals--as estimated from species occurrences in 10x10-cm cells (as rows).
#-------------------------------------------------------------------------------  

# Load herbaceous plant census data in long form.
  dataSummer_long <- read.csv("HerbData2015_Summer_long_07012016.csv") 
  dataSummer_long$plot <- as.character(dataSummer_long$plot)
  dataSummer_long$species <- as.character(dataSummer_long$species)
  str(dataSummer_long)

# Subset to create a new data set for each productivity level
  dataLow1 <- dataSummer_long[dataSummer_long$productivity == "low", ]
  dataMed1 <- dataSummer_long[dataSummer_long$productivity == "med", ]
  dataHigh1 <- dataSummer_long[dataSummer_long$productivity == "high", ]

# Low productivity plots
  length(unique(dataLow1$species))  # 24 species
  length(unique(dataLow1$plot))     # 30 plots 

# Medium productivity plots
  length(unique(dataMed1$species))  # 27 species
  length(unique(dataMed1$plot))     # 30 plots 

# High productivity plots
  length(unique(dataHigh1$species))  # 44 species
  length(unique(dataHigh1$plot))     # 30 plots


#-------------------------------------------------------------------------------
# 6. Create empty vectors to store results from simulations
#-------------------------------------------------------------------------------

# Low productivity plots
  sim.species.abundance.low <- c()
  sim.plot.density.low <- c()
  sim.species.occurrence.low <- c()
  sim.plot.richness.low <- c()
  sim.bray.low <- c()  

# Medium productivity plots
  sim.species.abundance.med <- c()
  sim.plot.density.med <- c()
  sim.species.occurrence.med <- c()
  sim.plot.richness.med <- c()
  sim.bray.med <- c()  

# High productivity plots
  sim.species.abundance.high <- c()
  sim.plot.density.high <- c()
  sim.species.occurrence.high <- c()
  sim.plot.richness.high <- c()
  sim.bray.high <- c()  


#-------------------------------------------------------------------------------
# 7. Null models
#-------------------------------------------------------------------------------

# Set number of simulations (k). 
  k <- 2000

# The final vector of pairwise dissimilarities for analyses with N = 30 plots 
# will be: k*((30*30)-30)/2 = k*435

# Low Productivity--------------------------------------------------------------
  for (i in 1:k){
    a <- sample(dataLow1$plot, length(dataLow1$plot))
  
    sim.table.low <- tapply(rep(1,length(a)), 
                            list(a, as.vector(dataLow1$species)), sum)
  
    sim.table.low[is.na(sim.table.low)] <- 0
  
    sim.species.abundance.low <- append(sim.species.abundance.low,
                                        colSums(sim.table.low), 
                                        after = length(sim.species.abundance.low))
  
    sim.plot.density.low <- append(sim.plot.density.low, 
                                   rowSums(sim.table.low), 
                                   after = length(sim.plot.density.low))
  
    sim.species.occurrence.low <- append(sim.species.occurrence.low, 
                                         colSums(sim.table.low > 0), 
                                         after = length(sim.species.occurrence.low))
  
    sim.plot.richness.low <- append(sim.plot.richness.low,
                                    rowSums(sim.table.low > 0), 
                                    after = length(sim.plot.richness.low))
  
    sim.bray.low <- append(sim.bray.low, vegdist(sim.table.low, method = "bray"), 
                           after = length(sim.bray.low))
  }

  length(sim.bray.low)

# Med Productivity--------------------------------------------------------------
  for (i in 1:k){
    a <- sample(dataMed1$plot, length(dataMed1$plot))
  
    sim.table.med <- tapply(rep(1,length(a)), 
                            list(a, as.vector(dataMed1$species)), sum)
  
    sim.table.med[is.na(sim.table.med)] <- 0
  
    sim.species.abundance.med <- append(sim.species.abundance.med,
                                        colSums(sim.table.med), 
                                        after = length(sim.species.abundance.med))
  
    sim.plot.density.med <- append(sim.plot.density.med, 
                                   rowSums(sim.table.med), 
                                   after = length(sim.plot.density.med))
  
    sim.species.occurrence.med <- append(sim.species.occurrence.med, 
                                         colSums(sim.table.med > 0), 
                                         after = length(sim.species.occurrence.med))
  
    sim.plot.richness.med <- append(sim.plot.richness.med,
                                    rowSums(sim.table.med > 0), 
                                    after = length(sim.plot.richness.med))
  
    sim.bray.med <- append(sim.bray.med, vegdist(sim.table.med, method = "bray"), 
                           after = length(sim.bray.med))
  }

# High Productivity--------------------------------------------------------------
  for (i in 1:k){
    a <- sample(dataHigh1$plot, length(dataHigh1$plot))
  
    sim.table.high <- tapply(rep(1,length(a)), 
                             list(a, as.vector(dataHigh1$species)), sum)
  
    sim.table.high[is.na(sim.table.high)] <- 0
  
    sim.species.abundance.high <- append(sim.species.abundance.high,
                                         colSums(sim.table.high), 
                                         after = length(sim.species.abundance.high))
  
    sim.plot.density.high <- append(sim.plot.density.high, 
                                    rowSums(sim.table.high), 
                                    after = length(sim.plot.density.high))
  
    sim.species.occurrence.high <- append(sim.species.occurrence.high, 
                                          colSums(sim.table.high > 0), 
                                          after = length(sim.species.occurrence.high))
  
    sim.plot.richness.high <- append(sim.plot.richness.high,
                                     rowSums(sim.table.high > 0), 
                                     after = length(sim.plot.richness.high))
  
    sim.bray.high <- append(sim.bray.high, vegdist(sim.table.high, method = "bray"), 
                            after = length(sim.bray.high))
  }


#-------------------------------------------------------------------------------
# 8. Calculate values of mean simulated beta-diversity, save objects to text 
#    files, and create an object of class 'distance'.
#-------------------------------------------------------------------------------

# Low productivity--------------------------------------------------------------
  sim.bray.low.mean <- tapply(sim.bray.low, rep(1:length(bcLow), times = k), mean) 
  length(sim.bray.low.mean) 

#write.table(sim.bray.low.mean, "sim_bray_low_mean.txt", quote = T, sep = ",", 
#           row.names = F) 
  sim.bray.low.mean1 <- matrix(nrow = nrow(dataSummer_Low),  
                               ncol = nrow(dataSummer_Low)) 

  sim.bray.low.mean1[lower.tri(sim.bray.low.mean1)] <- as.vector(sim.bray.low.mean) 

  dimnames(sim.bray.low.mean1)[[1]] <- dimnames(dataSummer_Low)[[1]]

  sim.bray.low.mean2 <- as.dist(sim.bray.low.mean1)
  length(sim.bray.low.mean2)

# Medium productivity-----------------------------------------------------------
  sim.bray.med.mean <- tapply(sim.bray.med, rep(1:length(bcMed), times = k), mean)
  length(sim.bray.med.mean)

#write.table(sim.bray.med.mean, "sim_bray_med_mean.txt", quote = T, sep = ",", 
#           row.names = F)
  sim.bray.med.mean1 <- matrix(nrow = nrow(dataSummer_Med), 
                               ncol = nrow(dataSummer_Med))

  sim.bray.med.mean1[lower.tri(sim.bray.med.mean1)] <- as.vector(sim.bray.med.mean)

  dimnames(sim.bray.med.mean1)[[1]] <- dimnames(dataSummer_Med)[[1]]

  sim.bray.med.mean2 <- as.dist(sim.bray.med.mean1)
  length(sim.bray.med.mean2)  

# High productivity-------------------------------------------------------------
  sim.bray.high.mean <- tapply(sim.bray.high, rep(1:length(bcHigh), times = k), mean)
  length(sim.bray.high.mean)

#write.table(sim.bray.high.mean, "sim_bray_high_mean.txt", quote = T, sep = ",", 
#           row.names = F)
  sim.bray.high.mean1 <- matrix(nrow = nrow(dataSummer_High), 
                                ncol = nrow(dataSummer_High))

  sim.bray.high.mean1[lower.tri(sim.bray.high.mean1)] <- as.vector(sim.bray.high.mean)

  dimnames(sim.bray.high.mean1)[[1]] <- dimnames(dataSummer_High)[[1]]

  sim.bray.high.mean2 <- as.dist(sim.bray.high.mean1)
  length(sim.bray.high.mean2)   


#-------------------------------------------------------------------------------
# 9. Calculate distance-to-centroids for mean simulated values
#-------------------------------------------------------------------------------  

# Low productivity--------------------------------------------------------------
  betadisp.sim.low <- betadisper(sim.bray.low.mean2, rep(1, length(dataSummer_Low$plot)),
                                 type = "median", bias.adjust=FALSE) 
  betadisp.sim.low
  length(betadisp.sim.low$distance)
  plot(betadisp.sim.low$vectors[ , 1], betadisp.sim.low$vectors[ , 2])
  points(betadisp.sim.low$centroids[1], betadisp.sim.low$centroids[2], 
         col = "red", pch = 19)

# Medium productivity-----------------------------------------------------------
  betadisp.sim.med <- betadisper(sim.bray.med.mean2, rep(1, length(dataSummer_Med$plot)),
                                 type = "median", bias.adjust=FALSE) 
  betadisp.sim.med
  length(betadisp.sim.med$distance)
  plot(betadisp.sim.med$vectors[ , 1], betadisp.sim.med$vectors[ , 2])
  points(betadisp.sim.med$centroids[1], betadisp.sim.med$centroids[2], 
         col = "red", pch = 19) 

# High productivity-------------------------------------------------------------
  betadisp.sim.high <- betadisper(sim.bray.high.mean2, rep(1, length(dataSummer_High$plot)),
                                  type = "median", bias.adjust=FALSE) 
  betadisp.sim.high
  length(betadisp.sim.high$distance)
  plot(betadisp.sim.high$vectors[ , 1], betadisp.sim.high$vectors[ , 2])
  points(betadisp.sim.high$centroids[1], betadisp.sim.high$centroids[2], 
         col = "red", pch = 19) 


#-------------------------------------------------------------------------------
# 10. Calculate effect sizes (es) for each pair of quadrats as the difference 
#      between observed beta-diversity and the mean simulated beta-diversity.
#------------------------------------------------------------------------------- 

# Low productivity plots
  es.low <- bcLow - tapply(sim.bray.low, rep(1:length(bcLow), times = k),mean)
  plot(bcLow, es.low, xlab = "Observed beta-diversity", ylab = "Effect size", 
       cex.lab = 1.5, cex.axis = 1.5)
  abline(h = 0, lty = 3)
  abline(0, 1, col = "red")

# Medium productivity plots
  es.med <- bcMed - tapply(sim.bray.med, rep(1:length(bcMed), times = k), mean)
  plot(bcMed, es.med, xlab = "Observed beta-diversity", ylab = "Effect size", 
       cex.lab = 1.5, cex.axis = 1.5)
  abline(h = 0, lty = 3)
  abline(0, 1, col = "red")

# High productivity plots
  es.high <- bcHigh - tapply(sim.bray.high, rep(1:length(bcHigh), times = k), mean)
  plot(bcHigh, es.high, xlab = "Observed beta-diversity", ylab = "Effect size", 
       cex.lab = 1.5, cex.axis = 1.5)
  abline(h = 0, lty = 3)
  abline(0, 1, col = "red")


#-------------------------------------------------------------------------------
# 11. Calculate standardized effect sizes (ses) for each pair of quadrats as the 
#     effect size divided by the the standard deviation (sd) of the null 
#     distribution of simulated beta-diversity for each pair of quadrats. Save 
#     ses values to text files
#-------------------------------------------------------------------------------

# Low productivity plots
  ses.low <- es.low / tapply(sim.bray.low, rep(1:length(bcLow), times = k), sd)
  plot(es.low, ses.low, xlab = "Effect size", ylab = "Standardized effect size", 
       cex.lab = 1.5, cex.axis = 1.5)
  abline(0, 1, col = "red")
#write.table(as.vector(ses.low), "ses_bray_low.txt", quote = T, sep = ",", row.names = F)

# Medium productivity plots
  ses.med <- es.med / tapply(sim.bray.med, rep(1:length(bcMed), times = k), sd)
  plot(es.med, ses.med, xlab = "Effect size", ylab = "Standardized effect size", 
       cex.lab = 1.5, cex.axis = 1.5)
  abline(0, 1, col = "red")
#write.table(as.vector(ses.med), "ses_bray_med.txt", quote = T, sep = ",", row.names = F) 

# High productivity plots
  ses.high <- es.high / tapply(sim.bray.high, rep(1:length(bcHigh), times = k), sd)
  plot(es.high, ses.high, xlab = "Effect size", ylab = "Standardized effect size", 
       cex.lab = 1.5, cex.axis = 1.5)
  abline(0, 1, col = "red")
#write.table(as.vector(ses.high), "ses_bray_high.txt", quote = T, sep = ",", row.names = F)  


#-------------------------------------------------------------------------------
# 12. Check for negative standardized effect sizes. If there are are negative 
#     effect sizes, add a constant (the smallest value) to make the response 
#     non-negative.
#-------------------------------------------------------------------------------

# Check for negative standardized effect sizes
  summary(ses.low)   # negative effect sizes
  summary(ses.med)   # negative effect sizes
  summary(ses.high)  # negative effect sizes

# Add a constant to make standardized effect sizes positive (these values will
# be used for the beta-dispersion analysis which requires non-negetive values)
  ses.low.positive <- ses.low + (-min(ses.low))
  summary(ses.low.positive)
  length(ses.low.positive)

  ses.med.positive <- ses.med + (-min(ses.med))
  summary(ses.med.positive)
  length(ses.med.positive)

  ses.high.positive <- ses.high + (-min(ses.high))
  summary(ses.high.positive)
  length(ses.high.positive)


#-------------------------------------------------------------------------------
# 13. Calculate distance-to-medians for mean simulated values
#-------------------------------------------------------------------------------

# Low productivity plots
  betadisp.ses.low <- betadisper(ses.low.positive, rep(1, length(dataSummer_Low$plot)),
                                 type = "median", bias.adjust=FALSE) 
  length(betadisp.ses.low$distance)
  plot(betadisp.ses.low$vectors[,1], betadisp.ses.low$vectors[,2])
  points(betadisp.ses.low$centroids[1], betadisp.ses.low$centroids[2], 
         col = "red", pch = 19)

# Medium productivity plots
  betadisp.ses.med <- betadisper(ses.med.positive, rep(1, length(dataSummer_Med$plot)),
                                 type = "median", bias.adjust=FALSE) 
  length(betadisp.ses.med$distance)
  plot(betadisp.ses.med$vectors[,1], betadisp.ses.med$vectors[,2])
  points(betadisp.ses.med$centroids[1], betadisp.ses.med$centroids[2], 
         col = "red", pch = 19)  

# High productivity plots
  betadisp.ses.high <- betadisper(ses.high.positive, rep(1, length(dataSummer_High$plot)),
                                  type = "median", bias.adjust=FALSE) 
  length(betadisp.ses.high$distance)
  plot(betadisp.ses.high$vectors[,1], betadisp.ses.high$vectors[,2])
  points(betadisp.ses.high$centroids[1], betadisp.ses.high$centroids[2], 
         col = "red", pch = 19)   


#-------------------------------------------------------------------------------
# 14. Plot observed beta-diversity, simulated beta-diversity, standardized effect 
#     sizes, and distance-to-spatial medians
#-------------------------------------------------------------------------------

  par(mfrow=c(2,3), las=1, pty="s")

  boxplot(bcLow, bcMed, bcHigh, ylim = c(0, 1), main = "Observed", 
          ylab = "Bray-Curtis dissimilarity", names = c("Low", "Medium", "High"), 
          xlab = "Productivity", col = c("red", "yellow", "green"))

  boxplot(sim.bray.low.mean, sim.bray.med.mean, sim.bray.high.mean, 
          ylim = c(0, 1), main = "Simulated", ylab = "Bray-Curtis dissimilarity", 
          names = c("Low", "Medium", "High"), col = c("red", "yellow", "green"),
          xlab = "Productivity")

  boxplot(ses.low, ses.med, ses.high, ylim = c(-10, 20), 
          main = "Standardized effect size", ylab = "SES", 
          names = c("Low", "Medium", "High"), col = c("red", "yellow", "green"),
          xlab = "Productivity")
  abline(h = 0, lty = 2)

  par(mfrow = c(1,3))
  boxplot(betadisp.bcLow$distance, betadisp.bcMed$distance, betadisp.bcHigh$distance,
          ylim = c(0, 1), 
          ylab = "Distance-to-median", main = "Observed beta-diversity",
          names = c("Low", "Medium", "High"), col = c("red", "yellow", "green"),
          xlab = "Productivity")

  boxplot(betadisp.sim.low$distance, betadisp.sim.med$distance, 
          betadisp.sim.high$distance, main = "Simulated beta-diversity", ylim = c(0, 1),  
          ylab = "Distance-to-median", names = c("Low","Medium", "High"), 
          xlab = "Productiviy", col = c("red", "yellow", "green"))
  
  boxplot(betadisp.ses.low$distance, betadisp.ses.med$distance, 
          betadisp.ses.high$distance, main = "SES beta-dversity", ylim = c(0, 15), 
          ylab = "Distance-to-medium", 
          names = c("Low", "Medium", "High"), col = c("red", "yellow", "green"),
          xlab = "Productivity")


#-------------------------------------------------------------------------------
# 15. Combine distance-to-centroids across treatments for ANOVA $ permANOVA
#-------------------------------------------------------------------------------

  dist.centroid.ses.low.med <- as.vector(c(betadisp.ses.low$distance, 
                                           betadisp.ses.med$distance))
  dist.centroid.ses.low.high <- as.vector(c(betadisp.ses.low$distance, 
                                            betadisp.ses.high$distance))
  dist.centroid.ses.med.high <- as.vector(c(betadisp.ses.med$distance, 
                                            betadisp.ses.high$distance))

  prod.low.med <- c(rep("low", length(betadisp.bcLow$distance)), 
                    rep("med", length(betadisp.bcMed$distance)))
  prod.low.high <- c(rep("low", length(betadisp.bcLow$distance)), 
                     rep("high", length(betadisp.bcHigh$distance)))
  prod.med.high <- c(rep("med", length(betadisp.bcMed$distance)), 
                     rep("high", length(betadisp.bcHigh$distance)))


#-------------------------------------------------------------------------------
# 16. Classic ANOVA on distance-to-centroids
#-------------------------------------------------------------------------------

# low vs med productivity 
  ANOVA.ses.low.med <- aov(dist.centroid.ses.low.med ~ prod.low.med)
  anova(ANOVA.ses.low.med) ## P < 0.1231

## Calculate observed F-value
  obs.F.ses.low.med <- anova(ANOVA.ses.low.med)$"F value"[1]
  obs.F.ses.low.med ## F = 2.447724

# low vs high productivity 
  ANOVA.ses.low.high <- aov(dist.centroid.ses.low.high ~ prod.low.high)
  anova(ANOVA.ses.low.high) ## P < 0.0001

## Calculate observed F-value
  obs.F.ses.low.high <- anova(ANOVA.ses.low.high)$"F value"[1]
  obs.F.ses.low.high ## F = 40.356

# med vs high productivity 
  ANOVA.ses.med.high <- aov(dist.centroid.ses.med.high ~ prod.med.high)
  anova(ANOVA.ses.med.high) ## P < 0.0001

## Calculate observed F-value
  obs.F.ses.med.high <- anova(ANOVA.ses.med.high)$"F value"[1]
  obs.F.ses.med.high ## F = 63.00007

#-------------------------------------------------------------------------------
# 17. Permutational ANOVA on distance-to-centroids
#-------------------------------------------------------------------------------

# Set number of randomizations
  k <- 9999

# Create empty vectors to store the values of the F statistic
  null.F.ses.low.med <- c()
  null.F.ses.low.high <- c()
  null.F.ses.med.high <- c()

# low vs med productivity-------------------------------------------------------

# Run the loop for the permutations
  for (i in 1:k){
    a.low.med <- as.factor(sample(prod.low.med, length(prod.low.med))) 
    permANOVA.ses.low.med <- aov(dist.centroid.ses.low.med ~ a.low.med) 
    null.F.ses.low.med <- append(null.F.ses.low.med, 
                                 anova(permANOVA.ses.low.med)$"F value"[1], 
                                 after = length(null.F.ses.low.med)) 
  }

# Plot results
  hist(c(null.F.ses.low.med, obs.F.ses.low.med), 
       main = "F permutation test, 95 percentile = dotted line, seserved = red line", 
       xlab="F", cex.axis=1.5, cex.lab=1.5)
  abline(v = obs.F.ses.low.med, col = "red") 
  abline(v = quantile(null.F.ses.low.med, prob = 0.95), lty = 3)

# Calculate P-value
  sum(c(null.F.ses.low.med, obs.F.ses.low.med) >= obs.F.ses.low.med)/length(c(null.F.ses.low.med, obs.F.ses.low.med))
  ## P = 0.1234

# low vs high productivity-------------------------------------------------------

# Run the loop for the permutations
  for (i in 1:k){
    a.low.high <- as.factor(sample(prod.low.high, length(prod.low.high))) 
    permANOVA.ses.low.high <- aov(dist.centroid.ses.low.high ~ a.low.high) 
    null.F.ses.low.high <- append(null.F.ses.low.high, 
                                  anova(permANOVA.ses.low.high)$"F value"[1], 
                                  after = length(null.F.ses.low.high)) 
  }

# Plot results
  hist(c(null.F.ses.low.high, obs.F.ses.low.high), 
       main = "F permutation test, 95 percentile = dotted line, seserved = red line", 
       xlab="F", cex.axis=1.5, cex.lab=1.5)
  abline(v = obs.F.ses.low.high, col = "red") 
  abline(v = quantile(null.F.ses.low.high, prob = 0.95), lty = 3)

# Calculate P-value
  sum(c(null.F.ses.low.high, obs.F.ses.low.high) >= obs.F.ses.low.high)/length(c(null.F.ses.low.high, obs.F.ses.low.high))
  ## P < 0.0001


# med vs high productivity-------------------------------------------------------

# Run the loop for the permutations
  for (i in 1:k){
    a.med.high <- as.factor(sample(prod.med.high, length(prod.med.high))) 
    permANOVA.ses.med.high <- aov(dist.centroid.ses.med.high ~ a.med.high) 
    null.F.ses.med.high <- append(null.F.ses.med.high, 
                                  anova(permANOVA.ses.med.high)$"F value"[1], 
                                  after = length(null.F.ses.med.high)) 
  }

# Plot results
  hist(c(null.F.ses.med.high, obs.F.ses.med.high), 
       main = "F permutation test, 95 percentile = dotted line, seserved = red line", 
       xlab="F", cex.axis=1.5, cex.lab=1.5)
  abline(v = obs.F.ses.med.high, col = "red") 
  abline(v = quantile(null.F.ses.med.high, prob = 0.95), lty = 3)

# Calculate P-value
  sum(c(null.F.ses.med.high, obs.F.ses.med.high) >= obs.F.ses.med.high)/length(c(null.F.ses.med.high, obs.F.ses.med.high))
  ## P < 0.0001


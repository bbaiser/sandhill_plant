library(vegan)
library(lme4)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(DHARMa)

source("http://highstat.com/Books/BGS/GAMM/RCodeP2/HighstatLibV6.R")
################################################
# load data
################################################

# Species data in species-by-site matrix

plants <- t(read.csv("data/FocalPlants.csv",row=1)) #transpose so it is site-by-species

#environmental variables
Env<-read.csv("data/Env_variables.csv",row=1 )

#z-scale predictor variables:
Env_Z<-scale(Env[,2:15])
Env_Z2<-data.frame(Env_Z,Env$Site)#join site column back
colnames(Env_Z2)[15]<-"Site"#rename column
corvif(Env_Z)
#Join Env_Z and plants
plant_env<-merge(Env_Z2, plants, by=0, all=TRUE) 

          
#make into long format 
final_dat = gather(plant_env, key = "SPP", value = "COVER", BALANG:STYPAT)

colnames(final_dat)[1]<- 'Plot' #change column name for plots to "Plot"
################################################
# end load data
################################################


#AIC Model Selection (try all possible models)

#make polynomial terms...
final_dat$soil_moisture2 <- final_dat$soil_moisture ^ 2
final_dat$soil_moisture2 <- final_dat$soil_moisture ^ 2
final_dat$soil_moisture2 <- final_dat$soil_moisture ^ 2
final_dat$soil_moisture2 <- final_dat$soil_moisture ^ 2
final_dat$soil_moisture2 <- final_dat$soil_moisture ^ 2

#Full model- no check for VIF
full_mod<-lmer(COVER~soil_moisture+(0 + soil_moisture | SPP)
          + herb_cover+ (0 + herb_cover | SPP) 
          + shrub_cover+ (0 + shrub_cover | SPP)
          + LAI_tree_layer +(0+ LAI_tree_layer|SPP)
          + LAI_shrub_layer + (0+LAI_shrub_layer|SPP)
          + LAI_total + (0+ LAI_total|SPP)
          + Total_Fires + (0+ Total_Fires)
          + Fire_Simpson_Div + (0+Fire_Simpson_Div|SPP)
          + fire_return + (0+fire_return|SPP)
          + yrs_in_rotation+ (0+yrs_in_rotation|SPP)
          + spring_sum_fire + (0+spring_sum_fire|SPP)
          + fall_wint_fire + (0+fall_wint_fire|SPP)
          + (1|SPP)
          +(1|Site),
          data = final_dat)



#NO random SLOPE and vars removed based on VIF

full<-lmer(COVER~soil_moisture
               + herb_cover
               + shrub_cover
               + LAI_tree_layer 
               + LAI_shrub_layer 
               + Total_Fires
               + Fire_Simpson_Div 
               + yrs_in_rotation+
               +(1|SPP)
               +(1|Site),
               data = final_dat)
car::vif(full)
summary(full)
sjstats::r2(full) 

#with random slopes
full2<-lmer(COVER~soil_moisture+(0 + soil_moisture | SPP) 
           + herb_cover
           + shrub_cover+(0 + shrub_cover | SPP) 
           + LAI_tree_layer 
           + LAI_shrub_layer 
           + Total_Fires
           + Fire_Simpson_Div +(0 + Fire_Simpson_Div  | SPP) 
           + yrs_in_rotation+(0 + yrs_in_rotation  | SPP) 
           +(1|SPP)
           +(1|Site),
           ziformula = ~1
           data = final_dat)
car::vif(full2)
summary(full2)
performance::r2(full2) 
drop1(full2)


#explore residuals with dHarma
simulationOutput <- simulateResiduals(fittedModel = full2, plot = F)
res<-residuals(simulationOutput)
hist(res, breaks=100)
plot(simulationOutput)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)
hist(final_dat$COVER, breaks=100)

#test fewervaraibles after error message
test_mod<-lmer(COVER~soil_moisture+(0 + soil_moisture | SPP)
               +  (1|SPP)
               +(1|Site),
               data = final_dat)
str(Env)
summary(test_mod)
sjstats::r2(full)             
 + herb_cover+ (0 + herb_cover | SPP) 
               + shrub_cover+ (0 + shrub_cover | SPP)
               + LAI_tree_layer +(0+ LAI_tree_layer|SPP)
               + LAI_shrub_layer + (0+LAI_shrub_layer|SPP)
               + LAI_total + (0+ LAI_total|SPP)
               + Total_Fires + (0+ Total_Fires)
               + Fire_Simpson_Div + (0+Fire_Simpson_Div|SPP)
               + fire_return + (0+fire_return|SPP)
               + yrs_in_rotation+ (0+yrs_in_rotation|SPP)
               + spring_sum_fire + (0+spring_sum_fire|SPP)
               + fall_wint_fire + (0+fall_wint_fire|SPP)
               + (1|SPP)
               +(1|Site),
               data = final_dat)
best <- glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + (0 + Elevation | SPP) +
                Herb + Herb2 + (0 + Herb | SPP) + TreeBA + ACSA + LITU + (0 + LITU | SPP) +
                Ca + Ca2 + (0 + Ca | SPP) + P + (0 + P | SPP),
              family = binomial,
              data = hh)

summary(best)

ranef(best)$SPP$Elevation + 0.825
fixef(best)

nsite <- dim(h)[1]
nx <- 10 #number of fixed effects
nspp <- 14 #number of species

################################################
# analyze data: lmer (MLM)

z = glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + (0 + Elevation | SPP) +
            Herb + Herb2 + (0 + Herb | SPP) + TreeBA + ACSA + LITU + (0 + LITU | SPP) +
            Ca + Ca2 + (0 + Ca | SPP) + P + (0 + P | SPP),
          family = binomial,
          data = hh
)

# compute ranef pvalues
#Without elevation
z1 = glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 +
             Herb + Herb2 + (0 + Herb | SPP) + TreeBA + ACSA + LITU + (0 + LITU | SPP) +
             Ca + Ca2 + (0 + Ca | SPP) + P + (0 + P | SPP),
           family = binomial,
           data = hh
)

#Without Herb
z2 = glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + (0 + Elevation | SPP) +
             Herb + Herb2 + TreeBA + ACSA + LITU + (0 + LITU | SPP) +
             Ca + Ca2 + (0 + Ca | SPP) + P + (0 + P | SPP),
           family = binomial,
           data = hh
)

#Without LITU
z3 = glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + (0 + Elevation | SPP) +
             Herb + Herb2 + (0 + Herb | SPP) + TreeBA + ACSA + LITU +
             Ca + Ca2 + (0 + Ca | SPP) + P + (0 + P | SPP),
           family = binomial,
           data = hh
)

#Without Ca
z4 = glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + (0 + Elevation | SPP) +
             Herb + Herb2 + (0 + Herb | SPP) + TreeBA + ACSA + LITU + (0 + LITU | SPP) +
             Ca + Ca2 +  P + (0 + P | SPP),
           family = binomial,
           data = hh
)

#Without P
z5 = glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + (0 + Elevation | SPP) +
             Herb + Herb2 + (0 + Herb | SPP) + TreeBA + ACSA + LITU + (0 + LITU | SPP) +
             Ca + Ca2 + (0 + Ca | SPP) + P,
           family = binomial,
           data = hh
)

LL <- c(deviance(z), deviance(z1), deviance(z2), deviance(z3), deviance(z4),
        deviance(z5))

#these are chi-square distribution functions). 
mlm.pvals <- c(
  1 - pchisq(LL[2] - LL[1], 1),
  1 - pchisq(LL[3] - LL[1], 1),
  1 - pchisq(LL[4] - LL[1], 1),
  1 - pchisq(LL[5] - LL[1], 1),
  1 - pchisq(LL[6] - LL[1], 1)
)

mlm.pvals

# compute residual effect of random effects environmental variables
z.r <- glmer(PRESENCE ~ (1 | SPP) + Elevation + Elevation2 + Herb + Herb2 +
               TreeBA + ACSA + LITU + Ca + Ca2 + P,
             data = hh,
             family = binomial
)

MLM.fitted <- array(fitted(z) - fitted(z.r), c(nsite, nspp))

rownames(MLM.fitted) = c(1:54)
colnames(MLM.fitted) = names(h)[26:39]

# standardize over spp
MLM.fitted.standard <- MLM.fitted
for (j in 1:nspp)
  MLM.fitted.standard[, j] <-
  (MLM.fitted[, j] - mean(MLM.fitted[, j])) /
  sd(MLM.fitted[, j])
# another way to do this??

# browseURL("http://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca")
ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[, 1:2]

# environmental variables (only those with random effects)
envir.vars <- cbind(hh$Elevation, hh$Herb, hh$LITU, hh$Ca, hh$P)
mlm.envir <- NULL
for (j in 1:5)
  mlm.envir <-
  cbind(mlm.envir, envir.vars[, j] * mlm.fit[, 1], envir.vars[, j] * mlm.fit[, 2])

envir.points <- matrix(colMeans(mlm.envir), nrow = ncol(envir.vars), ncol = 2, byrow = T)

# plot mlm
par(mfcol = c(1, 1))
plot(-mlm.fit,
     xlab = "PC1",
     ylab = "PC2",
     type = "n")
text(-mlm.fit, label = 1:nsite, cex = .5)

arrow.coordMLM <- cbind(array(0, dim(envir.points)), -envir.points)

arrows(
  arrow.coordMLM[, 1],
  arrow.coordMLM[, 2],
  arrow.coordMLM[, 3],
  arrow.coordMLM[, 4],
  col = "black",
  length = 0.1
)

text(
  1.3 * -envir.points,
  label = c("Elevation", "Herb", "LITU", "Ca", "P"),
  cex = .7
)


################################################
# analyze data: CCA
comm.matrix = h[, 26:39]

envir.matrix <- h[c("Elevation", "Herb", "TreeBA", "ACSA", "LITU", "Ca", "P")]
cbind(hh$Elevation, hh$Herb, hh$TreeBA, hh$ACSA, hh$LITU, hh$Ca, hh$P)[1:nsite, ]
colnames(envir.matrix) = c("Elevation", "Herb", "TreeBA", "ACSA", "LITU", "Ca", "P")

cca.fit <- cca(comm.matrix, envir.matrix, scale = TRUE)

nx.cca=7
cca.partial.pvals <- NULL
for (j in 1:nx.cca) {
  # partial cca
  Z <- envir.matrix[, -j]
  cca.fit.partial <- cca(comm.matrix, envir.matrix[, j], Z, scale = TRUE)
  anova.cca <- permutest(cca.fit.partial, permutations = 10000)
  pval <- mean(anova.cca$F.0 < anova.cca$F.perm)
  cca.partial.pvals <- c(cca.partial.pvals, pval)
}
cca.partial.pvals

################################################
# analyze data: RDA
rda.fit <- rda(comm.matrix, envir.matrix, scale = TRUE)

rda.partial.pvals <- NULL
for (j in 1:nx.cca) {
  # partial rda
  Z <- envir.matrix[, -j]
  rda.fit.partial <- rda(comm.matrix, envir.matrix[, j], Z, scale = TRUE)
  anova.rda <- permutest(rda.fit.partial, permutations = 10000)
  pval <- mean(anova.rda$F.0 < anova.rda$F.perm)
  rda.partial.pvals <- c(rda.partial.pvals, pval)
}
rda.partial.pvals

################################################
# analyze data: NMDS
nmds.fit <-
  metaMDS(
    comm.matrix,
    distance = "bray",
    k = 2,
    trymin = 200,
    trymax = 300,
    zerodist = "add",
    trace = 0
  )
nmds.fit$stress

nmds.envir <- envfit(nmds.fit, envir.matrix, permutations = 10000)
nmds.pvals <- nmds.envir$vectors[4]$pvals
nmds.pvals

################################################
# analyze data: lm
z.lm <- lm(PRESENCE ~ SPP * Elevation + SPP * Herb + SPP * TreeBA + SPP * ACSA + SPP *
             LITU + SPP * Ca + SPP * P, data = hh)

MLM.fitted.lm <- array(z.lm$fitted.values, c(nsite, nspp))
colnames(MLM.fitted.lm) = colnames(comm.matrix)
rownames(MLM.fitted.lm) = c(1:54)
mlm.fit.lm <- rda(MLM.fitted.lm, scale = TRUE)
# mlm.envir.lm
mlm.envir.lm = rda(MLM.fitted.lm, Y = envir.matrix, scale = T)

################################################
# plot results
par(mfrow = c(2, 2))

# plot mlm
plot(-mlm.fit,
     xlab = "PC1",
     ylab = "PC2",
     type = "n")
text(-mlm.fit, label = 1:nsite, cex = .5)
arrows(
  arrow.coordMLM[, 1],
  arrow.coordMLM[, 2],
  arrow.coordMLM[, 3],
  arrow.coordMLM[, 4],
  col = "black",
  length = 0.05
)

text(
  1.3 * -envir.points,
  label = c("Elevation", "Herb", "LITU", "Ca", "P"),
  cex = .7
)

# plot rda
plot(rda.fit, type = "n")
text(rda.fit,
     display = "lc",
     cex = 0.5,
     col = "black")
text(rda.fit,
     display = "bp",
     cex = 0.7,
     select = c(1, 4:7))

# plot cca
plot(cca.fit, type = "n", xlim = c(-8, 5))
text(cca.fit,
     display = "lc",
     cex = 0.5,
     col = "black")
text(
  cca.fit,
  display = "bp",
  cex = 0.7,
  select = c(1, 2, 5, 6, 7)
)

# plot nmds
plot(nmds.fit, type = "n")
text(nmds.fit,
     display = c("sites"),
     col = "black",
     cex = 0.5)

arrow.list = which(nmds.pvals < 0.05) # choose significant vectors
arrow.rvalues = sqrt(nmds.envir$vector$r[arrow.list])
# correlations with species determine length of the arrow
arrow.endpt = nmds.envir$vectors$arrows[arrow.list, ] * arrow.rvalues
arrow.endpt = 0.7 * arrow.endpt / max(arrow.endpt[1, ])
# adjust length of arrows to fit plot
# 0.7 is a scaling coefficient; change if desired
arrow.coord = cbind(0, 0, arrow.endpt)
arrow.labels = row.names(arrow.coord)
arrows(
  arrow.coord[, 1],
  arrow.coord[, 2],
  arrow.coord[, 3],
  arrow.coord[, 4],
  col = "black",
  length = 0.05
)
text(1.1 * arrow.coord[, 3:4], label = arrow.labels, cex = 0.7)

###############################################
dev.set(4)
par(mfrow = c(1, 2))

plot(mlm.fit.lm, type = "n")
text(mlm.fit.lm,
     display = "sites",
     cex = 0.5,
     col = "black")
text(mlm.envir.lm, display = "bp", cex = 0.7)

plot(rda.fit, type = "n")
text(rda.fit, cex = 0.5, col = "black")
text(rda.fit, display = "bp", cex = 0.7)

################################################
# method comparisons
################################################

# access site scores
mlm.scores <- mlm.fit
rda.scores <- scores(rda.fit, display = "lc")
cca.scores <- scores(cca.fit, display = "lc")
nmds.scores <- nmds.fit$points

# procrustes
procrustes(mlm.scores,
           rda.scores,
           scale = TRUE,
           symmetric = TRUE) #0.1857
procrustes(mlm.scores,
           cca.scores,
           scale = TRUE,
           symmetric = TRUE) #0.1872
procrustes(mlm.scores,
           nmds.scores,
           scale = TRUE,
           symmetric = TRUE) #0.63

procrustes(rda.scores,
           cca.scores,
           scale = TRUE,
           symmetric = TRUE) #0.1373
procrustes(rda.scores,
           nmds.scores,
           scale = TRUE,
           symmetric = TRUE) #0.7258
procrustes(cca.scores,
           nmds.scores,
           scale = TRUE,
           symmetric = TRUE) #0.7075


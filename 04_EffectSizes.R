## Script to calculate effect sizes

library(Hmisc)
library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished


# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")

## LOAD THE DATA
dataDir <- "Data/February2022"

DataCleaned <- read.csv(file.path(dataDir, "processed", "0_2_alldata.csv"), header = TRUE) # load the cleaned data (Change the name)

table(DataCleaned$Measurement)


# Sometimes zeros were missing that should actually be zeros
DataCleaned$Control_SD[which(DataCleaned$Control_mean == 0 & is.na(DataCleaned$Control_SD))] <- 0
DataCleaned$Treatment_SD[which(DataCleaned$Treatment_mean == 0 & is.na(DataCleaned$Treatment_SD))] <- 0

## There's studies that rounded their SD's so that they were zero
DataCleaned$Control_SD[which(DataCleaned$ID == 757 & DataCleaned$Control_SD == 0)] <- 0.01
DataCleaned$Treatment_SD[which(DataCleaned$ID == 757 & DataCleaned$Treatment_SD == 0)] <- 0.01

DataCleaned$Control_SD[which(DataCleaned$ID == 1857 & DataCleaned$Control_SD == 0)] <- 0.01
DataCleaned$Treatment_SD[which(DataCleaned$ID == 1857 & DataCleaned$Treatment_SD == 0)] <- 0.01


## One case has an NA for SD, when it is actually a zero (according to the paper)
DataCleaned$Treatment_SD[which(DataCleaned$ID == 129 & is.na(DataCleaned$Treatment_SD))] <- 0



## Which measurements are there that have nas in Sds
table(DataCleaned$Measurement[which(!(is.na(DataCleaned$Control_mean)) & is.na(DataCleaned$Control_SD))])
## That's fine (11 abundance and 1 biomass)



# there are some climate cases that have no SDs
DataCleaned <- DataCleaned[which(!(is.na(DataCleaned$Control_SD))),] 
DataCleaned <- DataCleaned[which(!(is.na(DataCleaned$Treatment_SD))),] # 3440



## How many cases (and which drivers) now have zeros
nrow(DataCleaned[which(DataCleaned$Control_mean == 0),]) # 127
nrow(DataCleaned[which(DataCleaned$Control_SD   == 0),]) # 159

table(DataCleaned$driver[which(DataCleaned$Control_SD   == 0)]) # Mainly lui, nutrient enrichment and pollution


nrow(DataCleaned[which(DataCleaned$Treatment_mean  == 0),]) # 161
nrow(DataCleaned[which(DataCleaned$Treatment_SD  == 0),]) # 197
table(DataCleaned$driver[which(DataCleaned$Treatment_SD   == 0)]) # Mainly lui,  nutrient enrichment and pollution




hedges <- escalc(measure = "SMD", # log response ratio ("ROM" in metafor)
                      m2i = Control_mean, # group 2 corresponds to the control group
                      sd2i = Control_SD,
                      n2i = Control_N,
                      m1i = Treatment_mean, # group 1 is the treatment group
                      sd1i = Treatment_SD,
                      n1i = Treatment_N,
                      data = DataCleaned)
## get a warning message here because of the zeros, which gets dealt with shortly
# 3440

## 99 NAs in hedges
# When theres no Ns, or when there are zeros, in the Sds?


hedges <- hedges[which(!(is.na(hedges$yi))),]


# Rename Effect sizes
hedges$effect <- hedges$yi
hedges$var <- hedges$vi

hist(hedges$effect)
summary(hedges$effect)

hist(hedges$var)
summary(hedges$var)

## What have the high effect sizes and var

head(hedges[which(hedges$effect  > 75),]) # one case, and is legit



head(hedges[which(hedges$var > 2),])
table(hedges$Measurement[which(hedges$var > 10)])
head(hedges[which(hedges$var > 10),])


# TODO: check the high vars with Lea



# Save data
write.csv(hedges, "Data/03_Data/HedgesData.csv")
# 3341


hedges$UniqueID <- paste(hedges$ID, hedges$Case_ID, hedges$driver)


## TRYING TO ESTABLISH WHETHER ALL EFFECT SIZES CAN BE USED ----------


measurement.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

qqnorm(residuals(measurement.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(measurement.mod.1,type="pearson"),col="red")
# not bad


## Everything is significantly different.
## So split everything via measurement type



## Make biomass the baseline
# Maybe biomass can be with richness??
biomass_baseline <- hedges
biomass_baseline$Measurement <- as.factor(biomass_baseline$Measurement)
biomass_baseline$Measurement <- relevel(biomass_baseline$Measurement, ref = "Biomass")

biomass.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=biomass_baseline)
## As I suspected, Richness (and Shannon) are not significantly different from Biomass



## Let's look at bodysize

bodysize.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Body.Size, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)

summary(bodysize.mod.1)


## FULL MODEL ---------




# If fittiung with REML (which is recommended) then can't compare models 
# (https://stats.stackexchange.com/questions/48671/what-is-restricted-maximum-likelihood-and-when-should-it-be-used)
# https://stats.stackexchange.com/questions/517857/testing-the-effect-of-moderators-in-metafor-package
# so test moderators using the anova(singleMod) way




# Remove indices with little data
hedges <- hedges[!(hedges$Body.Size == ""),]
hedges <- hedges[!(hedges$Measurement == "FunctionalRichness"),]
hedges <- hedges[!(hedges$Measurement == "GroupAbundance"),]
hedges <- hedges[!(hedges$Measurement == "Simpson"),]
hedges <- hedges[!(hedges$Measurement == "Evenness"),]



hedges$Body.Size[which(hedges$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
hedges$Body.Size[which(hedges$Body.Size %in% c("Macro-arthropods", "Macro-invertebrates"))] <- "Macro-fauna"
hedges$Body.Size[which(hedges$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"


# Make macrofauna the baseline

hedges$Body.Size <- as.factor(hedges$Body.Size)
hedges$Body.Size <- relevel(hedges$Body.Size, ref = "Macro-fauna")


write.csv(hedges, "Data/03_Data/HedgesData_cleaned.csv")



## quick test for variance covariance structure

mod.un<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * Body.Size + Measurement, ## 
  random= ~ as.factor(ID)|UniqueID, # need to understand better about these inner and outer terms
  struct="UN",
  method="REML",
  digits=4,
  data=hedges)

# beep()


# Testing crossed random effects
mod.crossed<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()
anova(mod.crossed, btt = ":") # not significant

mod.crossed2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)
beep()
anova(mod.crossed2, btt = "Body") # not significant
anova(mod.crossed2, btt = "driver") #  significant


mod.crossed_sigma0<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges,
  sigma2 =  c(NA, NA,0))
anova(mod.crossed, mod.crossed_sigma0)
beep()

anova(mod.1, mod.crossed_sigma0)







mod.1a<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)


mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

mod.1







mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

anova(mod.1, btt = ":") # should be testing the levels that are the interactions
## not significant
anova(mod.1, btt = "Measurement") #  significant

# Ultimately, this interaction is not needed


mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)
beep()


anova(mod.1b, btt = "driver") # significant
anova(mod.1b, btt = "Body") # not significant
anova(mod.1b, btt = "Measurement") #  significant


qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
# not bad


mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

anova(mod.2, btt = "driver") # significant
anova(mod.2, btt = "Measurement") #  significant




## we stick with mod.2
saveRDS(mod.2, file = "Models/MainMod.rds")



# https://stats.stackexchange.com/questions/109841/presenting-results-of-a-meta-analysis-with-multiple-moderators
### compute the estimated/predicted correlation for each combination
test <-predict(mod.2, newmods=rbind(c(0,0,0,0,0,0,0,0),
                             c(0,0,0,0,0,1,0,0),
                             c(0,0,0,0,0,0,1,0),
                             c(0,0,0,0,0,0,0,1),  # climate with all different measurement types
                             c(1,0,0,0,0,0,0,0),
                             c(1,0,0,0,0,1,0,0),
                             c(1,0,0,0,0,0,1,0),
                             c(1,0,0,0,0,0,0,1),# habitat 
                             c(0,1,0,0,0,0,0,0),
                             c(0,1,0,0,0,1,0,0),
                             c(0,1,0,0,0,0,1,0),
                             c(0,1,0,0,0,0,0,1),  # invasive
                             c(0,0,1,0,0,0,0,0),
                             c(0,0,1,0,0,1,0,0),
                             c(0,0,1,0,0,0,1,0),
                             c(0,0,1,0,0,0,0,1),  # lui
                             c(0,0,0,1,0,0,0,0),
                             c(0,0,0,1,0,1,0,0),
                             c(0,0,0,1,0,0,1,0),
                             c(0,0,0,1,0,0,0,1), # nutrient
                             c(0,0,0,0,1,0,0,0),
                             c(0,0,0,0,1,1,0,0),
                             c(0,0,0,0,1,0,1,0),
                             c(0,0,0,0,1,0,0,1)# pollution
                             ), addx=TRUE, digits=2) #

slabs <- rep(c("abundance", "biomass", "richness", "shannon"), times = 6)
par(mar=c(3, 8, 1, 1))
forest(test$pred, sei=test$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))
abline(h=4.5, b=0, lty= 2)
abline(h=8.5, b=0, lty= 2)
abline(h=12.5, b=0, lty= 2)
abline(h=16.5, b=0, lty= 2)
abline(h=20.5, b=0, lty= 2)


mtext("Climate", side = 2, line = 3, at = 22.5)
mtext("Fragmentation", side = 2, line = 3, at = 18.5, cex = 0.55)
mtext("Invasives", side = 2, line = 3, at = 14.4, cex = 0.8)
mtext("LUI", side = 2, line = 3, at = 10.7)
mtext("Nutrient", side = 2, line = 3, at = 6.5)
mtext("Pollution", side = 2, line = 3, at = 2)






## SPLITTING DRIVERS --------

climate <- hedges[which(hedges$driver == "Climate"),] # 462

table(climate$GCDType)
table(climate$GCDType, climate$Measurement)


## going to remove vegetation as a gcd type. 2 lines are thickness of moss, and the 10 other lines are canopy gaps due to ice storm
climate <- climate[which(climate$GCDType != "Vegetation"),]

# For now will also remove the precip+temp
climate <- climate[which(climate$GCDType != "Precipitation+Temp"),]

# Checking body size
table(climate$Body.Size)
# removing the all sizes
climate <- droplevels(climate[which(climate$Body.Size != "All sizes"),])


climate.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

summary(climate.mod.1)

anova(climate.mod.1, btt = ":") # should be testing the levels that are the interactions
# Not significant


climate.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

anova(climate.mod.1b, btt = "GCD") # significant
anova(climate.mod.1b, btt = "Measurement") # not significant
anova(climate.mod.1b, btt = "Size") # not significant


climate.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)


## Use climate.mod.5
saveRDS(climate.mod.5, file = "Models/ClimateMod.rds")


qqnorm(residuals(climate.mod.5,type="pearson"),main="QQ plot: residuals")
qqline(residuals(climate.mod.5,type="pearson"),col="red") # fine




climatedat <-predict(climate.mod.5, newmods=rbind(c(0,0,0,0),
                                    c(1,0,0,0),
                                    c(0,1,0,0),
                                    c(0,0,1,0),  
                                    c(0,0,0,1)
), addx=TRUE, digits=2) #

slabs <- c("Gas", "Temperature",
           "UVB Radiation", "Water Availability-Drought", "Water Availability-Flood")
par(mar=c(3, 8, 1, 1))
forest(climatedat$pred, sei=climatedat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))


## LUI -------

lui <- hedges[which(hedges$driver == "LUI"),] # 911

table(lui$GCDType)
lui$GCDType[which(lui$GCDType %in% c("Defoliation"))] <- "Grazing"
lui$GCDType[which(lui$GCDType %in% c("Intensity", "Landscape", "Management", "Weeding", "Planting", "Intensification"))] <- "Management"
lui$GCDType[which(lui$GCDType %in% c("Logging"))] <- "Harvesting"
lui$GCDType[which(lui$GCDType %in% c("Water"))] <- "Irrigation"
lui$GCDType[which(lui$GCDType %in% c("degradation", "Disturbance"))] <- "Degradation"

lui <- lui[which(lui$GCDType != "Human population"),]
lui <- lui[which(lui$GCDType != "Mono- versus poly-culture"),] # too little
lui <- lui[which(lui$GCDType != "Irrigation"),] # too little



table(lui$GCDType, lui$Measurement)



## Checking body size
table(lui$Body.Size)



lui.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=lui)


anova(lui.mod.1, btt = ":") # not significant
anova(lui.mod.1, btt = "Measurement") # significant


summary(lui.mod.1)



lui.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=lui)



anova(lui.mod.1b, btt = "GCD") # significant
anova(lui.mod.1b, btt = "Measurement") # significant
anova(lui.mod.1b, btt = "Size") #  significant at 0.05 level. but we will have level at 0.01

summary(lui.mod.1b)


lui.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=lui)

anova(lui.mod.2, btt = "GCD") # significant
anova(lui.mod.2, btt = "Measurement") # significant
summary(lui.mod.2)



## final model
saveRDS(lui.mod.2, file = "Models/LUIMod.rds")




luidat <-predict(lui.mod.2, newmods=rbind(c(0,0,0,0,0,0,0,0,0),
                                              c(0,0,0,0,0,0,1,0,0),
                                              c(0,0,0,0,0,0,0,1,0),
                                              c(0,0,0,0,0,0,0,0,1),  # intercept
                                              
                                              c(1,0,0,0,0,0,0,0,0),
                                              c(1,0,0,0,0,0,1,0,0),
                                              c(1,0,0,0,0,0,0,1,0),
                                              c(1,0,0,0,0,0,0,0,1),# fire
                                              
                                              c(0,1,0,0,0,0,0,0,0),
                                              c(0,1,0,0,0,0,1,0,0),
                                              c(0,1,0,0,0,0,0,1,0),
                                              c(0,1,0,0,0,0,0,0,1),  # grazing   
                                              
                                              c(0,0,1,0,0,0,0,0,0),
                                              c(0,0,1,0,0,0,1,0,0),
                                              c(0,0,1,0,0,0,0,1,0),
                                              c(0,0,1,0,0,0,0,0,1),  # Harvesting
                                              
                                              
                                              c(0,0,0,1,0,0,0,0,0),
                                              c(0,0,0,1,0,0,1,0,0),
                                              c(0,0,0,1,0,0,0,1,0),
                                              c(0,0,0,1,0,0,0,0,1),# Management
                                              
                                              
                                              c(0,0,0,0,1,0,0,0,0),
                                              c(0,0,0,0,1,0,1,0,0),
                                              c(0,0,0,0,1,0,0,1,0),
                                              c(0,0,0,0,1,0,0,0,1), # organic
                                              
                                              c(0,0,0,0,0,1,0,0,0),
                                              c(0,0,0,0,0,1,1,0,0),
                                              c(0,0,0,0,0,1,0,1,0),
                                              c(0,0,0,0,0,1,0,0,1)# tillage
), addx=TRUE, digits=2) #




slabs <- rep(c("abundance", "biomass", "richness", "shannon"), times = 7)
par(mar=c(3, 8, 1, 1))
forest(luidat$pred, sei=luidat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))
abline(h=4.5, b=0, lty= 2)
abline(h=8.5, b=0, lty= 2)
abline(h=12.5, b=0, lty= 2)
abline(h=16.5, b=0, lty= 2)
abline(h=20.5, b=0, lty= 2)
abline(h=24.5, b=0, lty= 2)


mtext("Degradation", side = 2, line = 3, at = 26.5, cex = 0.5)
mtext("Fire", side = 2, line = 3, at = 22.5, cex = 0.7)
mtext("Grazing", side = 2, line = 3, at = 18.4, cex = 0.7)
mtext("Harvesting", side = 2, line = 3, at = 14.4, cex = 0.5)
mtext("Management", side = 2, line = 3,at = 10.4, cex = 0.5)
mtext("Inorganic", side = 2, line = 3, at = 6.4, cex = 0.6)
mtext("Tillage", side = 2, line = 3, at = 2.4, cex = 0.7)

luidat_abundance <- luidat[c(seq(1, 28, 4)),]

slabs <- c("Degradation", "Fire", "Grazing", "Harvesting", "Management", "Inorganic", "Tillage")
par(mar=c(3, 8, 1, 1))
forest(luidat_abundance$pred, sei=luidat_abundance$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))


#3 Just showing the impact of measurement

measurement_coefs <- data.frame(meas = c("abundance", "biomass", "richness", "shannon"),
                                coef = c(0, lui.mod.2$beta[8:10]), 
                                ses = c(0, lui.mod.2$se[8:10]))


measurement_coefs <- measurement_coefs[c(4,3,2,1),]

errbar(x =measurement_coefs$meas, y = measurement_coefs$coef, 
      yplus = measurement_coefs$coef + measurement_coefs$ses,
       yminus=measurement_coefs$coef-measurement_coefs$ses, cap=0.015)
abline(v=0, lty =2)

## NUTRIENT ENRICHMENT -------



nutri <- hedges[which(hedges$driver == "NutrientEnrichment"),] # 820

table(nutri$GCDType)


table(nutri$GCDType, nutri$Measurement)


# Removing some GCDs and relevelling
nutri <- nutri[which(nutri$GCDType != "Biochar"),]
nutri$GCDType <- as.factor(nutri$GCDType)
nutri$GCDType <- relevel(nutri$GCDType, ref = "Synthetic Fertilizers")


## Check body size
table(nutri$Body.Size)



nutri.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(nutri.mod.1)


anova(nutri.mod.1, btt = ":") # not significant
anova(nutri.mod.1, btt = "Measurement") # not significant


nutri.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

anova(nutri.mod.1b, btt = "GCD") #  significant
anova(nutri.mod.1b, btt = "Measurement") # not significant at 0.01 level
anova(nutri.mod.1b, btt = "Size") # not significant


nutri.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(nutri.mod.2)

saveRDS(nutri.mod.2, file = "Models/nutriMod.rds")



nutridat <-predict(nutri.mod.2, newmods=rbind(c(0,0,0,0,0,0,0),
                               # intercept (synthetic Fertilizers)
                                    
                                    
                                    c(1,0,0,0,0,0,0),
                         #  "Ca-liming + Wood ash"                                         
                                    
                                    c(0,1,0,0,0,0,0),
                                   # Compost
                                    
                                    c(0,0,1,0,0,0,0),
                              # Manure + Slurry 
                                    
                                    c(0,0,0,1,0,0,0),
                               # OMixture
                                    
                                    c(0,0,0,0,1,0,0),
                             # rOther Organic fertilisers
                                    
                                    c(0,0,0,0,0,1,0),
                        # Residue 
                            
                                    c(0,0,0,0,0,0,1)# Sludge 
), addx=TRUE, digits=2) #


slabs <- c("Synthetic Fertilizers", "Ca-liming + Wood ash", "Compost", "Manure + Slurry", "Mixture", 
           "Other Organic fertilisers", "Residue + Mulch", "Sludge")
par(mar=c(3, 8, 1, 1))
forest(nutridat$pred, sei=nutridat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))





## Invasives

invas <- hedges[which(hedges$driver == "Invasives"),] # 188

## GCD type
table(invas$GCDType)
# Probably need to remove plants-mixture
invas <- invas[which(invas$GCDType != "Plants-mixture"),]


# Measurement type

table(invas$GCDType, invas$Measurement)
invas <- invas[which(invas$Measurement == "Abundance"),] # just use abundance

## Body size
table(invas$GCDType, invas$Body.Size)



invas.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=invas)
anova(invas.mod.1, btt = ":") # not  significant


invas.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=invas)



anova(invas.mod.2, btt = "GCD") # not  significant
anova(invas.mod.2, btt = "Size") # not  significant




invas.mod.3<-rma.mv(
  yi=effect,
  V=var, 
 #  mods=~GCDType + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=invas)

summary(invas.mod.3) # not significant intercapt
saveRDS(invas.mod.3, file = "Models/invasiveMod.rds")




### POLLUTION
poll <- hedges[which(hedges$driver == "Pollution"),] # 850

# pollution type
table(poll$GCDType)

poll <- poll[which(poll$GCDType %in% c("Metals", "Pesticides")),] 


# measurement
table(poll$GCDType, poll$Measurement)

# body size
table(poll$Body.Size)
table(poll$GCDType, poll$Body.Size)



poll.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType *Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

anova(poll.mod.1, btt = ":") #  not significant
anova(poll.mod.1, btt = "Measurement") #  not significant


poll.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

anova(poll.mod.1b, btt = "GCD") #  nearly significant
anova(poll.mod.1b, btt = "Size") #  not significant



poll.mod.1c<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=poll)
anova(poll.mod.1c, btt = "GCD") #  nearly significant






poll.mod.1d<-rma.mv(
  yi=effect,
  V=var, 
  # mods=~GCDType, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

summary(poll.mod.1d)


saveRDS(poll.mod.1d, file = "Models/pollutionMod.rds")


## ### THINKING ABOUT OTHER COVARIATES



mod.gsba <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + GSBA, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

summary(mod.gsba)
anova(mod.gsba, btt = "GSBA") # NOT SIGNIFICANT



# just large groups

taxa_dat <- hedges[hedges$GSBA %in% c("Acari", "Collembola",  "Earthworms", "Nematodes"),]
table(taxa_dat$driver, taxa_dat$GSBA)



mod.gsba.taxa <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * GSBA + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=taxa_dat)
anova(mod.gsba.taxa, btt = ":") #  SIGNIFICANT
anova(mod.gsba.taxa, btt = "Measurement") #  SIGNIFICANT

beep()
summary(mod.gsba.taxa)


saveRDS(mod.gsba.taxa, file = "Models/GSBAMod.rds")


t_dat <-predict(mod.gsba.taxa, newmods=rbind(c(0,0,0,0,0, 0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,0, 1,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,0, 0,1,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,0, 0,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                              # intercept (climate change)
                                              
                                              
                                                c(1,0,0,0,0, 0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(1,0,0,0,0, 1,0,0, 0,0,0, 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(1,0,0,0,0, 0,1,0, 0,0,0, 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
                                                c(1,0,0,0,0, 0,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
                                              #  HabitatLoss                                        
                                              
                                                c(0,1,0,0,0, 0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,1,0,0,0, 1,0,0, 0,0,0, 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,1,0,0,0, 0,1,0, 0,0,0, 0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                                                c(0,1,0,0,0, 0,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
                                              # Invasives 
                                              
                                                c(0,0,1,0,0, 0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,1,0,0, 1,0,0, 0,0,0, 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,1,0,0, 0,1,0, 0,0,0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                                                c(0,0,1,0,0, 0,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
                                              # LUI
                                              
                                                c(0,0,0,1,0, 0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,1,0, 1,0,0, 0,0,0, 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,1,0, 0,1,0, 0,0,0, 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
                                                c(0,0,0,1,0, 0,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
                                              # NutrientEnrichment
                                              
                                                c(0,0,0,0,1, 0,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,1, 1,0,0, 0,0,0, 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,1, 0,1,0, 0,0,0, 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                                                c(0,0,0,0,1, 0,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)),
                                              # Pollution 
                   addx=TRUE, digits=2) #

gcds <- rep(c("Climate", "habitat", "invasives", "lui", "nutrient", "pollution"), each = 4)
taxa <- rep(c("acari","collembola", "earthworms", "nematodes"), times = 6)

slabs <- paste(gcds, taxa)

par(mar=c(3, 8, 1, 1))
forest(t_dat$pred, sei=t_dat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))


## ### THINKING ABOUT OTHER COVARIATES - habitat

system_dat <- hedges[which(hedges$System != ""),]

mod.system <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * System + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=system_dat)
anova(mod.system, btt = ":") #  not significant
anova(mod.system, btt = "Measurement") #   significant


## Script to calculate effect sizes


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
  method="ML",
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

# Remove indices with little data
hedges <- hedges[!(hedges$Body.Size == ""),]
hedges <- hedges[!(hedges$Measurement == "FunctionalRichness"),]
hedges <- hedges[!(hedges$Measurement == "GroupAbundance"),]
hedges <- hedges[!(hedges$Measurement == "Simpson"),]
hedges <- hedges[!(hedges$Measurement == "Evenness"),]



hedges$Body.Size[which(hedges$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
hedges$Body.Size[which(hedges$Body.Size %in% c("Macro-arthropods", "Macro-invertebrates"))] <- "Macro-fauna"


write.csv(hedges, "Data/03_Data/HedgesData_cleaned.csv")



mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver:Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)

beep()


mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)
beep()
# mod.1 and mod.1b are significantly different. Indicating we need the interaction


qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
# not bad


mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)


anova(mod.1, mod.2) # not significantly different
beep()



mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)
anova(mod.2, mod.3) # that is significantly different


mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)
anova(mod.2, mod.4) # that is significantly different


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
                             c(0,0,0,1,0,1,0,0),
                             c(0,0,0,1,0,0,0,0),
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
climate$Body.Size[which(climate$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"


climate.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)

summary(climate.mod.1)

qqnorm(residuals(climate.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(climate.mod.1,type="pearson"),col="red")


climate.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)


anova(climate.mod.1, climate.mod.2) # not diff


climate.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)

anova(climate.mod.1, climate.mod.3) # not diff, but lower p-val

climate.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)

anova(climate.mod.1, climate.mod.4) # that's significant



climate.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)
anova(climate.mod.2, climate.mod.5) # not significant


climate.mod.6<-rma.mv(
  yi=effect,
  V=var, 
#  mods=~GCDType  , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)
anova(climate.mod.1, climate.mod.6) #  significant at 0.05 level



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



table(lui$GCDType, lui$Measurement)



## Checking body size
table(lui$Body.Size)
lui$Body.Size[which(lui$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"



lui.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)

summary(lui.mod.1)

qqnorm(residuals(lui.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(lui.mod.1,type="pearson"),col="red")


lui.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)


anova(lui.mod.1, lui.mod.2) #  signi. different (at 0.05 level) need body size?


lui.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)


anova(lui.mod.1, lui.mod.3) #  definitely significantly different. Need measurement


lui.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)


anova(lui.mod.1, lui.mod.4) #  definitely significantly different. Need driver type


## final model has all three
saveRDS(lui.mod.1, file = "Models/LUIMod.rds")



## NO IDEA HOW TO PRESENT THIS YET

table(lui$Body.Size) # macrofauna has the greatest numbers
# maybe show only macrofauna for the different measurement types


luidat <-predict(lui.mod.1, newmods=rbind(c(0,0,0,0,0,0,0,0,0,0,0,1,0,0),
                                              c(0,0,0,0,0,0,0,0,1,0,0,1,0,0),
                                              c(0,0,0,0,0,0,0,0,0,1,0,1,0,0),
                                              c(0,0,0,0,0,0,0,0,0,0,1,1,0,0),  # intercept
                                              
                                              c(1,0,0,0,0,0,0,0,0,0,0,1,0,0),
                                              c(1,0,0,0,0,0,0,0,1,0,0,1,0,0),
                                              c(1,0,0,0,0,0,0,0,0,1,0,1,0,0),
                                              c(1,0,0,0,0,0,0,0,0,0,1,1,0,0),# fire
                                              
                                              c(0,1,0,0,0,0,0,0,0,0,0,1,0,0),
                                              c(0,1,0,0,0,0,0,0,1,0,0,1,0,0),
                                              c(0,1,0,0,0,0,0,0,0,1,0,1,0,0),
                                              c(0,1,0,0,0,0,0,0,0,0,1,1,0,0),  # grazing   
                                              
                                              c(0,0,1,0,0,0,0,0,0,0,0,1,0,0),
                                              c(0,0,1,0,0,0,0,0,1,0,0,1,0,0),
                                              c(0,0,1,0,0,0,0,0,0,1,0,1,0,0),
                                              c(0,0,1,0,0,0,0,0,0,0,1,1,0,0),  # Harvesting
                                              
                                              c(0,0,0,1,0,0,0,0,0,0,0,1,0,0),
                                              c(0,0,0,1,0,0,0,0,1,0,0,1,0,0),
                                              c(0,0,0,1,0,0,0,0,0,1,0,1,0,0),
                                              c(0,0,0,1,0,0,0,0,0,0,1,1,0,0), # Irrigation 
                                              
                                              c(0,0,0,0,1,0,0,0,0,0,0,1,0,0),
                                              c(0,0,0,0,1,0,0,0,1,0,0,1,0,0),
                                              c(0,0,0,0,1,0,0,0,0,1,0,1,0,0),
                                              c(0,0,0,0,1,0,0,0,0,0,1,1,0,0),# Management
                                              
                                              c(0,0,0,0,0,1,0,0,0,0,0,1,0,0),
                                              c(0,0,0,0,0,1,0,0,1,0,0,1,0,0),
                                              c(0,0,0,0,0,1,0,0,0,1,0,1,0,0),
                                              c(0,0,0,0,0,1,0,0,0,0,1,1,0,0),  # Mono
                                              
                                              c(0,0,0,0,0,0,1,0,0,0,0,1,0,0),
                                              c(0,0,0,0,0,0,1,0,1,0,0,1,0,0),
                                              c(0,0,0,0,0,0,1,0,0,1,0,1,0,0),
                                              c(0,0,0,0,0,0,1,0,0,0,1,1,0,0), # organic
                                              
                                              c(0,0,0,0,0,0,0,1,0,0,0,1,0,0),
                                              c(0,0,0,0,0,0,0,1,1,0,0,1,0,0),
                                              c(0,0,0,0,0,0,0,1,0,1,0,1,0,0),
                                              c(0,0,0,0,0,0,0,1,0,0,1,1,0,0)# tillage
), addx=TRUE, digits=2) #




slabs <- rep(c("abundance", "biomass", "richness", "shannon"), times = 9)
par(mar=c(3, 8, 1, 1))
forest(luidat$pred, sei=luidat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))
abline(h=4.5, b=0, lty= 2)
abline(h=8.5, b=0, lty= 2)
abline(h=12.5, b=0, lty= 2)
abline(h=16.5, b=0, lty= 2)
abline(h=20.5, b=0, lty= 2)
abline(h=24.5, b=0, lty= 2)
abline(h=28.5, b=0, lty= 2)
abline(h=32.5, b=0, lty= 2)


mtext("Degradation", side = 2, line = 3, at = 34.5, cex = 0.5)
mtext("Fire", side = 2, line = 3, at = 30.5, cex = 0.7)
mtext("Grazing", side = 2, line = 3, at = 26.4, cex = 0.7)
mtext("Harvesting", side = 2, line = 3, at = 22.4, cex = 0.5)
mtext("Irrigation", side = 2, line = 3, at = 18.4, cex = 0.7)
mtext("Management", side = 2, line = 3,at = 14.4, cex = 0.5)
mtext("Mono-\nculture", side = 2, line = 3, at = 10.4, cex = 0.7)
mtext("Inorganic", side = 2, line = 3, at = 6.4, cex = 0.6)
mtext("Tillage", side = 2, line = 3, at = 2.4, cex = 0.7)






## NUTRIENT ENRICHMENT -------



nutri <- hedges[which(hedges$driver == "NutrientEnrichment"),] # 820

table(nutri$GCDType)


table(nutri$GCDType, nutri$Measurement)


## Check body size
table(nutri$Body.Size)
nutri$Body.Size[which(nutri$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"



nutri.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

summary(nutri.mod.1)


nutri.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

anova(nutri.mod.1, nutri.mod.2) #  not different (at 0.05 level)


nutri.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)
anova(nutri.mod.1, nutri.mod.3) #  significantly different (at 0.05 level)


nutri.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)
anova(nutri.mod.1, nutri.mod.4) #  significantly different (at 0.05 level)
# move to mod2


nutri.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)
anova(nutri.mod.2, nutri.mod.5) #  not significantly different (at 0.05 level)



nutri.mod.6<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

anova(nutri.mod.2, nutri.mod.6) #  significantly different



anova(nutri.mod.1, nutri.mod.5) # but that is sign. different. So that makes no sense

## Need mod 5 in theory. But it doesn't quite make sense
# Clearly as both reductions were on the cusp of being significant
# the full reduciton is too much
## we will go with mod2


saveRDS(nutri.mod.2, file = "Models/nutriMod.rds")



nutridat <-predict(nutri.mod.2, newmods=rbind(c(0,0,0,0,0,0,0,0,0,0,0),
                                    c(0,0,0,0,0,0,0,0,1,0,0),
                                    c(0,0,0,0,0,0,0,0,0,1,0),
                                    c(0,0,0,0,0,0,0,0,0,0,1),  # intercept
                                    
                                    c(1,0,0,0,0,0,0,0,0,0,0),
                                    c(1,0,0,0,0,0,0,0,1,0,0),
                                    c(1,0,0,0,0,0,0,0,0,1,0),
                                    c(1,0,0,0,0,0,0,0,0,0,1),# Ca-liming + Wood ash 
                                    
                                    c(0,1,0,0,0,0,0,0,0,0,0),
                                    c(0,1,0,0,0,0,0,0,1,0,0),
                                    c(0,1,0,0,0,0,0,0,0,1,0),
                                    c(0,1,0,0,0,0,0,0,0,0,1),  # Compost   
                                    
                                    c(0,0,1,0,0,0,0,0,0,0,0),
                                    c(0,0,1,0,0,0,0,0,1,0,0),
                                    c(0,0,1,0,0,0,0,0,0,1,0),
                                    c(0,0,1,0,0,0,0,0,0,0,1),  # Manure + Slurry
                                    
                                    c(0,0,0,1,0,0,0,0,0,0,0),
                                    c(0,0,0,1,0,0,0,0,1,0,0),
                                    c(0,0,0,1,0,0,0,0,0,1,0),
                                    c(0,0,0,1,0,0,0,0,0,0,1), # Mixture 
                                    
                                    c(0,0,0,0,1,0,0,0,0,0,0),
                                    c(0,0,0,0,1,0,0,0,1,0,0),
                                    c(0,0,0,0,1,0,0,0,0,1,0),
                                    c(0,0,0,0,1,0,0,0,0,0,1),# Other Organic fertilisers
                                    
                                    c(0,0,0,0,0,1,0,0,0,0,0),
                                    c(0,0,0,0,0,1,0,0,1,0,0),
                                    c(0,0,0,0,0,1,0,0,0,1,0),
                                    c(0,0,0,0,0,1,0,0,0,0,1),  # residue + Mulch
                                    
                                    c(0,0,0,0,0,0,1,0,0,0,0),
                                    c(0,0,0,0,0,0,1,0,1,0,0),
                                    c(0,0,0,0,0,0,1,0,0,1,0),
                                    c(0,0,0,0,0,0,1,0,0,0,1), # Sludge
                                    
                                    c(0,0,0,0,0,0,0,1,0,0,0),
                                    c(0,0,0,0,0,0,0,1,1,0,0),
                                    c(0,0,0,0,0,0,0,1,0,1,0),
                                    c(0,0,0,0,0,0,0,1,0,0,1)# Synthetic Fertilizers
), addx=TRUE, digits=2) #


slabs <- rep(c("abundance", "biomass", "richness", "shannon"), times = 9)
par(mar=c(3, 8, 1, 1))
forest(nutridat$pred, sei=nutridat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))
abline(h=4.5, b=0, lty= 2)
abline(h=8.5, b=0, lty= 2)
abline(h=12.5, b=0, lty= 2)
abline(h=16.5, b=0, lty= 2)
abline(h=20.5, b=0, lty= 2)
abline(h=24.5, b=0, lty= 2)
abline(h=28.5, b=0, lty= 2)
abline(h=32.5, b=0, lty= 2)


mtext("Biochar", side = 2, line = 3, at = 34.5, cex = 0.7)
mtext("Liming", side = 2, line = 3, at = 30.5, cex = 0.7)
mtext("Compost", side = 2, line = 3, at = 26.4, cex = 0.7)
mtext("Manure \n+ \nSlurry", side = 2, line = 3, at = 22.4, cex = 0.7)
mtext("Mixture", side = 2, line = 3, at = 18.4, cex = 0.7)
mtext("Organic \nfertilisers", side = 2, line = 3,at = 14.4, cex = 0.7)


mtext("Residue \n+ Mulch", side = 2, line = 3, at = 10.4, cex = 0.7)
mtext("Sludge", side = 2, line = 3, at = 6.4, cex = 0.7)
mtext("Synthetic \nFertilizers", side = 2, line = 3, at = 2.4, cex = 0.7)








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
invas$Body.Size[which(invas$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"



invas.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)


invas.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)

anova(invas.mod.1, invas.mod.2) ## That is not quite significant


invas.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)

anova(invas.mod.1, invas.mod.3) ## That is also not significant



invas.mod.4<-rma.mv(
  yi=effect,
  V=var, 
 #  mods=~Body.Size , ## intercept only
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)

anova(invas.mod.1, invas.mod.4) ## That is  not significant


summary(invas.mod.4)
saveRDS(invas.mod.4, file = "Models/invasiveMod.rds")




### POLLUTION
poll <- hedges[which(hedges$driver == "Pollution"),] # 850

# pollution type
table(poll$GCDType)

poll$GCDType[which(poll$GCDType %in% c("Metals, Radionuclides", "Metals; PAH", "Metals; PCBs; PAHs", "Pesticides,Metals"))] <- "Mixture"
poll <- poll[which(poll$GCDType != ""),] 
poll <- poll[which(poll$GCDType != "Antibiotics"),] 
poll <- poll[which(poll$GCDType != "Chloride, Sodium"),] 
poll <- poll[which(poll$GCDType != "endocrine disruptors"),] 
poll <- poll[which(poll$GCDType != "Salinization"),] 
poll <- poll[which(poll$GCDType != "Sulphate"),] 


# measurement
table(poll$GCDType, poll$Measurement)

# body size
table(poll$Body.Size)


poll$Body.Size[which(poll$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"



poll.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)

summary(poll.mod.1)


poll.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)

anova(poll.mod.1, poll.mod.2) # not Significant at 0.05 level

poll.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.1, poll.mod.3) # not Significantly different


poll.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Body.Size + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)

anova(poll.mod.1, poll.mod.4) # not Significantly different


# carry on with mod3



poll.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.3, poll.mod.5) # that is not Significantly different



poll.mod.6<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.3, poll.mod.6) # that is  Significantly different

## carry on with  poll.mod.6


poll.mod.7<-rma.mv(
  yi=effect,
  V=var, 
  # mods=~Measurement, ## intercept only
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.6, poll.mod.7) # not significantly different


anova(poll.mod.1, poll.mod.7) # not different at all



## The type of GCD doesn't matter, just the measurement

summary(poll.mod.7)
# the itnercept is significant. 
# So this is a different result than compared to the invasive species


saveRDS(poll.mod.7, file = "Models/pollutionMod.rds")








## The code created a colours qqplot
## Coloured qq plot https://stackoverflow.com/questions/42678858/q-q-plot-with-ggplot2stat-qq-colours-single-group
library(broom) ## for augment()
dda <- cbind(augment(mod.1),f=EffectSizes$driver)
dda = cbind(dda, setNames(qqnorm(dda$.resid, plot.it=FALSE), c("Theoretical", "Sample")))

ggplot(dda) + 
  geom_point(aes(x=Theoretical, y=Sample, colour=f))



## This code for creating a basic plot

y<-summary(mod.1)$b
ci_l<-summary(mod.1)$ci.lb
ci_h<-summary(mod.1)$ci.ub

fg1<-data.frame(cbind(y,ci_l,ci_h))
colnames(fg1)[1]<-"y"
colnames(fg1)[2]<-"ci_l"
colnames(fg1)[3]<-"ci_h"
fg1$GCD<-c("Climate change","Fragmentation", "Invasive species",
           "LUI", "Nutrient enrichment", "Pollution", rep(NA, 7))
fg1$GCD<-as.factor(fg1$GCD)


p <- ggplot(fg1, aes(x=GCD, y=y, ymin=ci_l, ymax=ci_h))+
  geom_point(aes(size = 1)) +
  geom_pointrange()+
  geom_hline(yintercept = 0, linetype=2)+
  coord_flip()+
  theme_bw() +
  theme(axis.text=element_text(size=16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab('') +
  ylab ('Effect Size')
p


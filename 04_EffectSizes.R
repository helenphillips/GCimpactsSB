

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




hedges <- escalc(measure = "SMD", # 
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
hedges[which(hedges$effect < -50),]
hedges[which(hedges$effect > 50),]

## remove outliers as they do cause issues, especially in pub bias
hedges <- hedges[which(hedges$effect > -50),]
hedges <- hedges[which(hedges$effect < 50),]
## removes only four cases



head(hedges[which(hedges$var > 2),])
table(hedges$Measurement[which(hedges$var > 10)])
head(hedges[which(hedges$var > 10),])





# Save data
write.csv(hedges, "Data/03_Data/HedgesData.csv")
# 3341


hedges$UniqueID <- paste(hedges$ID, hedges$Case_ID, hedges$driver)
hedges$ID <- as.factor(hedges$ID) # to be on the safe side






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


## TRYING TO ESTABLISH WHETHER ALL EFFECT SIZES CAN BE USED ----------


measurement.mod.1 <-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

anova(measurement.mod.1, btt = "Measurement") # is significant

saveRDS(measurement.mod.1, file = "Models/MeasurementMod.rds")


## But all responses are negative


table(hedges$Measurement, hedges$driver)
table(hedges$GCDType, hedges$Measurement)


## the fact that the data points are spare across other metrics (compared to abundance)
## means I am hesitant to interpret the results if we have measurement as a fixed
## effect.
## The additive effect doesn't really show what may be happening
## But as there is an impact of the metric type, having it as a random effect,,
## thus averaging over the impact (keeping in mind, that all are negative anyway), 
## may be the best compromise.

qqnorm(residuals(measurement.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(measurement.mod.1,type="pearson"),col="red")
# not bad




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














## PUBLICATION BIAS --------


test.egger <- rma.mv(yi,vi, mod = vi, random = list(~1|ID/UniqueID, ~ 1 | Measurement), data = hedges)
saveRDS(test.egger, file = "Models/Eggersmod.rds")


## from lea's code
# standardized effect sizes (E*) 
std.es.B <- hedges$effect/sqrt(hedges$var)

# standard errors (SE)
se.es.B <-  sqrt(hedges$var)

# precision (1/SE)
precision.es.B <- 1/se.es.B

# Egger's regressions : E* vs. 1/SE and incoporating the main covariate
eggerreg <- rma.mv(std.es.B, var, 
                   random =  list(~1|ID/UniqueID, ~ 1 | Measurement),
                   data = hedges,
                   mods =~ driver * precision.es.B)
saveRDS(eggerreg, file = "Models/Eggersmod_lea.rds")

plot(std.es.B ~ precision.es.B)






plot(precision.es.B ~ hedges$effect)





## Using Lea's code for each of the GCD individual models.

eggersdat <- hedges
eggersdat$std.es.B <- eggersdat$effect/sqrt(eggersdat$var)
eggersdat$se.es.B <-  sqrt(eggersdat$var)
eggersdat$precision.es.B <- 1/eggersdat$se.es.B


# climate
climate_eggersdat <- eggersdat[which(eggersdat$driver == "Climate"),] # 462
climate_eggersdat <- climate_eggersdat[which(climate_eggersdat$GCDType != "Vegetation"),]
climate_eggersdat <- climate_eggersdat[which(climate_eggersdat$GCDType != "Precipitation+Temp"),]
climate_eggersdat <- climate_eggersdat[which(climate_eggersdat$GCDType != "UVB Radiation"),]
climate_eggersdat <- droplevels(climate_eggersdat[which(climate_eggersdat$Body.Size != "All sizes"),])


climate_eggerreg <- rma.mv(std.es.B, var, 
                   random =  list(~1|ID/UniqueID, ~ 1 | Measurement),
                   data = climate_eggersdat,
                   mods =~ GCDType * precision.es.B)
saveRDS(climate_eggerreg, file = "Models/ClimateEggersmod_lea.rds")
# climate_eggerreg <- readRDS("Models/ClimateEggersmod_lea.rds")

# lui
lui_eggersdat <- eggersdat[which(eggersdat$driver == "LUI"),] # 911


lui_eggersdat <- droplevels(lui_eggersdat[which(lui_eggersdat$Body.Size != "All sizes"),])
lui_eggersdat$GCDType[which(lui_eggersdat$GCDType %in% c("Defoliation"))] <- "Grazing"
lui_eggersdat$GCDType[which(lui_eggersdat$GCDType %in% c("Intensity", "Landscape", "Management", "Weeding", "Planting", "Intensification"))] <- "Management"
lui_eggersdat$GCDType[which(lui_eggersdat$GCDType %in% c("Logging"))] <- "Harvesting"
lui_eggersdat$GCDType[which(lui_eggersdat$GCDType %in% c("Water"))] <- "Irrigation"
lui_eggersdat$GCDType[which(lui_eggersdat$GCDType %in% c("degradation", "Disturbance"))] <- "Degradation"
lui_eggersdat <- lui_eggersdat[which(lui_eggersdat$GCDType != "Human population"),]
lui_eggersdat <- lui_eggersdat[which(lui_eggersdat$GCDType != "Mono- versus poly-culture"),] # too little
lui_eggersdat <- lui_eggersdat[which(lui_eggersdat$GCDType != "Irrigation"),] # too little
lui_eggersdat <- lui_eggersdat[which(lui_eggersdat$GCDType != "Management"),] # too little
lui_eggersdat <- lui_eggersdat[which(lui_eggersdat$GCDType != "Degradation"),] # too little
lui_eggersdat$GCDType <- relevel(as.factor(lui_eggersdat$GCDType), ref = "Grazing")



lui_eggerreg <- rma.mv(std.es.B, var, 
                           random =  list(~1|ID/UniqueID, ~ 1 | Measurement),
                           data = lui_eggersdat,
                           mods =~ GCDType * precision.es.B)
saveRDS(lui_eggerreg, file = "Models/LUIEggersmod_lea.rds")


# nutrient
nutri_eggersdat <- eggersdat[which(eggersdat$driver == "NutrientEnrichment"),] # 820

nutri_eggersdat <- nutri_eggersdat[which(nutri_eggersdat$GCDType != "Biochar"),]
nutri_eggersdat$GCDType <- as.factor(nutri_eggersdat$GCDType)
nutri_eggersdat$GCDType <- relevel(nutri_eggersdat$GCDType, ref = "Synthetic Fertilizers")



nutri_eggerreg <- rma.mv(std.es.B, var, 
                       random =  list(~1|ID/UniqueID, ~ 1 | Measurement),
                       data = nutri_eggersdat,
                       mods =~ GCDType * precision.es.B)
saveRDS(nutri_eggerreg, file = "Models/NutriEggersmod_lea.rds")


# invasives
invas_eggersdat <- eggersdat[which(eggersdat$driver == "Invasives"),] # 188
invas_eggersdat <- invas_eggersdat[which(invas_eggersdat$GCDType != "Plants-mixture"),]

invas_eggerreg <- rma.mv(std.es.B, var, 
                         random =  list(~1|ID/UniqueID, ~ 1 | Measurement),
                         data = invas_eggersdat,
                         mods =~ GCDType * precision.es.B)
saveRDS(invas_eggerreg, file = "Models/InvasEggersmod_lea.rds")


# pollution
poll_eggersdat <- eggersdat[which(eggersdat$driver == "Pollution"),] # 850


fulltext <- read.csv("Data/April2022/Full Text Screening - Sheet1.csv")
fulltext <- fulltext[,c('PaperID', 'PollutionSource')]

poll_eggersdat <- merge(poll_eggersdat, fulltext, by.x = "ID", by.y = "PaperID", all.x = TRUE)

poll_eggersdat <- poll_eggersdat[which(poll_eggersdat$GCDType %in% c("Metals", "Pesticides")),] 


## previously missing
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 105)] <- "Agricultural"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 2011)] <- "Agricultural"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 3372)] <- "Others"

## changing to remove multiple sources
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 395)] <- "Mining/Smelting"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 2433 )] <- "Waste/sewage"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 1648)] <- "Mining/Smelting"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 1850)] <- "Mining/Smelting"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 564)] <- "Industrial"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 2509)] <- "Industrial"
poll_eggersdat$PollutionSource[which(poll_eggersdat$ID == 2715)] <- "Urban/transport"


poll_eggersdat$PollutionSource[which(poll_eggersdat$PollutionSource == "Agricultural")] <- "Farming"
poll_eggersdat$PollutionSource[which(poll_eggersdat$PollutionSource == "Agricultural/livestock")] <- "Farming"
poll_eggersdat$PollutionSource[which(poll_eggersdat$PollutionSource == "Military/wars")] <- "Others"
poll_eggersdat$PollutionSource[which(poll_eggersdat$PollutionSource == "Natural/geogenic")] <- "Others"


poll_eggersdat$GCDTypeSource <- paste(poll_eggersdat$GCDType, "-", poll_eggersdat$PollutionSource)
poll_eggersdat$GCDTypeSource[which(poll_eggersdat$GCDTypeSource == "Pesticides - Industrial")] <- "Pesticides - Others"

poll_eggerreg <- rma.mv(std.es.B, var, 
                         random =  list(~1|ID/UniqueID, ~ 1 | Measurement),
                         data = poll_eggersdat,
                         mods =~ GCDTypeSource * precision.es.B)
saveRDS(poll_eggerreg, file = "Models/PollEggersmod_lea.rds")



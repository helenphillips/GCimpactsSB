
library(Hmisc)
library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished
library(emmeans)
library(dplyr)

setwd("C:/Users/helenp/WORK/GCimpactsSB")

hedges <- read.csv("Data/03_Data/HedgesData_cleaned_June2023.csv")




## LUI -------

lui <- hedges[which(hedges$driver == "LUI"),] 
lui <- droplevels(lui[which(lui$Body.Size != "All sizes"),])
lui$GCDType[which(lui$GCDType %in% c("Defoliation"))] <- "Grazing"
lui$GCDType[which(lui$GCDType %in% c("Intensity", "Landscape", "Management", "Weeding", "Planting", "Intensification"))] <- "Management"
lui$GCDType[which(lui$GCDType %in% c("Logging"))] <- "Harvesting"
lui$GCDType[which(lui$GCDType %in% c("Water"))] <- "Irrigation"
lui$GCDType[which(lui$GCDType %in% c("degradation", "Disturbance"))] <- "Degradation"
lui <- lui[which(lui$GCDType != "Human population"),]
lui <- lui[which(lui$GCDType != "Mono- versus poly-culture"),] # too little
lui <- lui[which(lui$GCDType != "Irrigation"),] # too little
lui <- lui[which(lui$GCDType != "Management"),] # too little
lui <- lui[which(lui$GCDType != "Degradation"),] # too little

lui$GCDType <- relevel(as.factor(lui$GCDType), ref = "Grazing")



lui.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=lui)


anova(lui.mod.1, btt = ":") # not significant





lui.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=lui)



anova(lui.mod.2, btt = "GCD") # significant
anova(lui.mod.2, btt = "Size") #  significant

summary(lui.mod.2)


## All the results are consistent with the original model / method
# microfauna still diff from the macro (meso not)
# everything is causing a decline
# grazing, inorganic, and tillage significant
# organic the most
# essentially the same

## final model
saveRDS(lui.mod.2, file = "Models/LUIMod_redo_june2023_KH.rds")



## Climate

climate <- hedges[which(hedges$driver == "Climate"),] # 

table(climate$GCDType)
table(climate$GCDType, climate$Measurement)


## going to remove vegetation as a gcd type. 2 lines are thickness of moss, and the 10 other lines are canopy gaps due to ice storm
climate <- climate[which(climate$GCDType != "Vegetation"),]

# For now will also remove the precip+temp
climate <- climate[which(climate$GCDType != "Precipitation+Temp"),]
climate <- climate[which(climate$GCDType != "UVB Radiation"),]
climate <- climate[which(climate$GCDType != "Gas - O3"),]

# Checking body size
table(climate$Body.Size)
# removing the all sizes
climate <- droplevels(climate[which(climate$Body.Size != "All sizes"),])


table(climate$GCDType, climate$GSBA) # to check for taxonomic model
# Acari, collembola and nematodes are definitrely fine. Earthworms less so



climate.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=climate)

summary(climate.mod.1)

anova(climate.mod.1, btt = ":") # should be testing the levels that are the interactions
# Not significant


climate.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=climate)

anova(climate.mod.1b, btt = "GCD") # significant
anova(climate.mod.1b, btt = "Size") # not significant


climate.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="ML",
  test = "t",
  digits=4,
  data=climate)
anova(climate.mod.5, btt = "GCD") # significant


summary(climate.mod.5) #No difference

## Use climate.mod.5
saveRDS(climate.mod.5, file = "Models/ClimateMod_june2023_KH.rds")




## Nutrient Enrichment

nutri <- hedges[which(hedges$driver == "NutrientEnrichment"),] # 820

table(nutri$GCDType)


table(nutri$GCDType, nutri$Measurement)


# Removing some GCDs and relevelling
nutri <- nutri[which(nutri$GCDType != "Biochar"),]
nutri$GCDType <- as.factor(nutri$GCDType)
nutri$GCDType <- relevel(nutri$GCDType, ref = "Synthetic Fertilizers")

table(nutri$GCDType, nutri$GSBA)

## Check body size
table(nutri$Body.Size)

table(nutri$GCDType, nutri$GSBA) # to check for taxonomic model
# acari, earthworms and nematodes are definitely fine. collembola could be better, but not bad

nutri.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=nutri)


anova(nutri.mod.1, btt = ":") # not significant


nutri.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=nutri)

anova(nutri.mod.1b, btt = "GCD") #  significant
anova(nutri.mod.1b, btt = "Size") # not significant


nutri.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=nutri)

summary(nutri.mod.2)
anova(nutri.mod.2, btt = "GCD") #  significant


saveRDS(nutri.mod.2, file = "Models/nutriMod_june2023_KH.rds")

## no difference in the two methods



## POLLUTION



poll <- readRDS(file = "Models/pollutionDataFrame_june2023.rds")


poll.mod.20<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=poll)
anova(poll.mod.20, btt = ":") #   not significant


poll.mod.21<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=poll)

anova(poll.mod.21, btt = "Size") # not  significant
anova(poll.mod.21, btt = "GCD") #   significant


poll.mod.22<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=poll)
anova(poll.mod.22, btt = "GCD") #   significant

summary(poll.mod.22)


saveRDS(poll.mod.22, file = "Models/pollutionMod_june2023_KH.rds")



## Invasives

invas <- hedges[which(hedges$driver == "Invasives"),] # 188

invas <- invas[which(invas$GCDType != "Plants-mixture"),]



invas.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=invas)
anova(invas.mod.1, btt = ":") # not  significant


invas.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=invas)



anova(invas.mod.2, btt = "GCD") # not  significant
anova(invas.mod.2, btt = "Size") # not  significant




invas.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  #  mods=~GCDType + Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=invas)

summary(invas.mod.3) # not significant intercapt
saveRDS(invas.mod.3, file = "Models/invasiveMod_june2023_KH.rds")



## Taxonomic model

taxa_dat <- hedges[hedges$GSBA %in% c("Acari", "Collembola",  "Earthworms", "Nematodes"),]



mod.gsba.taxa <-rma.mv(
  yi=effect,
  V=var, 
  mods=~ driver * GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=taxa_dat)
anova(mod.gsba.taxa, btt = ":") #  SIGNIFICANT
beep()

summary(mod.gsba.taxa)


saveRDS(mod.gsba.taxa, file = "Models/GSBAMod_june2023_KH.rds")


## ph mod
phdat <- hedges[which(!(is.na(hedges$ph_fixed))),] # 1715
phdat <- phdat[which(phdat$driver != "HabitatLoss"),] # 1712



mod.ph <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * ph_fixed, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=phdat)

mod.ph2 <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + ph_fixed, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=phdat)
saveRDS(mod.ph2, file = "Models/pHMod_june2023_KH.rds")


## habitat type


system_dat <- hedges[which(hedges$System != ""),]

mod.system <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * System, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=system_dat)
anova(mod.system, btt = ":") #  not significant



mod.system2 <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + System, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test = "t",
  digits=4,
  data=system_dat)
anova(mod.system2, btt = "System") #  not significant
anova(mod.system2, btt = "driver") #   significant



saveRDS(mod.system2, file = "Models/SystemMod_june2023_KH.rds")


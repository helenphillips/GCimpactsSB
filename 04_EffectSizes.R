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



head(hedges[which(hedges$var > 2),])
table(hedges$Measurement[which(hedges$var > 10)])
head(hedges[which(hedges$var > 10),])


# TODO: check the high vars with Lea



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


measurement.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

anova(measurement.mod.1, btt = "Measurement") # is significant



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


## FULL MODEL ---------




# If fittiung with REML (which is recommended) then can't compare models 
# (https://stats.stackexchange.com/questions/48671/what-is-restricted-maximum-likelihood-and-when-should-it-be-used)
# https://stats.stackexchange.com/questions/517857/testing-the-effect-of-moderators-in-metafor-package
# so test moderators using the anova(singleMod) way




write.csv(hedges, "Data/03_Data/HedgesData_cleaned.csv")



## Publication bias




############################
## FULL ANALYSIS WITH CROSSED RANDOM EFFECTS
########################





mod.1<-rma.mv(
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

anova(mod.1, btt = ":") # should be testing the levels that are the interactions
## not significant

# Ultimately, this interaction is not needed


mod.1b<-rma.mv(
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


anova(mod.1b, btt = "driver") # significant
anova(mod.1b, btt = "Body") # not significant



qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
# not bad


mod.2<-rma.mv(
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

anova(mod.2, btt = "driver") # significant




## we stick with mod.2
saveRDS(mod.2, file = "Models/MainMod.rds")



## Can use this to look at publication bias

test <- hedges[which(hedges$vi < 5),]
funnel(x = test$effect, vi = test$var, ni = test$Control_N,
       yaxis="sei")

funnel(x = test$effect, vi = test$var, ni = test$Control_N,
       yaxis="vi")




funnel(mod.2, main="Standard Error")
funnel(mod.2, yaxis="vi", main="Sampling Variance")
funnel(mod.2, yaxis="seinv", main="Inverse Standard Error")
funnel(mod.2, yaxis="vinv", main="Inverse Sampling Variance")



# https://stats.stackexchange.com/questions/109841/presenting-results-of-a-meta-analysis-with-multiple-moderators
### compute the estimated/predicted correlation for each combination
test <-predict(mod.2, newmods=rbind(c(0,0,0,0,0),
  # climate with all different measurement types
                             c(1,0,0,0,0),# habitat 
                             c(0,1,0,0,0), # invasive
                             c(0,0,1,0,0),# lui
                             c(0,0,0,1,0),# nutrient
                             c(0,0,0,0,1)# pollution
                             ), addx=TRUE, digits=2) #

slabs <- c("Climate", "Fragmentation", "Invasives", "LUI", "Nutrient", "Pollution")
par(mar=c(3, 8, 1, 1))
forest(test$pred, sei=test$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))





## SPLITTING DRIVERS --------

climate <- hedges[which(hedges$driver == "Climate"),] # 462

table(climate$GCDType)
table(climate$GCDType, climate$Measurement)


## going to remove vegetation as a gcd type. 2 lines are thickness of moss, and the 10 other lines are canopy gaps due to ice storm
climate <- climate[which(climate$GCDType != "Vegetation"),]

# For now will also remove the precip+temp
climate <- climate[which(climate$GCDType != "Precipitation+Temp"),]
climate <- climate[which(climate$GCDType != "UVB Radiation"),]

# Checking body size
table(climate$Body.Size)
# removing the all sizes
climate <- droplevels(climate[which(climate$Body.Size != "All sizes"),])



climate.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
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
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
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
  digits=4,
  data=climate)


## Use climate.mod.5
saveRDS(climate.mod.5, file = "Models/ClimateMod.rds")


qqnorm(residuals(climate.mod.5,type="pearson"),main="QQ plot: residuals")
qqline(residuals(climate.mod.5,type="pearson"),col="red") # fine




climatedat <-predict(climate.mod.5, newmods=rbind(c(0,0,0),
                                    c(1,0,0),
                                    c(0,1,0),
                                    c(0,0,1)
), addx=TRUE, digits=2) #

slabs <- c("Gas", "Temperature",
           "Water Availability-Drought", "Water Availability-Flood")
par(mar=c(3, 8, 1, 1))
forest(climatedat$pred, sei=climatedat$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))


## LUI -------

lui <- hedges[which(hedges$driver == "LUI"),] # 911


lui <- droplevels(lui[which(lui$Body.Size != "All sizes"),])



table(lui$GCDType)
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



table(lui$GCDType, lui$Measurement)


lui$GCDType <- relevel(as.factor(lui$GCDType), ref = "Grazing")



## Checking body size
table(lui$Body.Size)



lui.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
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
  digits=4,
  data=lui)



anova(lui.mod.2, btt = "GCD") # significant
anova(lui.mod.2, btt = "Size") #  significant

summary(lui.mod.1b)





## final model
saveRDS(lui.mod.2, file = "Models/LUIMod.rds")




luidat <-predict(lui.mod.2, newmods=rbind(c(0,0,0,0,0,0),  # intercept
                                          c(0,0,0,0,1,0),
                                          c(0,0,0,0,0,1),
                                          
                                              c(1,0,0,0,0,0),# fire
                                          c(1,0,0,0,1,0),
                                          c(1,0,0,0,0,1),
                                      
                                              c(0,1,0,0,0,0), 
                                          c(0,1,0,0,1,0),
                                          c(0,1,0,0,0,1),# Harvesting   
                                       
                                              c(0,0,1,0,0,0),
                                          c(0,0,1,0,1,0),
                                          c(0,0,1,0,0,1),# organic
                                              
                                              c(0,0,0,1,0,0),
                                          c(0,0,0,1,1,0),
                                          c(0,0,0,1,0,1)# tillage
), addx=TRUE, digits=2) #

macro <- c(1, 4, 7, 10, 13) #3 just the macrofauna


slabs <- c("Grazing", "Fire", "Harvesting", "Inorganic", "Tillage")
par(mar=c(3, 8, 1, 1))
forest(luidat$pred[macro], sei=luidat$se[macro], slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))



#3 Just showing the impact of body size

measurement_coefs <- data.frame(meas = c("Macro-fauna", "Meso-fauna", "Micro-fauna"),
                                coef = c(0, lui.mod.2$beta[6:7]), 
                                ses = c(0, lui.mod.2$se[6:7]))


measurement_coefs <- measurement_coefs[c(3,2,1),]

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
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(nutri.mod.1)


anova(nutri.mod.1, btt = ":") # not significant


nutri.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
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
  digits=4,
  data=nutri)

summary(nutri.mod.2)
anova(nutri.mod.2, btt = "GCD") #  significant


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



## Body size
table(invas$GCDType, invas$Body.Size)



invas.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
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
  digits=4,
  data=invas)

summary(invas.mod.3) # not significant intercapt
saveRDS(invas.mod.3, file = "Models/invasiveMod.rds")




### POLLUTION
poll <- hedges[which(hedges$driver == "Pollution"),] # 850


fulltext <- read.csv("Data/April2022/Full Text Screening - Sheet1.csv")
fulltext <- fulltext[,c('PaperID', 'PollutionSource')]

poll <- merge(poll, fulltext, by.x = "ID", by.y = "PaperID", all.x = TRUE)

# pollution type
table(poll$GCDType)

poll$GCDType[which(poll$GCDType %in% c("Chloride, Sodium", "Metals, Radionuclides", "Metals; PAH", 
                          "Metals; PCBs; PAHs", "Pesticides,Metals"))] <- "Mixture"
  

aggregate(poll$ID, by = list(GCD = poll$GCDType), function(x) length(unique(x)))


poll <- poll[which(poll$GCDType != ""),]
poll <- poll[which(poll$GCDType != "Antibiotics"),]
poll <- poll[which(poll$GCDType != "endocrine disruptors"),]
poll <- poll[which(poll$GCDType != "Salinization"),]
poll <- poll[which(poll$GCDType != "Sulphate"),]
poll <- poll[which(poll$GCDType != "Nanoparticles"),]


poll <- poll[which(poll$Body.Size != "All sizes"),]


# poll <- poll[which(poll$GCDType %in% c("Metals", "Pesticides")),] 

poll$Body.Size <- as.factor(poll$Body.Size)
poll$Body.Size <- relevel(poll$Body.Size, ref = "Meso-fauna")


aggregate(poll$ID, 
          by = list(body = poll$Body.Size, GCD = poll$GCDType), 
          function(x) length(unique(x)))



# pollution type
table(poll$GCDType)




# measurement
table(poll$GCDType, poll$Measurement)

# body size
table(poll$Body.Size)
table(poll$GCDType, poll$Body.Size)



poll.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

anova(poll.mod.1, btt = ":") #  significant
summary(poll.mod.1)

poll.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
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
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)
anova(poll.mod.1c, btt = "GCD") #  nearly significant






poll.mod.1d<-rma.mv(
  yi=effect,
  V=var, 
  # mods=~GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

summary(poll.mod.1d)


saveRDS(poll.mod.1d, file = "Models/pollutionMod.rds")








poll2 <- poll[which(poll$GCDType %in% c("Metals", "Pesticides")),]

poll2.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll2)

anova(poll2.mod.1, btt = ":") #  not significant
summary(poll.mod.1)


## ### THINKING ABOUT OTHER COVARIATES



mod.gsba <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
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
  mods=~driver * GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=taxa_dat)
anova(mod.gsba.taxa, btt = ":") #  SIGNIFICANT

beep()
summary(mod.gsba.taxa)


saveRDS(mod.gsba.taxa, file = "Models/GSBAMod.rds")


t_dat <-predict(mod.gsba.taxa, newmods=rbind(c(0,0,0,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,0, 1,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,0, 0,1,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                              # intercept (climate change)
                                              
                                              
                                                c(1,0,0,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(1,0,0,0,0, 1,0,0,  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(1,0,0,0,0, 0,1,0,  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
                                                c(1,0,0,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
                                              #  HabitatLoss                                        
                                              
                                                c(0,1,0,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,1,0,0,0, 1,0,0,  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,1,0,0,0, 0,1,0,  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                                                c(0,1,0,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
                                              # Invasives 
                                              
                                                c(0,0,1,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,1,0,0, 1,0,0,  0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,1,0,0, 0,1,0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                                                c(0,0,1,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
                                              # LUI
                                              
                                                c(0,0,0,1,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,1,0, 1,0,0,  0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,1,0, 0,1,0,  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
                                                c(0,0,0,1,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
                                              # NutrientEnrichment
                                              
                                                c(0,0,0,0,1, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,1, 1,0,0,  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                                                c(0,0,0,0,1, 0,1,0,  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                                                c(0,0,0,0,1, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)),
                                              # Pollution 
                   addx=TRUE, digits=2) #

gcds <- rep(c("Climate", "habitat", "invasives", "lui", "nutrient", "pollution"), each = 4)
taxa <- rep(c("acari","collembola", "earthworms", "nematodes"), times = 6)

slabs <- paste(gcds, taxa)

# exclude habitat frag and invasives (little data)
subs <- c(1:4, 13:24)

par(mar=c(3, 8, 1, 1))
forest(t_dat$pred[subs], sei=t_dat$se[subs], slab=slabs[subs],  xlab="Effect Size", xlim=c(-.4,.7))


## ### THINKING ABOUT OTHER COVARIATES - habitat

system_dat <- hedges[which(hedges$System != ""),]

mod.system <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * System, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
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
  digits=4,
  data=system_dat)
anova(mod.system2, btt = "System") #  not significant
anova(mod.system2, btt = "driver") #   significant



saveRDS(mod.system2, file = "Models/SystemMod.rds")



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

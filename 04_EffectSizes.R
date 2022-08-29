

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


## FULL MODEL ---------




# If fittiung with REML (which is recommended) then can't compare models 
# (https://stats.stackexchange.com/questions/48671/what-is-restricted-maximum-likelihood-and-when-should-it-be-used)
# https://stats.stackexchange.com/questions/517857/testing-the-effect-of-moderators-in-metafor-package
# so test moderators using the anova(singleMod) way




write.csv(hedges, "Data/03_Data/HedgesData_cleaned.csv")
# hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")



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




funnel(mod.2, main="Standard Error")
funnel(mod.2, yaxis="vi", main="Sampling Variance")
funnel(mod.2, yaxis="seinv", main="Inverse Standard Error")
funnel(mod.2, yaxis="vinv", main="Inverse Sampling Variance")


## Using publication bias test from Nakagawa
hedges$sei <- sqrt(hedges$vi)

## Using the all in approach to check for effect of time

meta <- read.csv("Data/February2022/processed/metadata.csv")
meta$year <- sapply(strsplit(meta$NameOfPDF,'_'), "[", 2)
hedges <- merge(hedges, meta[,c('ID', 'year')], by = "ID", all.x = TRUE)
hedges$year <- as.integer(hedges$year)
hedges$year[which(is.na(hedges$year))] <- 2015 # just a study that hadn't been inputted fully
hedges$year.c <- as.vector(scale(hedges$year, scale = F))




# extracting the mean and 95% confidence intervals
# functions from here: https://github.com/elmacartney/EE_stress_MA/blob/7a46862ea5b015e7dc5d97d9f6952569135e8c8c/R/Old%20files/functions.R
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}
# custom function for extracting mean and CI for emmeans (marginalized means)
estimates.CI2 <- function(res){
  db.mf <- data.frame(summary(res)[,2],row.names = 1:length( summary(res)[,2]))
  db.mf <- cbind(db.mf,summary(res)[,5],summary(res)[,6],paste0(names(res@levels),summary(res)[,1] ))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}



# Application of Equation 21 from the main text
publication.bias.model.r.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + driver + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)


# print(publication.bias.model.r.se,digits=3)

# extracting the mean and 95% confidence intervals
estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.r.se)


## Indicates bias
## Can account for that bias in the estimates



publication.bias.model.time.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + driver + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.time.allin, data = hedges, at = list(var = 0, year.c = 0))
# marginalized overall mean at vi = 0 and year.c = 0; also weights = "prop" or "cells" average things over proportionally. if not specified, all groups (levels) get the same weights
#overall.res.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = ~1, df = 104 - 7, weights = "prop") # using effect size - 7 

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "driver")

# comparing with results without correcting for publication bias
# this model (without the intercept) are whats reported in the manuscript
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + driver, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

summary(publication.bias.model.r.v.1b)


estimates.publication.bias.model.r.v.1 <- estimates.CI2(mm.publication.bias.model.r.v.1)
estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)

table.comparing.captivity.levels <- merge(estimates.publication.bias.model.r.v.1,
                                          estimates.publication.bias.model.r.v.1b,
                                          by="estimate",
                                          all.x=T)

# rounding estimates
table.comparing.captivity.levels <- table.comparing.captivity.levels %>% mutate(across(where(is.numeric), round, 2))


table.comparing.captivity.levels <- data.frame(driver = table.comparing.captivity.levels[,1],
                                               adjusted.mean=table.comparing.captivity.levels[,2],               adjusted.CI=paste0("[",table.comparing.captivity.levels[,3],",",table.comparing.captivity.levels[,4],"]"),
                                               unadjusted.mean=table.comparing.captivity.levels[,5],                                               unadjusted.CI=paste0("[",table.comparing.captivity.levels[,6],",",table.comparing.captivity.levels[,7],"]"))













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


## climate pub bias

# recentering pub year
climate$year.c <- as.vector(scale(climate$year, scale = F))



publication.bias.model.climate.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + GCDType + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.climate.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.climate.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + GCDType + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.climate.allin, data = climate, at = list(var = 0, year.c = 0))

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "GCDType")

# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

# also use this model for coefficient sin the manuscript

estimates.publication.bias.model.r.v.1 <- estimates.CI2(mm.publication.bias.model.r.v.1)
estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)

table.comparing.captivity.levels <- merge(estimates.publication.bias.model.r.v.1,
                                          estimates.publication.bias.model.r.v.1b,
                                          by="estimate",
                                          all.x=T)

# rounding estimates
table.comparing.captivity.levels <- table.comparing.captivity.levels %>% mutate(across(where(is.numeric), round, 2))


table.comparing.captivity.levels <- data.frame(driver = table.comparing.captivity.levels[,1],
                                               adjusted.mean=table.comparing.captivity.levels[,2],               adjusted.CI=paste0("[",table.comparing.captivity.levels[,3],",",table.comparing.captivity.levels[,4],"]"),
                                               unadjusted.mean=table.comparing.captivity.levels[,5],                                               unadjusted.CI=paste0("[",table.comparing.captivity.levels[,6],",",table.comparing.captivity.levels[,7],"]"))





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


table(lui$GCDType, lui$GSBA) # to check for taxonomic model
# all four are fine



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

summary(lui.mod.2)





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



## lui pub bias
lui$year.c <- as.vector(scale(lui$year, scale = F))



publication.bias.model.lui.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + GCDType + Body.Size + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.lui.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.lui.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + GCDType + Body.Size + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.lui.allin, data = lui, at = list(var = 0, year.c = 0, Body.Size = "Macro-fauna"))

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "GCDType")

# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)

summary(publication.bias.model.r.v.1b)
# also using these coefficients in the manuscript

estimates.publication.bias.model.r.v.1 <- estimates.CI2(mm.publication.bias.model.r.v.1)
estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)

table.comparing.captivity.levels <- merge(estimates.publication.bias.model.r.v.1,
                                          estimates.publication.bias.model.r.v.1b,
                                          by="estimate",
                                          all.x=T)

# rounding estimates
table.comparing.captivity.levels <- table.comparing.captivity.levels %>% mutate(across(where(is.numeric), round, 2))


table.comparing.captivity.levels <- data.frame(driver = table.comparing.captivity.levels[,1],
                                               adjusted.mean=table.comparing.captivity.levels[,2],               adjusted.CI=paste0("[",table.comparing.captivity.levels[,3],",",table.comparing.captivity.levels[,4],"]"),
                                               unadjusted.mean=table.comparing.captivity.levels[,5],                                               unadjusted.CI=paste0("[",table.comparing.captivity.levels[,6],",",table.comparing.captivity.levels[,7],"]"))









## NUTRIENT ENRICHMENT -------



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


## nut enrich pub bias
nutri$year.c <- as.vector(scale(nutri$year, scale = F))



publication.bias.model.nutri.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + GCDType + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.nutri.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.nutri.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + GCDType + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.nutri.allin, data = nutri,
                                         at = list(var = 0, year.c = 0))

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "GCDType")

# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(publication.bias.model.r.v.1b) # for the manuscript 

estimates.publication.bias.model.r.v.1 <- estimates.CI2(mm.publication.bias.model.r.v.1)
estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)

table.comparing.captivity.levels <- merge(estimates.publication.bias.model.r.v.1,
                                          estimates.publication.bias.model.r.v.1b,
                                          by="estimate",
                                          all.x=T)

# rounding estimates
table.comparing.captivity.levels <- table.comparing.captivity.levels %>% mutate(across(where(is.numeric), round, 2))


table.comparing.captivity.levels <- data.frame(driver = table.comparing.captivity.levels[,1],
                                               adjusted.mean=table.comparing.captivity.levels[,2],               adjusted.CI=paste0("[",table.comparing.captivity.levels[,3],",",table.comparing.captivity.levels[,4],"]"),
                                               unadjusted.mean=table.comparing.captivity.levels[,5],                                               unadjusted.CI=paste0("[",table.comparing.captivity.levels[,6],",",table.comparing.captivity.levels[,7],"]"))







## Invasives

invas <- hedges[which(hedges$driver == "Invasives"),] # 188

## GCD type
table(invas$GCDType)
# Probably need to remove plants-mixture
invas <- invas[which(invas$GCDType != "Plants-mixture"),]



## Body size
table(invas$GCDType, invas$Body.Size)
table(invas$GCDType, invas$GSBA) # to check for taxonomic model
# as good as can be




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



## invasives pub bias
invas$year.c <- as.vector(scale(invas$year, scale = F))



publication.bias.model.invas.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.invas.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.invas.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)



# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~1, ## intercept only model 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)

summary(publication.bias.model.r.v.1b)

estimates.publication.bias.model.r.v.1 <- estimates.CI(publication.bias.model.invas.allin)

estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)
# Most previous code doesn't work, because its an intercept only model





### POLLUTION
poll <- hedges[which(hedges$driver == "Pollution"),] # 850


fulltext <- read.csv("Data/April2022/Full Text Screening - Sheet1.csv")
fulltext <- fulltext[,c('PaperID', 'PollutionSource')]

poll <- merge(poll, fulltext, by.x = "ID", by.y = "PaperID", all.x = TRUE)

# pollution type
table(poll$GCDType)

poll <- poll[which(poll$GCDType %in% c("Metals", "Pesticides")),] # the only two that have
# enough data/studies across the different pollutant types and body sizes



poll$Body.Size <- as.factor(poll$Body.Size)
poll$Body.Size <- relevel(poll$Body.Size, ref = "Meso-fauna")




# pollution type
table(poll$GCDType, poll$PollutionSource)




## previously missing
poll$PollutionSource[which(poll$ID == 105)] <- "Agricultural"
poll$PollutionSource[which(poll$ID == 2011)] <- "Agricultural"
poll$PollutionSource[which(poll$ID == 3372)] <- "Others"

## changing to remove multiple sources
poll$PollutionSource[which(poll$ID == 395)] <- "Mining/Smelting"
poll$PollutionSource[which(poll$ID == 2433 )] <- "Waste/sewage"
poll$PollutionSource[which(poll$ID == 1648)] <- "Mining/Smelting"
poll$PollutionSource[which(poll$ID == 1850)] <- "Mining/Smelting"
poll$PollutionSource[which(poll$ID == 564)] <- "Industrial"
poll$PollutionSource[which(poll$ID == 2509)] <- "Industrial"
poll$PollutionSource[which(poll$ID == 2715)] <- "Urban/transport"


poll$PollutionSource[which(poll$PollutionSource == "Agricultural")] <- "Farming"
poll$PollutionSource[which(poll$PollutionSource == "Agricultural/livestock")] <- "Farming"
poll$PollutionSource[which(poll$PollutionSource == "Military/wars")] <- "Others"
poll$PollutionSource[which(poll$PollutionSource == "Natural/geogenic")] <- "Others"



poll$GCDTypeSource <- paste(poll$GCDType, "-", poll$PollutionSource)
aggregate(poll$ID, by = list(GCD = poll$GCDTypeSource), function(x) length(unique(x)))
poll$GCDTypeSource[which(poll$GCDTypeSource == "Pesticides - Industrial")] <- "Pesticides - Others"

table(poll$GCDTypeSource, poll$Body.Size)


table(poll$GCDTypeSource, poll$GSBA) # to check for taxonomic model
# all four are not bad



# measurement
table(poll$GCDType, poll$Measurement)

# body size
table(poll$Body.Size)
table(poll$GCDType, poll$Body.Size)



poll <- droplevels(poll)

saveRDS(poll, file = "Models/pollutionDataFrame.rds")




poll.mod.10<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDTypeSource * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)


anova(poll.mod.10, btt = ":") #  not significant

poll.mod.11<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDTypeSource + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)


anova(poll.mod.11, btt = "GCD") #  nearly significant
anova(poll.mod.11, btt = "Size") #  not significant


poll.mod.12<-rma.mv(
  yi=effect,
  V=var, 
 mods=~GCDTypeSource, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)
anova(poll.mod.12, btt = "GCD") #  nearly significant

poll.mod.13<-rma.mv(
  yi=effect,
  V=var, 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

summary(poll.mod.13)

saveRDS(poll.mod.13, file = "Models/pollutionMod.rds")




## pollution pub bias
poll$year.c <- as.vector(scale(poll$year, scale = F))



publication.bias.model.poll.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.poll.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.poll.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)



# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~1, ## intercept only model 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)


estimates.publication.bias.model.r.v.1 <- estimates.CI(publication.bias.model.poll.allin)

estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)
# Most previous code doesn't work, because its an intercept only model




## Would the pollutant type make a difference if we accounted for bias?
publication.bias.model.poll.allin_fakemodel <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + var + GCDTypeSource, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

anova(publication.bias.model.poll.allin_fakemodel, btt = "GCD") #  nearly significant




## Looking at the raw values

plot(poll$effect ~ as.factor(poll$GCDTypeSource))
summary(poll$effect ~ as.factor(poll$GCDTypeSource))

tapply(poll$effect, as.factor(poll$GCDTypeSource), summary)


poll[which(poll$effect < -20),]
poll[which(poll$effect > 20),]

poll_cut <- poll[which(poll$effect > -20),]
poll_cut <- poll_cut[which(poll_cut$effect < 20),]
plot(poll_cut$effect ~ as.factor(poll_cut$GCDTypeSource))
stripchart(poll_cut$effect ~ as.factor(poll_cut$GCDTypeSource),
           vertical=TRUE,
           method = "jitter",
           pch = 19,
           col = 1:7,
           add = TRUE)


### HABITAT FRAG
frag <- hedges[which(hedges$driver == "HabitatLoss"),] # 100

table(frag$GCDType)
frag <- frag[which(frag$GCDType != "Fragmentation per se"),]

table(frag$Measurement)
table(frag$Measurement, frag$GCDType)


table(frag$Body.Size)
frag <- frag[which(frag$Body.Size != "All sizes"),]


table(frag$Body.Size, frag$GCDType) # just not enough data for an interaction
 


frag.mod1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=frag)

anova(frag.mod1, btt = "GCD") #  not significant
anova(frag.mod1, btt = "Size") #  not significant

frag.mod2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=frag)

anova(frag.mod2, btt = "GCD") #  not significant


frag.mod3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~1, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=frag)


summary(frag.mod3)

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
  mods=~- 1 + driver * GSBA, ## 
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



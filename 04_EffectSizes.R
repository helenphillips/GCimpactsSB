

## Script to calculate effect sizes

library(Hmisc)
library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished



# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")

## LOAD THE DATA
dataDir <- "Data/September2022"

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
## That's fine (8 abundance and 1 biomass)



# there are some climate cases that have no SDs
DataCleaned <- DataCleaned[which(!(is.na(DataCleaned$Control_SD))),] 
DataCleaned <- DataCleaned[which(!(is.na(DataCleaned$Treatment_SD))),] # 3440



## How many cases (and which drivers) now have zeros
nrow(DataCleaned[which(DataCleaned$Control_mean == 0),]) # 127
nrow(DataCleaned[which(DataCleaned$Control_SD   == 0),]) # 159

table(DataCleaned$driver[which(DataCleaned$Control_SD   == 0)]) # Mainly lui, nutrient enrichment and pollution


nrow(DataCleaned[which(DataCleaned$Treatment_mean  == 0),]) # 161
nrow(DataCleaned[which(DataCleaned$Treatment_SD  == 0),]) # 195
table(DataCleaned$driver[which(DataCleaned$Treatment_SD   == 0)]) # Mainly lui,  nutrient enrichment and pollution




hedges <- escalc(measure = "SMD", # 
                 m2i = Control_mean, # group 2 corresponds to the control group
                 sd2i = Control_SD,
                 n2i = Control_N,
                 m1i = Treatment_mean, # group 1 is the treatment group
                 sd1i = Treatment_SD,
                 n1i = Treatment_N,
                 data = DataCleaned)
## get a warning message here because of the zeros, 
# 3422

## 99 NAs in hedges
## when both the control AND the treatment is zero
# When theres no Ns, or when there are zeros, in the Sds?


hedges <- hedges[which(!(is.na(hedges$yi))),]
# 3323

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
## 3319


head(hedges[which(hedges$var > 2),])
table(hedges$Measurement[which(hedges$var > 10)])
head(hedges[which(hedges$var > 10),])




## Adding in variables that will be needed for publication bias analyses
## Using publication bias test from Nakagawa
hedges$sei <- sqrt(hedges$vi)

## Using the all in approach to check for effect of time

meta <- read.csv("Data/September2022/processed/metadata.csv")
meta$year <- sapply(strsplit(meta$NameOfPDF,'_'), "[", 2)
hedges <- merge(hedges, meta[,c('ID', 'year')], by = "ID", all.x = TRUE)
hedges$year <- as.integer(hedges$year)
hedges$year[which(is.na(hedges$year))] <- 2015 # just a study that hadn't been inputted fully
hedges$year.c <- as.vector(scale(hedges$year, scale = F))





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


## merge in others to species richness
hedges$Measurement[which(hedges$Measurement == "SpeciesTaxaRichness")] <- "Richness"


hedges$Body.Size[which(hedges$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
hedges$Body.Size[which(hedges$Body.Size %in% c("Macro-arthropods", "Macro-invertebrates"))] <- "Macro-fauna"
hedges$Body.Size[which(hedges$Body.Size %in% c("Micro-arthropods"))] <- "Meso-fauna"


# Make macrofauna the baseline

hedges$Body.Size <- as.factor(hedges$Body.Size)
hedges$Body.Size <- relevel(hedges$Body.Size, ref = "Macro-fauna")



## SAVE ------
# 3173

write.csv(hedges, "Data/03_Data/HedgesData_cleaned.csv")






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


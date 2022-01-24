# Script to get Metadata for Phil Martin
# He is doing a meta-analysis on disturbances on soil biodiversity
# I am sharing data


## LIBRARIES -----------------------------------------------

setwd("~/WORK/GCimpactsSB")

### LOAD LIBRARIES ------------------------

library(ggplot2)
library(metafor)
library(gtools)



### LOAD DATA -----------------------------

dataDir <- "Data/September2021/processed"
alldat <- read.csv(file.path(dataDir, "alldata.csv"))
meta <- read.csv(file.path(dataDir, "metadata.csv"))

### Filter data ----------------------------
alldat <- alldat[which(alldat$System == "Woody"),]

alldat <- alldat[,c("ï..ID","GCDType", "ChangeType", "TaxaGroup","TaxaBodysize","driver")]

alldat <- alldat[!duplicated(alldat), ]


## Bit of cleaning -------------------------
names(alldat)[which(names(alldat) == "ï..ID")] <- "ID"
alldat$TaxaGroup <- tolower(alldat$TaxaGroup)
alldat$TaxaGroup <- trimws(alldat$TaxaGroup, which = "both")



## Merge with the meta-data
alldat <- merge(alldat, meta, by = "ID", all.x = TRUE)

alldat$Sampled_Soil_Available_notpH <- NULL


write.csv(alldat, file.path(dataDir, "forPhilMartin.csv"), row.names = FALSE)

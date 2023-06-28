## A script to figure out some of the things that reviewers asked for


setwd("C:/Users/helenp/WORK/GCimpactsSB")
## LOAD THE DATA
dataDir <- "Data/September2022"
dat <- read.csv(file = file.path(dataDir, "processed", "alldata.csv"))



### PH --------------------

# First: Can we use pH in a model
# We did collect this, but never processed it. Not sure how much is available


## soilPh column is a character, as if ph was available for the control and 
## the treatment site, then both were collected, with a semi-colon in the middle
length(which(is.na(dat$SoilPH))) ## 77
# plus the empty character strings

length(which(dat$SoilPH == "")) # 1538
1538 + 77 # 1615

dat$UniqueID <- paste(dat$ID, dat$Case_ID, dat$driver)

ph <- dat[,c('UniqueID', 'SoilPH')]

# the final dataset
allDat <- read.csv("Data/03_Data/HedgesData_cleaned.csv")
head(allDat)

allDat2 <- merge(allDat, ph, by = 'UniqueID', all.x = TRUE)

length(which(is.na(allDat2$SoilPH))) ## 75
length(which(allDat2$SoilPH == "")) # 1376

75 + 1376 # 1451
table(allDat2$driver)
table(allDat2$driver[which(allDat2$SoilPH == "")])

# how many cases have ph
table(allDat2$driver) - table(allDat2$driver[which(allDat2$SoilPH == "")])



 
### LAB OR FIELD EXPERIMENT ---------------

meta <- read.csv("Data/February2022/processed/metadata.csv")
table(meta$SpatialExtent.km2.) # this is the number of papers




### COORDINATES -------------
# again, this is a character column, due to format of coordinates in dms
# also, some definitely seem to have country names in


length(which(dat$Latitude %in% c("Devecser, Hungary", 
"Kosogorsky near Tula, Russia",
"Shenyang Experimental Station of Ecology",
"Derio, Spain",
"Ghent (nearby)",
"Biesbosch", 
"Atcala de Henares",
"Naples",
"Krompachy, Slovakia",
"Donana Ntional Park, Spain",
"Mount Paddusstieva",                                                                                                         
"Mount Slattatjakka",
"Tianjin Normal University campus.",                                                                                                   
"Mortagne-du-Nord, France",
"Research Institute for Soil Science and Agricultural Chemistry (RISSAC) of the Hungarian Academy of Sciences at Nagyhorcsok, Hungary",
"Copper smelting works, Glogow, poland",                                                                                               
"Hygum site, Jutland",   
"germany")))   
# 179


length(which(dat$Latitude == "")) # 683
length(which(is.na(dat$Latitude))) # 110

## I think something can be done with these, so taking them into a new script
## (10) to work on
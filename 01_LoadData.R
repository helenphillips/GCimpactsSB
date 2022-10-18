## SetWD

## LIBRARIES -----------------------------------------------

# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")

### LOAD LIBRARIES ------------------------

library(ggplot2)
library(metafor)
library(gtools)

### FUNCTIONS -------

cleandata <- function(dat, IDCol){
  
  id <- which(names(dat) == IDCol)
  
  
  #remove blank rows
  dat <- dat[which(!(is.na(dat[,id]))),]
  dat <- dat[rowSums(is.na(dat)) != ncol(dat),]
  
  # remove blank columns
  dat[,grep("^X", names(dat))] <- NULL
  
  # Make sure the ID col is called 'ID'
  if(IDCol != "ID"){names(dat)[id] <- "ID"}
  return(dat)
}




### LOAD DATA -----------------------------

dataDir <- "Data/September2022"




## Metadata files first
meta <- read.csv(file.path(dataDir, "metadata.csv"))
meta <- cleandata(dat = meta, IDCol = "ID")

length(unique(meta$ID)) == nrow(meta) # true


## Main data files
frag <- read.csv(file.path(dataDir, "fragmentation.csv"))
frag <- cleandata(dat = frag, IDCol = "ID") # 117 # now only 111, where did they go? 
length(unique(paste(frag$ID, frag$Case_ID, sep = "_"))) == nrow(frag) # true



climate <- read.csv(file.path(dataDir, "climatechange.csv"))
climate <- cleandata(dat = climate, IDCol = "ID") # 507 # now 504, where did they go?
length(unique(paste(climate$ID, climate$Case_ID, sep = "_"))) == nrow(climate) # true





dat_inv <- read.csv(file.path(dataDir, "invasives.csv"))
invasives <- cleandata(dat = dat_inv, IDCol = "ID") # 146 # 188 now
length(unique(paste(invasives$ID, invasives$Case_ID, sep = "_"))) == nrow(invasives) # true




lui <- read.csv(file.path(dataDir, "lui.csv"))
lui <- cleandata(dat = lui, IDCol = "ID") # 801 # 991
length(unique(paste(lui$ID, lui$Case_ID, sep = "_"))) == nrow(lui) # true



nutrient <- read.csv(file.path(dataDir, "nutrient.csv"))
nutrient <- cleandata(dat = nutrient, IDCol = "ID") # 659 # 889
length(unique(paste(nutrient$ID, nutrient$Case_ID,  sep = "_"))) == nrow(nutrient) # true



dat_poll <- read.csv(file.path(dataDir, "pollution.csv"))
pollution <- cleandata(dat = dat_poll, IDCol = "ID") # 939 #926 , where did they go?
length(unique(paste(pollution$ID, pollution$Case_ID, sep = "_"))) == nrow(pollution) # true


## RENAME SOME COLUMNS ----------------------------------------------

names(pollution)[which(names(pollution) == "PollutionType")] <- "GCDType"
names(pollution)[which(names(pollution) == "PollutantClass")] <- "ChangeType"

names(climate)[which(names(climate) == "ClimateChangeVar")] <- "GCDType"
names(climate)[which(names(climate) == "Varied")] <- "ChangeType"
names(climate)[which(names(climate) == "Delta_Value")] <- "Delta_TreatmentChange"


names(lui)[which(names(lui) == "Intensification_Var")] <- "GCDType"
names(lui)[which(names(lui) == "Varied")] <- "ChangeType"

names(frag)[which(names(frag) == "FragmentationDesign")] <- "GCDType"
# names(HabitatLoss)[which(names(HabitatLoss) == "PollutantClass")] <- "ChangeType"

names(invasives)[which(names(invasives) == "InvasiveSpecies")] <- "GCDType"
# names(Invasives)[which(names(Invasives) == "PollutantClass")] <- "ChangeType"
names(invasives)[which(names(invasives) == "SpeciesName")] <- "InvasiveSpeciesName"


# names(nutrient)[which(names(nutrient) == "EnrichmentType")] <- "GCDType"
names(nutrient)[which(names(nutrient) == "Type")] <- "GCDType"


pollution$driver <- "Pollution"
climate$driver <- "Climate" 
lui$driver <- "LUI"
frag$driver <- "HabitatLoss" 
invasives$driver <- "Invasives"
nutrient$driver <- "NutrientEnrichment" 


## MERGE DATA SETS -----------------------------------------------------

## full text
fts <-  read.csv(file.path(dataDir, "Full Text Screening - Sheet1.csv"))
keep <- c("NameOfPDF","PaperID","DOI") #, "Screener","Extractor",
           # "DataExtracted","Extraction-Suitable","ExtractionComments",
           #"Method","ExperimentObservation")
fts <- fts[,names(fts) %in% keep]



## metadata
keep <- c("ID","Author", # "Title","Year",
          # "CodingFinished", "Checked_fullyfactorial",
          "Country",
          "Season.s.","SpatialExtent.km2.", "Sampled_Soil_Available_notpH") #, "FullyFactorial")                             

meta <- meta[,names(meta) %in% keep]

allMD <- merge(meta, fts, by.x = "ID", by.y = "PaperID", all.x = TRUE) #600
length(unique(allMD$ID)) == nrow(allMD)



# nrow pollution, climate, lui, frag, invasives, nutrients
924 + 508 + 966 + 115 + 188 + 896
## 3597


## The merge
dat <- smartbind(pollution, climate, lui, frag, invasives, nutrient)
nrow(dat) # 3597



## Bind with AllMd

#toKeep <- unique(dat$ID)
#allMD <- allMD[!(is.na(allMD$ID)),]
#allMD <- allMD[allMD$ID %in% toKeep,]
#dat <- merge(allMD, dat, by = "ID", all.y = TRUE)


write.csv(allMD, file = file.path(dataDir, "processed", "metadata.csv"), row.names = FALSE)
write.csv(dat, file = file.path(dataDir, "processed", "alldata.csv"), row.names = FALSE)



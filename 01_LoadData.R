## SetWD

## LIBRARIES -----------------------------------------------

setwd("~/WORK/GCimpactsSB")

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

dataDir <- "Data/September2021"




## Metadata files first
meta <- read.csv(file.path(dataDir, "metadata.csv"))
meta <- cleandata(dat = meta, IDCol = "ID")

length(unique(meta$ID)) == nrow(meta)


## Main data files
frag <- read.csv(file.path(dataDir, "fragmentation.csv"))
frag <- cleandata(dat = frag, IDCol = "ID") # 117

climate <- read.csv(file.path(dataDir, "climatechange.csv"))
climate <- cleandata(dat = climate, IDCol = "ID") # 507

dat_inv <- read.csv(file.path(dataDir, "invasives.csv"))
invasives <- cleandata(dat = dat_inv, IDCol = "ID") # 146

lui <- read.csv(file.path(dataDir, "lui.csv"))
lui <- cleandata(dat = lui, IDCol = "ID") # 801

nutrient <- read.csv(file.path(dataDir, "nutrient.csv"))
nutrient <- cleandata(dat = nutrient, IDCol = "ID") # 659

dat_poll <- read.csv(file.path(dataDir, "pollution.csv"))
pollution <- cleandata(dat = dat_poll, IDCol = "ID") # 939


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


names(nutrient)[which(names(nutrient) == "EnrichmentType")] <- "GCDType"
names(nutrient)[which(names(nutrient) == "Type")] <- "ChangeType"


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
          # "CodingFinished",
          "Country",
          "Season(s)","SpatialExtent(km2)", "Sampled_Soil_Available_notpH") #, "FullyFactorial")                             

meta <- meta[,names(meta) %in% keep]

allMD <- merge(meta, fts, by.x = "ID", by.y = "PaperID", all.x = TRUE)


# colOrder <- c("ID", "NameOfPDF", "Author", "DOI" , "CodingFinished", "Screener", "Extractor",
#               "DataExtracted", "Extraction-Suitable", "ExtractionComments",
#               "Sampled_Soil_Available_notpH", "FullyFactorial", "SpatialExtent(km2)",
#               "Country", "Season(s)", "ExperimentObservation" , "Method")
#                      
# allMD <- allMD[,colOrder]


## The merge
dat <- smartbind(pollution, climate, lui, frag, invasives, nutrient)

## Bind with AllMd

#toKeep <- unique(dat$ID)
#allMD <- allMD[!(is.na(allMD$ID)),]
#allMD <- allMD[allMD$ID %in% toKeep,]
#dat <- merge(allMD, dat, by = "ID", all.y = TRUE)


write.csv(allMD, file = file.path(dataDir, "processed", "metadata.csv"), row.names = FALSE)
write.csv(dat, file = file.path(dataDir, "processed", "alldata.csv"), row.names = FALSE)


## TODO:
# Script that ensures all columns are the right class
# Order of the columns
# enrichment - types need to be classified as organic or inorganic
# taxa groups need to be classified by body size
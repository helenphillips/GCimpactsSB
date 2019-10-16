## SetWD

## LIBRARIES -----------------------------------------------
library(googlesheets) # for getting data from Gdrive
library(gtools) # for merging dataframes



## CREATE FOLDERS ------------------------------------------------------
# One data folder
if(!dir.exists("Data")){
  dir.create("Data")
}
data_folder <- "Data"

# Data folders within it
n_Scripts <- 9 ## How many scripts are there

for(i in 1:n_Scripts){
  n <- (paste0("0", i, "_Data"))
  dir.create(file.path(data_folder, n))
}


## GET DATASETS ---------------------------------------------------------
##pull from the google
x <- gs_ls() ## Authentication
# gs_token <- gs_auth(cache = FALSE)

## Full text screening
fulltextscreening <- "Full Text Screening"
fts <- gs_title(fulltextscreening)
fts <- as.data.frame(gs_read(fts, ws = "Sheet1"))


## Metadata
dataEx <- "DataExtractionTable"
dataEx <- gs_title(dataEx)

metaData <- as.data.frame(gs_read(dataEx, ws = "Metadata"))
pollution <- as.data.frame(gs_read(dataEx, ws = "Pollution")) 
climate <- as.data.frame(gs_read(dataEx, ws = "ClimateChange")) 
LUI <- as.data.frame(gs_read(dataEx, ws = "LandUseIntensification")) 
HabitatLoss <- as.data.frame(gs_read(dataEx, ws = "HabitatLoss")) 
Invasives <- as.data.frame(gs_read(dataEx, ws = "Invasives")) 
NutrientEnrichment <- as.data.frame(gs_read(dataEx, ws = "NutrientEnrichment")) 

## RENAME SOME COLUMNS ----------------------------------------------

names(pollution)[which(names(pollution) == "PollutionType")] <- "GCDType"
names(pollution)[which(names(pollution) == "PollutantClass")] <- "ChangeType"

names(climate)[which(names(climate) == "ClimateChangeVar")] <- "GCDType"
names(climate)[which(names(climate) == "Varied")] <- "ChangeType"
names(climate)[which(names(climate) == "Delta_Value")] <- "Delta_TreatmentChange"


names(LUI)[which(names(LUI) == "Intensification_Var")] <- "GCDType"
names(LUI)[which(names(LUI) == "Varied")] <- "ChangeType"

names(HabitatLoss)[which(names(HabitatLoss) == "FragmentationDesign")] <- "GCDType"
# names(HabitatLoss)[which(names(HabitatLoss) == "PollutantClass")] <- "ChangeType"

names(Invasives)[which(names(Invasives) == "InvasiveSpecies")] <- "GCDType"
# names(Invasives)[which(names(Invasives) == "PollutantClass")] <- "ChangeType"
names(Invasives)[which(names(Invasives) == "SpeciesName")] <- "InvasiveSpeciesName"




names(NutrientEnrichment)[which(names(NutrientEnrichment) == "Type")] <- "GCDType"
names(NutrientEnrichment)[which(names(NutrientEnrichment) == "EnrichmentType")] <- "ChangeType"


pollution$driver <- "Pollution"
climate$driver <- "Climate" 
LUI$driver <- "LUI"
HabitatLoss$driver <- "HabitatLoss" 
Invasives$driver <- "Invasives"
NutrientEnrichment$driver <- "NutrientEnrichment" 


## MERGE DATA SETS -----------------------------------------------------

## full text
keep <- c("NameOfPDF","PaperID","DOI", "Screener","Extractor",
          "DataExtracted","Extraction-Suitable","ExtractionComments",
          "Method","ExperimentObservation")
fts <- fts[,names(fts) %in% keep]



## metadata
keep <- c("ID","Author", # "Title","Year",
          "CodingFinished", "Country",
          "Season(s)","SpatialExtent(km2)", "Sampled_Soil_Available_notpH", "FullyFactorial")                             

metaData <- metaData[,names(metaData) %in% keep]

allMD <- merge(metaData, fts, by.x = "ID", by.y = "PaperID")


colOrder <- c("ID", "NameOfPDF", "Author", "DOI" , "CodingFinished", "Screener", "Extractor",
              "DataExtracted", "Extraction-Suitable", "ExtractionComments",
              "Sampled_Soil_Available_notpH", "FullyFactorial", "SpatialExtent(km2)",
              "Country", "Season(s)", "ExperimentObservation" , "Method")
                     
allMD <- allMD[,colOrder]


## The merge
dat <- smartbind(pollution, climate, LUI, HabitatLoss, Invasives, NutrientEnrichment)

## Bind with AllMd

toKeep <- unique(dat$ID)
allMD <- allMD[!(is.na(allMD$ID)),]
allMD <- allMD[allMD$ID %in% toKeep,]


dat <- merge(allMD, dat, by = "ID", all.y = TRUE)

## TODO:
# Script that ensures all columns are the right class
# Order of the columns
# enrichment - types need to be classified as organic or inorganic
# taxa groups need to be classified by body size
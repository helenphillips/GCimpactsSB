## Clean the data

# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")



## LOAD THE DATA
dataDir <- "Data/February2022"

dat <- read.csv(file = file.path(dataDir, "processed", "alldata.csv"))



## CLEAN THE COLUMNS
tokeep <- c("ID","Case_ID","CaseDescription",
"GCDType","ChangeType",
  "TaxaGroup","TaxaBodysize",       
"Control_mean","Control_SD","Control_N"     ,         
"Treatment_mean","Treatment_SD","Treatment_N" ,           
"Measurement","MeasurementUnits","Error",                  
"Data_Source","driver")

dat <- dat[,which(names(dat) %in% tokeep)] # 3138 # 3609



###### CONVERSION OF COLUMN TYPES ----------------
dat$Control_SD_2 <- as.numeric(dat$Control_SD)
toquery <- which(is.na(dat$Control_SD_2) & !(is.na(dat$Control_SD)))
dat[toquery,]

dat$Treatment_SD_2 <- as.numeric(dat$Treatment_SD )
toquery <- which(is.na(dat$Treatment_SD_2) & !(is.na(dat$Treatment_SD)))
dat[toquery,]



dat$Treatment_SD_2 <- dat$Control_SD_2 <- NULL
dat$Control_SD <- as.numeric(dat$Control_SD)
dat$Treatment_SD <- as.numeric(dat$Treatment_SD )


###### CONVERSION OF ERROR -----------------------
# Should all be SDs

table(dat$Error)

## SEs
ses <- which(dat$Error == "SE")

dat$Control_SD[ses] <- dat$Control_SD[ses] * sqrt(dat$Control_N[ses])
dat$Treatment_SD[ses] <- dat$Treatment_SD[ses] * sqrt(dat$Treatment_N[ses])
dat$Error[ses] <- "SD"

## CI
## https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm#:~:text=The%20standard%20deviation%20for%20each,should%20be%20replaced%20by%205.15.
ci <- which(dat$Error == "CI")

dat$Control_SD[ci] <- sqrt(dat$Control_N[ci]) * ((2 * dat$Control_SD[ci])/3.92)
dat$Treatment_SD[ci] <- sqrt(dat$Treatment_N[ci]) * ((2 * dat$Treatment_SD[ci])/3.92)
dat$Error[ci] <- "SD"


ci <- which(dat$Error == "CI95")

dat$Control_SD[ci] <- sqrt(dat$Control_N[ci]) * ((2 * dat$Control_SD[ci])/3.92)
dat$Treatment_SD[ci] <- sqrt(dat$Treatment_N[ci]) * ((2 * dat$Treatment_SD[ci])/3.92)
dat$Error[ci] <- "SD"

ci <- which(dat$Error == "95CI")

dat$Control_SD[ci] <- sqrt(dat$Control_N[ci]) * ((2 * dat$Control_SD[ci])/3.92)
dat$Treatment_SD[ci] <- sqrt(dat$Treatment_N[ci]) * ((2 * dat$Treatment_SD[ci])/3.92)
dat$Error[ci] <- "SD"

ci <- which(dat$Error == "CI90")

dat$Control_SD[ci] <- sqrt(dat$Control_N[ci]) * ((2 * dat$Control_SD[ci])/3.29)
dat$Treatment_SD[ci] <- sqrt(dat$Treatment_N[ci]) * ((2 * dat$Treatment_SD[ci])/3.29)
dat$Error[ci] <- "SD"



## Convert CV to SD
cv <- which(dat$Error == "CV")

dat$Control_SD[cv] <- dat$Control_SD[cv] * dat$Control_mean[cv]
dat$Treatment_SD[cv] <- dat$Treatment_SD[cv] * dat$Treatment_mean[cv] # times the CV by the mean
dat$Error[cv] <- "SD"


## Get rid of DT
dat <- dat[-which(dat$Error == "DT"),]
dat <- droplevels(dat)

table(dat$Error)


## unknowns marked as SD (the largest error)
others <- which(dat$Error != "SD")
dat$Error[others] <- "SD"


## Keep blank if no SDs
dat$Error[which(is.na(dat$Treatment_SD))] <- NA



## 2) Clean values -------------------------------------------------------------

# GlobalChangeType
table(dat$GCDType)
# tioo difficult to fix at this stage 
dat$GCDType[which(dat$GCDType == "Organic versus conventional")] <-  "Organic versus Inorganic"
dat$GCDType[which(dat$GCDType == "Organic vs. Conventional")] <-  "Organic versus Inorganic"

dat$GCDType[which(dat$GCDType == "UVB-Radiation")] <-  "UVB Radiation"
dat$GCDType[which(dat$GCDType == "UVB")] <-  "UVB Radiation"

dat$GCDType[which(dat$GCDType == "Burning")] <-  "Fire"
dat$GCDType[which(dat$GCDType == "burning")] <-  "Fire"


## Move all Fire to LUI
dat$driver[which(dat$GCDType == "Fire")] <-  "LUI"


# BodySize and Taxonomic Group

# unique(dat$TaxaGroup[order(dat$TaxaGroup)])

taxa <- read.csv(file.path("Data","February2022", "taxonomic classification - version 3.csv"))
dat <- merge(dat, taxa, by.x = "TaxaGroup", by.y = "original_v2", all.x = TRUE)

unique(dat$TaxaGroup[which(!(is.na(dat$TaxaGroup)) & is.na(dat$Harmonised))])

dat[which(dat$TaxaGroup == "staphylinidae"),'Harmonised'] <- 'Rove beetle'
dat[which(dat$TaxaGroup == "staphylinidae"),'GSBA'] <- 'Coleoptera'
dat[which(dat$TaxaGroup == "staphylinidae"), 'Body.Size'] <- 'Macro-fauna'


dat[which(dat$TaxaGroup == "Enchytraeids "),'Harmonised'] <- 'Enchytraeids'
dat[which(dat$TaxaGroup == "Enchytraeids "),'GSBA'] <- 'Enchytraeids'
dat[which(dat$TaxaGroup == "Enchytraeids "), 'Body.Size'] <- 'Meso-fauna'

dat[which(dat$TaxaGroup == "Geophilida"),'Harmonised'] <- 'Chilopoda'
dat[which(dat$TaxaGroup == "Geophilida"),'GSBA'] <- 'Myriapoda'
dat[which(dat$TaxaGroup == "Geophilida"), 'Body.Size'] <- 'Macro-fauna'



## Add in body size when it is elsewhere


dat$Body.Size[which(!(is.na(dat$TaxaBodysize )) & is.na(dat$Body.Size))] <- dat$TaxaBodysize[which(!(is.na(dat$TaxaBodysize )) & is.na(dat$Body.Size))]



## Fix the column\
table(dat$Body.Size)

dat$Body.Size[which(dat$Body.Size %in% c("Macro-arthropod predators", "Macro-arthropod", "Macroarthropods"))] <-  "Macro-arthropods"
dat$Body.Size[which(dat$Body.Size %in% c("Microfauna"))] <-  "Micro-fauna"

dat$Body.Size[which(dat$Body.Size %in% c("Macrofauna"))] <-  "Macro-fauna"
dat$Body.Size[which(dat$Body.Size %in% c("Macroinvertebrates"))] <-  "Macro-invertebrates"
dat$Body.Size[which(dat$Body.Size %in% c("Mesofauna"))] <-  "Meso-fauna"

dat$Body.Size[which(dat$Body.Size %in% c("Meso and Macrofauna"))] <-  "Invertebrates (all sizes)"


# DiversityMetric


table(dat$Measurement)


dat$Measurement[which(dat$Measurement %in% c("FamilyRichness"    ,  
"GeneraRichness"           ,
"GenusRichness"         ,
"GroupRichness"      ,          
"Margalef" ,
"Margalef Richness index",
"MargalefRichness" ,               
"Richness"           ,      
"Species/Genus Richness",
"SpeciesRichness" ,
"TaxaRichness"    ,      
"Taxon richness",
"TaxonRichness"))] <- "Richness" # 325






dat$Measurement[which(dat$Measurement %in% c("Eveness"))] <- "Evenness"
       
                                             
dat$Measurement[which(dat$Measurement %in% c("Shannon Equitability"))] <- "Shannon"
                                                                                          

dat$Measurement[which(dat$Measurement %in% c("Simpson Equitability",  
                                             "Simpsons"))] <- "Simpson"
                                             



## SAVING --------
write.csv(dat, file = file.path(dataDir, "processed", "0_2_alldata.csv"), row.names = FALSE)



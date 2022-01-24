## Clean the data

setwd("~/WORK/GCimpactsSB")



## LOAD THE DATA
dataDir <- "Data/September2021"

dat <- read.csv(file = file.path(dataDir, "processed", "alldata.csv"))



## CLEAN THE COLUMNS
tokeep <- c("ID","Case_ID","CaseDescription",
"GCDType","ChangeType",
  "TaxaGroup","TaxaBodysize",       
"Control_mean","Control_SD","Control_N"     ,         
"Treatment_mean","Treatment_SD","Treatment_N" ,           
"Measurement","MeasurementUnits","Error",                  
"Data_Source","driver")

dat <- dat[,which(names(dat) %in% tokeep)] # 3138


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



## unknowns marked as SD (the largest error)
others <- which(dat$Error != "SD")
dat$Error[others] <- "SD"


## Keep blank if no SDs
dat$Error[which(is.na(dat$Treatment_SD))] <- NA

write.csv(dat, file = file.path(dataDir, "processed", "0_2_alldata.csv"), row.names = FALSE)


## 2) Clean values -------------------------------------------------------------

# GlobalChangeType

# ChangeType (amount/frequency/intensity)

# BodySize

# Taxonomic Group

# DiversityMetric

# System

# Location

# StudyType


# Study

# Case


## 2) Harmonizing the units-----------------------------------------------

# Abundance

# Richness

# Create mean soil pH

# Global Change Intensity (absolute values)





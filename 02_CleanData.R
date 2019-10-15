## Clean the data


# Load data from 01 folder
# CleanData <- read.csv("Data/01_LoadData/data.csv", h = TRUE) # change names to actual dataset name


## 1) Transform standard errors and confidence intervals to SD-----------

# functions to calculate SD from CI
myfun_SDfromCI95 <- function(CI, samplesize){
  sqrt(samplesize) * (CI * 2) / qt(c(.975), df=samplesize-1)
}

myfun_SDfromCI90 <- function(CI, samplesize){
  sqrt(samplesize) * (CI * 2) / qt(c(.95), df=samplesize-1)
}

# save initial error values
CleanData$Control_SDinit <- CleanData$Control_SD

# transform SE and CI into SD
CleanData$Control_SD <- ifelse(CleanData$Error == "SD", CleanData$Control_SDinit, 
                               # if SD is SD, do nothing
                               
                               ifelse(CleanData$Error == "SE", CleanData$Control_SDinit * sqrt(Control_n), 
                                      # if SE, calculate SD
                                      
                                      ifelse(CleanData$Error == "CI95", myfun_SDfromCI95(CI = CleanData$Control_SDinit,
                                                                                         samplesize = CleanData$Control_N),
                                            # If CI, use the custom functions depending on the limits (95% or 90% CI)
                                            
                                             ifelse(CleanData$Error == "CI90", myfun_SDfromCI90(CI = CleanData$Control_SDinit, 
                                                                                                samplesize = CleanData$Control_N), 
                                                    NA))))



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





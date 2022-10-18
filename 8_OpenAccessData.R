## Cleaning the data to make open access


setwd("C:/Users/helenp/WORK/GCimpactsSB")


hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")



## remove redundant or bad columns
hedges$X.1 <- NULL
hedges$X <- NULL
hedges$MeasurementUnits <- NULL
hedges$ChangeType <- NULL

hedges$TaxaGroup <- NULL
hedges$CaseDescription <- NULL
hedges$TaxaBodysize <- NULL
# 3173


## add in paper metadata

meta <- read.csv("Data/September2022/metadata.csv")
colstoKeep <- c("ID", "Author", "Title")
meta <- meta[,colstoKeep]
# 932

allDat <- merge(hedges, meta, by = "ID", all.x = TRUE)
# 3173

## add in DOI from elsewhere

fulltext <- read.csv("Data/September2022/Full Text Screening - Sheet1.csv")
fulltext <- fulltext[,c('PaperID', 'DOI')]

allDat <- merge(allDat, fulltext, by.x = "ID", by.y = "PaperID", all.x = TRUE )
# 3173


#### re-order the columns
NewOrder <- c("ID","Author","Title","year",
              "DOI","Case_ID","driver","UniqueID",
              "GCDType","System","Harmonised","GSBA",
              "Body.Size", "Control_mean","Control_SD","Control_N",
              "Treatment_mean", "Treatment_SD","Treatment_N",
              "Measurement","Error", "yi",  "vi","effect","var","sei",
              "year.c", "Data_Source")


allDat <- allDat[NewOrder]



## Save

saveFolder <- "Data/OpenAccessData"

write.csv(allDat, file = file.path(saveFolder, "PhillipsetalData.csv"), row.names = FALSE)        


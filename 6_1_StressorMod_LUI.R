
## Script to run main GC model

library(Hmisc)
library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished
library(emmeans)
library(dplyr)

setwd("C:/Users/helenp/WORK/GCimpactsSB")




# extracting the mean and 95% confidence intervals
# functions from here: https://github.com/elmacartney/EE_stress_MA/blob/7a46862ea5b015e7dc5d97d9f6952569135e8c8c/R/Old%20files/functions.R
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}
# custom function for extracting mean and CI for emmeans (marginalized means)
estimates.CI2 <- function(res){
  db.mf <- data.frame(summary(res)[,2],row.names = 1:length( summary(res)[,2]))
  db.mf <- cbind(db.mf,summary(res)[,5],summary(res)[,6],paste0(names(res@levels),summary(res)[,1] ))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}



## LUI -------
hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")

lui <- hedges[which(hedges$driver == "LUI"),] 


lui <- droplevels(lui[which(lui$Body.Size != "All sizes"),])



table(lui$GCDType)
lui$GCDType[which(lui$GCDType %in% c("Defoliation"))] <- "Grazing"
lui$GCDType[which(lui$GCDType %in% c("Intensity", "Landscape", "Management", "Weeding", "Planting", "Intensification"))] <- "Management"
lui$GCDType[which(lui$GCDType %in% c("Logging"))] <- "Harvesting"
lui$GCDType[which(lui$GCDType %in% c("Water"))] <- "Irrigation"
lui$GCDType[which(lui$GCDType %in% c("degradation", "Disturbance"))] <- "Degradation"

lui <- lui[which(lui$GCDType != "Human population"),]
lui <- lui[which(lui$GCDType != "Mono- versus poly-culture"),] # too little
lui <- lui[which(lui$GCDType != "Irrigation"),] # too little

lui <- lui[which(lui$GCDType != "Management"),] # too little
lui <- lui[which(lui$GCDType != "Degradation"),] # too little



table(lui$GCDType, lui$Measurement)


lui$GCDType <- relevel(as.factor(lui$GCDType), ref = "Grazing")



## Checking body size
table(lui$Body.Size)


table(lui$GCDType, lui$GSBA) # to check for taxonomic model
# all four are fine



lui.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)


anova(lui.mod.1, btt = ":") # not significant





lui.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)



anova(lui.mod.2, btt = "GCD") # significant
anova(lui.mod.2, btt = "Size") #  significant

summary(lui.mod.2)





## final model
saveRDS(lui.mod.2, file = "Models/LUIMod_redo.rds")





## lui pub bias
lui$year.c <- as.vector(scale(lui$year, scale = F))



publication.bias.model.lui.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + GCDType + Body.Size + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.lui.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.lui.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + GCDType + Body.Size + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.lui.allin, data = lui, at = list(var = 0, year.c = 0, Body.Size = "Macro-fauna"))

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "GCDType")

# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=lui)

summary(publication.bias.model.r.v.1b)
# also using these coefficients in the manuscript

estimates.publication.bias.model.r.v.1 <- estimates.CI2(mm.publication.bias.model.r.v.1)
estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)

table.comparing.captivity.levels <- merge(estimates.publication.bias.model.r.v.1,
                                          estimates.publication.bias.model.r.v.1b,
                                          by="estimate",
                                          all.x=T)

# rounding estimates
table.comparing.captivity.levels <- table.comparing.captivity.levels %>% mutate(across(where(is.numeric), round, 2))


table.comparing.captivity.levels <- data.frame(driver = table.comparing.captivity.levels[,1],
                                               adjusted.mean=table.comparing.captivity.levels[,2],               
                                               adjusted.CI=paste0("[",table.comparing.captivity.levels[,3],",",table.comparing.captivity.levels[,4],"]"),
                                               unadjusted.mean=table.comparing.captivity.levels[,5],                                               
                                               unadjusted.CI=paste0("[",table.comparing.captivity.levels[,6],",",table.comparing.captivity.levels[,7],"]"))




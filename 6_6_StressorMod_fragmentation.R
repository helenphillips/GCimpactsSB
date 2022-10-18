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




### HABITAT FRAG
hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")

frag <- hedges[which(hedges$driver == "HabitatLoss"),] # 100

table(frag$GCDType)
frag <- frag[which(frag$GCDType != "Fragmentation per se"),]

table(frag$Measurement)
table(frag$Measurement, frag$GCDType)


table(frag$Body.Size)
frag <- frag[which(frag$Body.Size != "All sizes"),]


table(frag$Body.Size, frag$GCDType) # just not enough data for an interaction



frag.mod1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=frag)

anova(frag.mod1, btt = "GCD") #  not significant
anova(frag.mod1, btt = "Size") #  not significant

frag.mod2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=frag)

anova(frag.mod2, btt = "GCD") #  not significant


frag.mod3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~1, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=frag)


summary(frag.mod3)


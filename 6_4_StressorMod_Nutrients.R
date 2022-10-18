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




## NUTRIENT ENRICHMENT -------
hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")


nutri <- hedges[which(hedges$driver == "NutrientEnrichment"),] # 820

table(nutri$GCDType)


table(nutri$GCDType, nutri$Measurement)


# Removing some GCDs and relevelling
nutri <- nutri[which(nutri$GCDType != "Biochar"),]
nutri$GCDType <- as.factor(nutri$GCDType)
nutri$GCDType <- relevel(nutri$GCDType, ref = "Synthetic Fertilizers")

table(nutri$GCDType, nutri$GSBA)

## Check body size
table(nutri$Body.Size)

table(nutri$GCDType, nutri$GSBA) # to check for taxonomic model
# acari, earthworms and nematodes are definitely fine. collembola could be better, but not bad

nutri.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(nutri.mod.1)


anova(nutri.mod.1, btt = ":") # not significant


nutri.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

anova(nutri.mod.1b, btt = "GCD") #  significant
anova(nutri.mod.1b, btt = "Size") # not significant


nutri.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(nutri.mod.2)
anova(nutri.mod.2, btt = "GCD") #  significant


saveRDS(nutri.mod.2, file = "Models/nutriMod.rds")


## nut enrich pub bias
nutri$year.c <- as.vector(scale(nutri$year, scale = F))



publication.bias.model.nutri.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + GCDType + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.nutri.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.nutri.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + GCDType + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.nutri.allin, data = nutri,
                                         at = list(var = 0, year.c = 0))

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "GCDType")

# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + GCDType, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=nutri)

summary(publication.bias.model.r.v.1b) # for the manuscript 

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





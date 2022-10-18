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




## Climate change

hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")


climate <- hedges[which(hedges$driver == "Climate"),] # 462

table(climate$GCDType)
table(climate$GCDType, climate$Measurement)


## going to remove vegetation as a gcd type. 2 lines are thickness of moss, and the 10 other lines are canopy gaps due to ice storm
climate <- climate[which(climate$GCDType != "Vegetation"),]

# For now will also remove the precip+temp
climate <- climate[which(climate$GCDType != "Precipitation+Temp"),]
climate <- climate[which(climate$GCDType != "UVB Radiation"),]

# Checking body size
table(climate$Body.Size)
# removing the all sizes
climate <- droplevels(climate[which(climate$Body.Size != "All sizes"),])


table(climate$GCDType, climate$GSBA) # to check for taxonomic model
# Acari, collembola and nematodes are definitrely fine. Earthworms less so

climate.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

summary(climate.mod.1)

anova(climate.mod.1, btt = ":") # should be testing the levels that are the interactions
# Not significant


climate.mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

anova(climate.mod.1b, btt = "GCD") # significant
anova(climate.mod.1b, btt = "Size") # not significant


climate.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="ML",
  digits=4,
  data=climate)
anova(climate.mod.5, btt = "GCD") # maybe not significant


## Use climate.mod.5
saveRDS(climate.mod.5, file = "Models/ClimateMod.rds")


qqnorm(residuals(climate.mod.5,type="pearson"),main="QQ plot: residuals")
qqline(residuals(climate.mod.5,type="pearson"),col="red") # fine




## climate pub bias

# recentering pub year
climate$year.c <- as.vector(scale(climate$year, scale = F))



publication.bias.model.climate.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + GCDType + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.climate.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.climate.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + GCDType + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=climate)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.climate.allin, data = climate, at = list(var = 0, year.c = 0))

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
  data=climate)

# also use this model for coefficient sin the manuscript

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






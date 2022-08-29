
## Script to run main GC model

library(Hmisc)
library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished

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

lui <- hedges[which(hedges$driver == "LUI"),] # 911


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
saveRDS(lui.mod.2, file = "Models/LUIMod.rds")




luidat <-predict(lui.mod.2, newmods=rbind(c(0,0,0,0,0,0),  # intercept
                                          c(0,0,0,0,1,0),
                                          c(0,0,0,0,0,1),
                                          
                                          c(1,0,0,0,0,0),# fire
                                          c(1,0,0,0,1,0),
                                          c(1,0,0,0,0,1),
                                          
                                          c(0,1,0,0,0,0), 
                                          c(0,1,0,0,1,0),
                                          c(0,1,0,0,0,1),# Harvesting   
                                          
                                          c(0,0,1,0,0,0),
                                          c(0,0,1,0,1,0),
                                          c(0,0,1,0,0,1),# organic
                                          
                                          c(0,0,0,1,0,0),
                                          c(0,0,0,1,1,0),
                                          c(0,0,0,1,0,1)# tillage
), addx=TRUE, digits=2) #

macro <- c(1, 4, 7, 10, 13) #3 just the macrofauna


slabs <- c("Grazing", "Fire", "Harvesting", "Inorganic", "Tillage")
par(mar=c(3, 8, 1, 1))
forest(luidat$pred[macro], sei=luidat$se[macro], slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))



#3 Just showing the impact of body size

measurement_coefs <- data.frame(meas = c("Macro-fauna", "Meso-fauna", "Micro-fauna"),
                                coef = c(0, lui.mod.2$beta[6:7]), 
                                ses = c(0, lui.mod.2$se[6:7]))


measurement_coefs <- measurement_coefs[c(3,2,1),]

errbar(x =measurement_coefs$meas, y = measurement_coefs$coef, 
       yplus = measurement_coefs$coef + measurement_coefs$ses,
       yminus=measurement_coefs$coef-measurement_coefs$ses, cap=0.015)
abline(v=0, lty =2)



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




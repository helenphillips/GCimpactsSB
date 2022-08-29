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





## Invasives

invas <- hedges[which(hedges$driver == "Invasives"),] # 188

## GCD type
table(invas$GCDType)
# Probably need to remove plants-mixture
invas <- invas[which(invas$GCDType != "Plants-mixture"),]



## Body size
table(invas$GCDType, invas$Body.Size)
table(invas$GCDType, invas$GSBA) # to check for taxonomic model
# as good as can be




invas.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType*Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)
anova(invas.mod.1, btt = ":") # not  significant


invas.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)



anova(invas.mod.2, btt = "GCD") # not  significant
anova(invas.mod.2, btt = "Size") # not  significant




invas.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  #  mods=~GCDType + Body.Size , ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)

summary(invas.mod.3) # not significant intercapt
saveRDS(invas.mod.3, file = "Models/invasiveMod.rds")



## invasives pub bias
invas$year.c <- as.vector(scale(invas$year, scale = F))



publication.bias.model.invas.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.invas.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.invas.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)



# comparing with results without correcting for publication bias
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~1, ## intercept only model 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=invas)

summary(publication.bias.model.r.v.1b)

estimates.publication.bias.model.r.v.1 <- estimates.CI(publication.bias.model.invas.allin)

estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)
# Most previous code doesn't work, because its an intercept only model


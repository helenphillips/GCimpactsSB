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






## ### THINKING ABOUT OTHER COVARIATES
hedges <- read.csv("Data/03_Data/HedgesData_cleaned_June2023.csv")



# mod.gsba <-rma.mv(
#   yi=effect,
#   V=var, 
#   mods=~driver + GSBA, ## 
#   random= list(~1|ID/UniqueID, 
#                ~ 1 | Measurement),
#   struct="CS",
#   method="REML",
#   digits=4,
#   data=hedges)
# 
# beep()
# 
# summary(mod.gsba)
# anova(mod.gsba, btt = "GSBA") # NOT SIGNIFICANT



# just large groups

taxa_dat <- hedges[hedges$GSBA %in% c("Acari", "Collembola",  "Earthworms", "Nematodes"),]
table(taxa_dat$driver, taxa_dat$GSBA)


## Reviewer commetns: re. using stressors rather than the main drivers


table(taxa_dat$GCDType[which(taxa_dat$driver == "Climate")], taxa_dat$GSBA[which(taxa_dat$driver == "Climate")])
table(taxa_dat$GCDType[which(taxa_dat$driver == "LUI")], taxa_dat$GSBA[which(taxa_dat$driver == "LUI")])
table(taxa_dat$GCDType[which(taxa_dat$driver == "NutrientEnrichment")], taxa_dat$GSBA[which(taxa_dat$driver == "NutrientEnrichment")])
table(taxa_dat$GCDType[which(taxa_dat$driver == "Pollution")], taxa_dat$GSBA[which(taxa_dat$driver == "Pollution")])


taxa_dat2 <- taxa_dat[which(taxa_dat$GCDType %in% c(
# climate
"Gas - CO2",
"Temperature",
"WaterAvailability-Drought",
"WaterAvailability-Flood",
# lui
"Fire",
"Grazing",
"Harvesting",
"Organic versus Inorganic",
"Tillage",
# nutrient
"Compost",
"Manure + Slurry",
"Other Organic fertilisers (NOT including compost and Urea)",
"Residue + Mulch",
"Sludge (including Biosolids)",
"Synthetic Fertilizers",
# pollution
"Metals",
"Pesticides")),]


table(taxa_dat2$GSBA, taxa_dat2$GCDType)
mod.gsba.taxa_stressor <-rma.mv(
  yi=effect,
  V=var, 
  mods=~ GCDType * GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=taxa_dat2)
anova(mod.gsba.taxa_stressor, btt = ":") #  SIGNIFICANT

summary(mod.gsba.taxa_stressor)

saveRDS(mod.gsba.taxa_stressor, file = "Models/GSBAMod_stressors.rds")




## The original model
mod.gsba.taxa <-rma.mv(
  yi=effect,
  V=var, 
  mods=~ driver * GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=taxa_dat)
anova(mod.gsba.taxa, btt = ":") #  SIGNIFICANT

beep()
summary(mod.gsba.taxa)


saveRDS(mod.gsba.taxa, file = "Models/GSBAMod_june2023.rds")








mod.gsba.taxa_nointercept <-rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + driver * GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=taxa_dat)
# for the manuscript

beep()


## ### THINKING ABOUT OTHER COVARIATES - habitat

system_dat <- hedges[which(hedges$System != ""),]

mod.system <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * System, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=system_dat)
anova(mod.system, btt = ":") #  not significant



mod.system2 <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + System, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=system_dat)
anova(mod.system2, btt = "System") #  not significant
anova(mod.system2, btt = "driver") #   significant



saveRDS(mod.system2, file = "Models/SystemMod.rds")


## kH method


## ph
phdat <- hedges[which(!(is.na(hedges$ph_fixed))),] # 1715
phdat <- phdat[which(phdat$driver != "HabitatLoss"),] # 1712


mod.ph <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * ph_fixed, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=phdat)

anova(mod.ph, btt = ":") #  not significant


mod.ph2 <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + ph_fixed, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=phdat)


anova(mod.ph2, btt = "ph_fixed") #  not significant
anova(mod.ph2, btt = "driver") #   significant


saveRDS(mod.ph2, file = "Models/pHMod.rds")

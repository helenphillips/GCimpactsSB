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



mod.gsba <-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

summary(mod.gsba)
anova(mod.gsba, btt = "GSBA") # NOT SIGNIFICANT



# just large groups

taxa_dat <- hedges[hedges$GSBA %in% c("Acari", "Collembola",  "Earthworms", "Nematodes"),]
table(taxa_dat$driver, taxa_dat$GSBA)



mod.gsba.taxa <-rma.mv(
  yi=effect,
  V=var, 
  mods=~- 1 + driver * GSBA, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=taxa_dat)
anova(mod.gsba.taxa, btt = ":") #  SIGNIFICANT

beep()
summary(mod.gsba.taxa)


saveRDS(mod.gsba.taxa, file = "Models/GSBAMod.rds")


t_dat <-predict(mod.gsba.taxa, newmods=rbind(c(0,0,0,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,0,0, 1,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,0,0, 0,1,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             # intercept (climate change)
                                             
                                             
                                             c(1,0,0,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(1,0,0,0,0, 1,0,0,  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(1,0,0,0,0, 0,1,0,  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
                                             c(1,0,0,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
                                             #  HabitatLoss                                        
                                             
                                             c(0,1,0,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,1,0,0,0, 1,0,0,  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,1,0,0,0, 0,1,0,  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                                             c(0,1,0,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
                                             # Invasives 
                                             
                                             c(0,0,1,0,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,1,0,0, 1,0,0,  0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,1,0,0, 0,1,0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                                             c(0,0,1,0,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
                                             # LUI
                                             
                                             c(0,0,0,1,0, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,1,0, 1,0,0,  0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,1,0, 0,1,0,  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
                                             c(0,0,0,1,0, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
                                             # NutrientEnrichment
                                             
                                             c(0,0,0,0,1, 0,0,0,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,0,1, 1,0,0,  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                                             c(0,0,0,0,1, 0,1,0,  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                                             c(0,0,0,0,1, 0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)),
                # Pollution 
                addx=TRUE, digits=2) #

gcds <- rep(c("Climate", "habitat", "invasives", "lui", "nutrient", "pollution"), each = 4)
taxa <- rep(c("acari","collembola", "earthworms", "nematodes"), times = 6)

slabs <- paste(gcds, taxa)

# exclude habitat frag and invasives (little data)
subs <- c(1:4, 13:24)

par(mar=c(3, 8, 1, 1))
forest(t_dat$pred[subs], sei=t_dat$se[subs], slab=slabs[subs],  xlab="Effect Size", xlim=c(-.4,.7))


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


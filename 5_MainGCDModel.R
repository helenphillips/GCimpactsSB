
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






## FULL MODEL ---------




# If fittiung with REML (which is recommended) then can't compare models 
# (https://stats.stackexchange.com/questions/48671/what-is-restricted-maximum-likelihood-and-when-should-it-be-used)
# https://stats.stackexchange.com/questions/517857/testing-the-effect-of-moderators-in-metafor-package
# so test moderators using the anova(singleMod) way



hedges <- read.csv("Data/03_Data/HedgesData_cleaned.csv")







############################
## FULL ANALYSIS WITH CROSSED RANDOM EFFECTS
########################





mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

anova(mod.1, btt = ":") # should be testing the levels that are the interactions
## not significant

# Ultimately, this interaction is not needed


mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)
beep()


anova(mod.1b, btt = "driver") # significant
anova(mod.1b, btt = "Body") # not significant



qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
# not bad


mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

beep()

anova(mod.2, btt = "driver") # significant




## we stick with mod.2
saveRDS(mod.2, file = "Models/MainMod.rds")




funnel(mod.2, main="Standard Error")
funnel(mod.2, yaxis="vi", main="Sampling Variance")
funnel(mod.2, yaxis="seinv", main="Inverse Standard Error")
funnel(mod.2, yaxis="vinv", main="Inverse Sampling Variance")


## Using publication bias test from Nakagawa


# Application of Equation 21 from the main text
publication.bias.model.r.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + driver + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)


# print(publication.bias.model.r.se,digits=3)

# extracting the mean and 95% confidence intervals
estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.r.se)


## Indicates bias
## Can account for that bias in the estimates



publication.bias.model.time.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~ -1 + year.c + driver + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)




# preparation to get marginalized mean (when vi = 0)
res.publication.bias.model.r.v.1 <- qdrg(object = publication.bias.model.time.allin, data = hedges, at = list(var = 0, year.c = 0))
# marginalized overall mean at vi = 0 and year.c = 0; also weights = "prop" or "cells" average things over proportionally. if not specified, all groups (levels) get the same weights
#overall.res.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = ~1, df = 104 - 7, weights = "prop") # using effect size - 7 

# marginalised means for different levels for driver
mm.publication.bias.model.r.v.1 <- emmeans(res.publication.bias.model.r.v.1, specs = "driver")

# comparing with results without correcting for publication bias
# this model (without the intercept) are whats reported in the manuscript
publication.bias.model.r.v.1b <-  rma.mv(
  yi=effect,
  V=var, 
  mods=~-1 + driver, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=hedges)

summary(publication.bias.model.r.v.1b)


estimates.publication.bias.model.r.v.1 <- estimates.CI2(mm.publication.bias.model.r.v.1)
estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)

table.comparing.captivity.levels <- merge(estimates.publication.bias.model.r.v.1,
                                          estimates.publication.bias.model.r.v.1b,
                                          by="estimate",
                                          all.x=T)

# rounding estimates
table.comparing.captivity.levels <- table.comparing.captivity.levels %>% mutate(across(where(is.numeric), round, 2))


table.comparing.captivity.levels <- data.frame(driver = table.comparing.captivity.levels[,1],
                                               adjusted.mean=table.comparing.captivity.levels[,2],               adjusted.CI=paste0("[",table.comparing.captivity.levels[,3],",",table.comparing.captivity.levels[,4],"]"),
                                               unadjusted.mean=table.comparing.captivity.levels[,5],                                               unadjusted.CI=paste0("[",table.comparing.captivity.levels[,6],",",table.comparing.captivity.levels[,7],"]"))













# https://stats.stackexchange.com/questions/109841/presenting-results-of-a-meta-analysis-with-multiple-moderators
### compute the estimated/predicted correlation for each combination
test <-predict(mod.2, newmods=rbind(c(0,0,0,0,0),
                                    # climate with all different measurement types
                                    c(1,0,0,0,0),# habitat 
                                    c(0,1,0,0,0), # invasive
                                    c(0,0,1,0,0),# lui
                                    c(0,0,0,1,0),# nutrient
                                    c(0,0,0,0,1)# pollution
), addx=TRUE, digits=2) #

slabs <- c("Climate", "Fragmentation", "Invasives", "LUI", "Nutrient", "Pollution")
par(mar=c(3, 8, 1, 1))
forest(test$pred, sei=test$se, slab=slabs,  xlab="Effect Size", xlim=c(-.4,.7))


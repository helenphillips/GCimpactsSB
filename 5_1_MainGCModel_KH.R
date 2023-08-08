
## Script to run main GC model
# But now to have test = t, which is analogous to the Knapp-HArtung method

library(Hmisc)
library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished
library(emmeans)
library(dplyr)

setwd("C:/Users/helenp/WORK/GCimpactsSB")



## FULL MODEL ---------


hedges <- read.csv("Data/03_Data/HedgesData_cleaned_June2023.csv")


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
  test="t",
  digits=4,
  data=hedges)

beep()

anova(mod.1, btt = ":") # should be testing the levels that are the interactions
## with test = t, still not significant

# Ultimately, this interaction is not needed


mod.1b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  test="t",
  digits=4,
  data=hedges)
beep()


anova(mod.1b, btt = "driver") # significant
anova(mod.1b, btt = "Body") # not significant # even eith the different test



mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  test="t",
  data=hedges)

beep()

anova(mod.2, btt = "driver") # with the t-test (F test statistic) still significant
summary(mod.2)

## All patterns are the same (i.e., pollution still strongest, then lui then climate)
# fragmentation and invasives not significant
# nutrient is now positive and significant
# but all coefficients are larger in magnitude (in their original direction)

## we would still stick with mod.2
saveRDS(mod.2, file = "Models/MainMod_rerun_June2023_KHtest.rds")



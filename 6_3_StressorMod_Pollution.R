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






### POLLUTION
poll <- hedges[which(hedges$driver == "Pollution"),] # 850


fulltext <- read.csv("Data/April2022/Full Text Screening - Sheet1.csv")
fulltext <- fulltext[,c('PaperID', 'PollutionSource')]

poll <- merge(poll, fulltext, by.x = "ID", by.y = "PaperID", all.x = TRUE)

# pollution type
table(poll$GCDType)

poll <- poll[which(poll$GCDType %in% c("Metals", "Pesticides")),] # the only two that have
# enough data/studies across the different pollutant types and body sizes



poll$Body.Size <- as.factor(poll$Body.Size)
poll$Body.Size <- relevel(poll$Body.Size, ref = "Meso-fauna")




# pollution type
table(poll$GCDType, poll$PollutionSource)




## previously missing
poll$PollutionSource[which(poll$ID == 105)] <- "Agricultural"
poll$PollutionSource[which(poll$ID == 2011)] <- "Agricultural"
poll$PollutionSource[which(poll$ID == 3372)] <- "Others"

## changing to remove multiple sources
poll$PollutionSource[which(poll$ID == 395)] <- "Mining/Smelting"
poll$PollutionSource[which(poll$ID == 2433 )] <- "Waste/sewage"
poll$PollutionSource[which(poll$ID == 1648)] <- "Mining/Smelting"
poll$PollutionSource[which(poll$ID == 1850)] <- "Mining/Smelting"
poll$PollutionSource[which(poll$ID == 564)] <- "Industrial"
poll$PollutionSource[which(poll$ID == 2509)] <- "Industrial"
poll$PollutionSource[which(poll$ID == 2715)] <- "Urban/transport"


poll$PollutionSource[which(poll$PollutionSource == "Agricultural")] <- "Farming"
poll$PollutionSource[which(poll$PollutionSource == "Agricultural/livestock")] <- "Farming"
poll$PollutionSource[which(poll$PollutionSource == "Military/wars")] <- "Others"
poll$PollutionSource[which(poll$PollutionSource == "Natural/geogenic")] <- "Others"



poll$GCDTypeSource <- paste(poll$GCDType, "-", poll$PollutionSource)
aggregate(poll$ID, by = list(GCD = poll$GCDTypeSource), function(x) length(unique(x)))
poll$GCDTypeSource[which(poll$GCDTypeSource == "Pesticides - Industrial")] <- "Pesticides - Others"

table(poll$GCDTypeSource, poll$Body.Size)


table(poll$GCDTypeSource, poll$GSBA) # to check for taxonomic model
# all four are not bad



# measurement
table(poll$GCDType, poll$Measurement)

# body size
table(poll$Body.Size)
table(poll$GCDType, poll$Body.Size)



poll <- droplevels(poll)

saveRDS(poll, file = "Models/pollutionDataFrame.rds")




poll.mod.10<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDTypeSource * Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)


anova(poll.mod.10, btt = ":") #  not significant

poll.mod.11<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDTypeSource + Body.Size, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)


anova(poll.mod.11, btt = "GCD") #  nearly significant
anova(poll.mod.11, btt = "Size") #  not significant


poll.mod.12<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDTypeSource, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)
anova(poll.mod.12, btt = "GCD") #  nearly significant

poll.mod.13<-rma.mv(
  yi=effect,
  V=var, 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

summary(poll.mod.13)

saveRDS(poll.mod.13, file = "Models/pollutionMod.rds")




## pollution pub bias
poll$year.c <- as.vector(scale(poll$year, scale = F))



publication.bias.model.poll.se <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + sei, ## 
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.poll.se)


## Indicates bias in sei
## Can account for that bias in the estimates


publication.bias.model.poll.allin <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + var, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)



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
  data=poll)


estimates.publication.bias.model.r.v.1 <- estimates.CI(publication.bias.model.poll.allin)

estimates.publication.bias.model.r.v.1b <- estimates.CI(publication.bias.model.r.v.1b)
# Most previous code doesn't work, because its an intercept only model




## Would the pollutant type make a difference if we accounted for bias?
publication.bias.model.poll.allin_fakemodel <- rma.mv(
  yi=effect,
  V=var, 
  mods=~1 + year.c + var + GCDTypeSource, ## note var not SEI # no intercept
  random= list(~1|ID/UniqueID, 
               ~ 1 | Measurement),
  struct="CS",
  method="REML",
  digits=4,
  data=poll)

anova(publication.bias.model.poll.allin_fakemodel, btt = "GCD") #  nearly significant




## Looking at the raw values

plot(poll$effect ~ as.factor(poll$GCDTypeSource))
summary(poll$effect ~ as.factor(poll$GCDTypeSource))

tapply(poll$effect, as.factor(poll$GCDTypeSource), summary)


poll[which(poll$effect < -20),]
poll[which(poll$effect > 20),]

poll_cut <- poll[which(poll$effect > -20),]
poll_cut <- poll_cut[which(poll_cut$effect < 20),]
plot(poll_cut$effect ~ as.factor(poll_cut$GCDTypeSource))
stripchart(poll_cut$effect ~ as.factor(poll_cut$GCDTypeSource),
           vertical=TRUE,
           method = "jitter",
           pch = 19,
           col = 1:7,
           add = TRUE)


## Script to calculate effect sizes


library(metafor)
library(ggplot2)
library(beepr) # to make a sound when a model has finished


# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")

## LOAD THE DATA
dataDir <- "Data/February2022"

DataCleaned <- read.csv(file.path(dataDir, "processed", "0_2_alldata.csv"), header = TRUE) # load the cleaned data (Change the name)

### Split to abundance and richness
table(DataCleaned$Measurement)


# Sometimes zeros were missing that should actually be zeros
DataCleaned$Control_SD[which(DataCleaned$Control_mean == 0 & is.na(DataCleaned$Control_SD))] <- 0
DataCleaned$Treatment_SD[which(DataCleaned$Treatment_mean == 0 & is.na(DataCleaned$Treatment_SD))] <- 0

## There's studies that rounded their SD's so that they were zero
DataCleaned$Control_SD[which(DataCleaned$ID == 757 & DataCleaned$Control_SD == 0)] <- 0.01
DataCleaned$Treatment_SD[which(DataCleaned$ID == 757 & DataCleaned$Treatment_SD == 0)] <- 0.01

DataCleaned$Control_SD[which(DataCleaned$ID == 1857 & DataCleaned$Control_SD == 0)] <- 0.01
DataCleaned$Treatment_SD[which(DataCleaned$ID == 1857 & DataCleaned$Treatment_SD == 0)] <- 0.01


## One case has an NA for SD, when it is actually a zero (according to the paper)
DataCleaned$Treatment_SD[which(DataCleaned$ID == 129 & is.na(DataCleaned$Treatment_SD))] <- 0



## Which measurements are there that have nas in Sds
table(DataCleaned$Measurement[which(!(is.na(DataCleaned$Control_mean)) & is.na(DataCleaned$Control_SD))])
## That's fine (11 abundance and 1 biomass)



# there are some climate cases that have no SDs
DataCleaned <- DataCleaned[which(!(is.na(DataCleaned$Control_SD))),] 
DataCleaned <- DataCleaned[which(!(is.na(DataCleaned$Treatment_SD))),] # 3440



## How many cases (and which drivers) now have zeros
nrow(DataCleaned[which(DataCleaned$Control_mean == 0),]) # 127
nrow(DataCleaned[which(DataCleaned$Control_SD   == 0),]) # 159

table(DataCleaned$driver[which(DataCleaned$Control_SD   == 0)]) # Mainly lui, nutrient enrichment and pollution


nrow(DataCleaned[which(DataCleaned$Treatment_mean  == 0),]) # 161
nrow(DataCleaned[which(DataCleaned$Treatment_SD  == 0),]) # 197
table(DataCleaned$driver[which(DataCleaned$Treatment_SD   == 0)]) # Mainly lui,  nutrient enrichment and pollution




hedges <- escalc(measure = "SMD", # log response ratio ("ROM" in metafor)
                      m2i = Control_mean, # group 2 corresponds to the control group
                      sd2i = Control_SD,
                      n2i = Control_N,
                      m1i = Treatment_mean, # group 1 is the treatment group
                      sd1i = Treatment_SD,
                      n1i = Treatment_N,
                      data = DataCleaned)
## get a warning message here because of the zeros, which gets dealt with shortly
# 3440

## 99 NAs in hedges
# When theres no Ns, or when there are zeros, in the Sds?


hedges <- hedges[which(!(is.na(hedges$yi))),]


# Rename Effect sizes
hedges$effect <- hedges$yi
hedges$var <- hedges$vi

hist(hedges$effect)
summary(hedges$effect)

hist(hedges$var)
summary(hedges$var)

## What have the high effect sizes and var

head(hedges[which(hedges$effect  > 75),]) # one case, and is legit



head(hedges[which(hedges$var > 2),])
table(hedges$Measurement[which(hedges$var > 10)])
head(hedges[which(hedges$var > 10),])


# TODO: check the high vars with Lea



# Save data
write.csv(hedges, "Data/03_Data/HedgesData.csv")
# 3307



hedges$UniqueID <- paste(hedges$ID, hedges$Case_ID, hedges$driver)


## TRYING TO ESTABLISH WHETHER ALL EFFECT SIZES CAN BE USED ----------

measurement.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)

qqnorm(residuals(measurement.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(measurement.mod.1,type="pearson"),col="red")
# not bad


## Everything is significantly different.
## So split everything via measurement type



## Make biomass the baseline
# Maybe biomass can be with richness??
biomass_baseline <- hedges
biomass_baseline$Measurement <- as.factor(biomass_baseline$Measurement)
biomass_baseline$Measurement <- relevel(biomass_baseline$Measurement, ref = "Biomass")

biomass.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=biomass_baseline)
## As I suspected, Richness (and Shannon) are not significantly different from Biomass



## Let's look at bodysize

bodysize.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Body.Size, ## Want to know if sign diff from Abundance
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)

summary(bodysize.mod.1)


## FULL MODEL ---------

hedges <- hedges[!(hedges$Body.Size == ""),]
hedges <- hedges[!(hedges$Measurement == "FunctionalRichness"),]
hedges <- hedges[!(hedges$Measurement == "GroupAbundance"),]


mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Measurement + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)

beep()


qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
# not bad


mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver + Measurement -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)


anova(mod.1, mod.2) # not significantly different
beep()



mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~driver  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)
anova(mod.2, mod.3) # that is significantly different


mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=hedges)
anova(mod.2, mod.4) # that is significantly different





## we stick with mod.2


### compute the estimated/predicted correlation for each combination
test <-predict(mod.2, newmods=rbind(c(1,0,0,0,0,0,0,0,0,0,0),
                             c(1,0,0,0,0,0,1,0,0,0,0),
                             c(1,0,0,0,0,0,0,1,0,0,0),
                             c(1,0,0,0,0,0,0,0,1,0,0),
                             c(1,0,0,0,0,0,0,0,0,1,0),
                             c(1,0,0,0,0,0,0,0,0,0,1), # climate with all different measurement types
                             c(0,1,0,0,0,0,0,0,0,0,0),
                             c(0,1,0,0,0,0,1,0,0,0,0),
                             c(0,1,0,0,0,0,0,1,0,0,0),
                             c(0,1,0,0,0,0,0,0,1,0,0),
                             c(0,1,0,0,0,0,0,0,0,1,0),
                             c(0,1,0,0,0,0,0,0,0,0,1), # habitat 
                             c(0,0,1,0,0,0,0,0,0,0,0),
                             c(0,0,1,0,0,0,1,0,0,0,0),
                             c(0,0,1,0,0,0,0,1,0,0,0),
                             c(0,0,1,0,0,0,0,0,1,0,0),
                             c(0,0,1,0,0,0,0,0,0,1,0),
                             c(0,0,1,0,0,0,0,0,0,0,1),  # invasive
                             c(0,0,0,1,0,0,0,0,0,0,0),
                             c(0,0,0,1,0,0,1,0,0,0,0),
                             c(0,0,0,1,0,0,0,1,0,0,0),
                             c(0,0,0,1,0,0,0,0,1,0,0),
                             c(0,0,0,1,0,0,0,0,0,1,0),
                             c(0,0,0,1,0,0,0,0,0,0,1), # lui
                             c(0,0,0,0,1,0,0,0,0,0,0),
                             c(0,0,0,0,1,0,1,0,0,0,0),
                             c(0,0,0,0,1,0,0,1,0,0,0),
                             c(0,0,0,0,1,0,0,0,1,0,0),
                             c(0,0,0,0,1,0,0,0,0,1,0),
                             c(0,0,0,0,1,0,0,0,0,0,1), # nutrient
                             c(0,0,0,0,0,1,0,0,0,0,0),
                             c(0,0,0,0,0,1,1,0,0,0,0),
                             c(0,0,0,0,0,1,0,1,0,0,0),
                             c(0,0,0,0,0,1,0,0,1,0,0),
                             c(0,0,0,0,0,1,0,0,0,1,0),
                             c(0,0,0,0,0,1,0,0,0,0,1) # pollution
                             ), addx=TRUE, digits=2) #

slab=slabs,  
forest(test$pred, sei=test$se, xlab="Effect Size", xlim=c(-.4,.7))

## SPLITTING DRIVERS --------

climate <- hedges[which(hedges$driver == "Climate"),] # 462

table(climate$GCDType)
table(climate$GCDType, climate$Measurement)


## Remove the very underrepresented measurement types
climate <- climate[which(climate$Measurement != "Evenness"),]
climate <- climate[which(climate$Measurement != "Simpson"),]


## going to remove vegetation as a gcd type. 2 lines are thickness of moss, and the 10 other lines are canopy gaps due to ice storm
climate <- climate[which(climate$GCDType != "Vegetation"),]

# For now will also remove the precip+temp
climate <- climate[which(climate$GCDType != "Precipitation+Temp"),]

# Checking body size
table(climate$Body.Size)
climate$Body.Size[which(climate$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
climate$Body.Size[which(climate$Body.Size %in% c("Macro-arthropods"))] <- "Macro-fauna"


climate.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)

summary(climate.mod.1)

qqnorm(residuals(climate.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(climate.mod.1,type="pearson"),col="red")


climate.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)


anova(climate.mod.1, climate.mod.2) # not diff


climate.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)

anova(climate.mod.1, climate.mod.3) # not diff, but lower p-val


climate.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=climate)
anova(climate.mod.3, climate.mod.4) # not diff, but lower p-val




qqnorm(residuals(climate.mod.4,type="pearson"),main="QQ plot: residuals")
qqline(residuals(climate.mod.4,type="pearson"),col="red") # fine


y<-summary(climate.mod.4)$b
ci_l<-summary(climate.mod.4)$ci.lb
ci_h<-summary(climate.mod.4)$ci.ub

fg1<-data.frame(cbind(y,ci_l,ci_h))
colnames(fg1)[1]<-"y"
colnames(fg1)[2]<-"ci_l"
colnames(fg1)[3]<-"ci_h"
fg1$GCD<-c("Fire","Gas", "Temperature",
           "UVB Radiation", "Water Availability-Drought", "Water Availability-Flood")
fg1$GCD<-as.factor(fg1$GCD)


p <- ggplot(fg1, aes(x=GCD, y=y, ymin=ci_l, ymax=ci_h))+
  geom_point(aes(size = 1)) +
  geom_pointrange()+
  geom_hline(yintercept = 0, linetype=2)+
  coord_flip()+
  theme_bw() +
  theme(axis.text=element_text(size=16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab('') +
  ylab ('Effect Size')
p


## LUI -------



lui <- hedges[which(hedges$driver == "LUI"),] # 911

table(lui$GCDType)
lui$GCDType[which(lui$GCDType %in% c("Defoliation"))] <- "Grazing"
lui$GCDType[which(lui$GCDType %in% c("Intensity", "Landscape", "Management", "Weeding", "Planting", "Intensification"))] <- "Management"
lui$GCDType[which(lui$GCDType %in% c("Logging"))] <- "Harvesting"
lui$GCDType[which(lui$GCDType %in% c("Water"))] <- "Irrigation"
lui$GCDType[which(lui$GCDType %in% c("degradation", "Disturbance"))] <- "Degradation"

lui <- lui[which(lui$GCDType != "Human population"),]



table(lui$GCDType, lui$Measurement)


## Remove the very underrepresented measurement types
lui <- lui[which(lui$Measurement != "FunctionalRichness"),]
lui <- lui[which(lui$Measurement != "Simpson"),]
lui <- lui[which(lui$Measurement != "Evenness"),]

## Checking body size
table(lui$Body.Size)
lui$Body.Size[which(lui$Body.Size %in% c("Arthropods (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
lui$Body.Size[which(lui$Body.Size %in% c("Macro-arthropods", "Macro-invertebrates"))] <- "Macro-fauna"



lui.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)

summary(lui.mod.1)

qqnorm(residuals(lui.mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(lui.mod.1,type="pearson"),col="red")
Anova(lui.mod.1)

lui.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)


anova(lui.mod.1, lui.mod.2) #  different (at 0.05 level)


lui.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=lui)


anova(lui.mod.1, lui.mod.3) #  definitely significantly different. Need measurement



## NUTRIENT ENRICHMENT -------



nutri <- hedges[which(hedges$driver == "NutrientEnrichment"),] # 820

table(nutri$GCDType)


table(nutri$GCDType, nutri$Measurement)


## Remove the very underrepresented measurement types
nutri <- nutri[which(nutri$Measurement != "Simpson"),]
nutri <- nutri[which(nutri$Measurement != "Evenness"),]

## Check body size
table(nutri$Body.Size)
nutri$Body.Size[which(nutri$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
nutri$Body.Size[which(nutri$Body.Size %in% c("Macro-arthropods", "Macro-invertebrates"))] <- "Macro-fauna"


nutri.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

summary(nutri.mod.1)


nutri.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

anova(nutri.mod.1, nutri.mod.2) #  not different (at 0.05 level)


nutri.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)
anova(nutri.mod.1, nutri.mod.3) #  significantly different (at 0.05 level). So move to mod2


nutri.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

anova(nutri.mod.2, nutri.mod.4) #  not quite significantly different (at 0.05 level)
anova(nutri.mod.1, nutri.mod.4) # but that is sign. different. So that makes no sense

nutri.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=nutri)

anova(nutri.mod.2, nutri.mod.5) #  Significantly different. 
## Need mod 4 in theory. But it doesn't quite make sense



## Invasives

invas <- hedges[which(hedges$driver == "Invasives"),] # 188

## GCD type
table(invas$GCDType)
# Probably need to remove plants-mixture
invas <- invas[which(invas$GCDType != "Plants-mixture"),]


# Measurement type

table(invas$GCDType, invas$Measurement)
invas <- invas[which(invas$Measurement == "Abundance"),] # just use abundance

## Body size
table(invas$GCDType, invas$Body.Size)

invas$Body.Size[which(invas$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"


invas.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  + Body.Size -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)


invas.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)

anova(invas.mod.1, invas.mod.2) ## That is significant


invas.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Body.Size  -1 , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=invas)

anova(invas.mod.1, invas.mod.3) ## That is also significant


summary(invas.mod.1)




### POLLUTION
poll <- hedges[which(hedges$driver == "Pollution"),] # 850

# pollution type
table(poll$GCDType)

poll$GCDType[which(poll$GCDType %in% c("Metals, Radionuclides", "Metals; PAH", "Metals; PCBs; PAHs", "Pesticides,Metals"))] <- "Mixture"
poll <- poll[which(poll$GCDType != ""),] 
poll <- poll[which(poll$GCDType != "Antibiotics"),] 
poll <- poll[which(poll$GCDType != "Chloride, Sodium"),] 
poll <- poll[which(poll$GCDType != "endocrine disruptors"),] 
poll <- poll[which(poll$GCDType != "Salinization"),] 
poll <- poll[which(poll$GCDType != "Sulphate"),] 


# measurement
table(poll$GCDType, poll$Measurement)

# body size
table(poll$Body.Size)


poll$Body.Size[which(poll$Body.Size %in% c("Arthropods (all sizes)", "Insects (all sizes)", "Invertebrates (all sizes)"))] <- "All sizes"
poll$Body.Size[which(poll$Body.Size %in% c("Macro-arthropods", "Macro-invertebrates"))] <- "Macro-fauna"



poll.mod.1<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement + Body.Size , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)

summary(poll.mod.1)


poll.mod.2<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType  + Body.Size, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)

anova(poll.mod.1, poll.mod.2) # Significant at 0.05 level

poll.mod.3<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType + Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.1, poll.mod.3) # not Significantly different



poll.mod.4<-rma.mv(
  yi=effect,
  V=var, 
  mods=~GCDType , ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.3, poll.mod.4) # that is not Significantly different



poll.mod.5<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
anova(poll.mod.3, poll.mod.5) # that is  Significantly different
## The type of GCD doesn't matter, just the measurement

summary(poll.mod.5)

poll.mod.5b<-rma.mv(
  yi=effect,
  V=var, 
  mods=~Measurement -1, ## 
  random= ~1|ID/UniqueID,
  struct="CS",
  method="ML",
  digits=4,
  data=poll)
summary(poll.mod.5b)

## Everything is significantly negative.
## (except simpson)





## The code created a colours qqplot
## Coloured qq plot https://stackoverflow.com/questions/42678858/q-q-plot-with-ggplot2stat-qq-colours-single-group
library(broom) ## for augment()
dda <- cbind(augment(mod.1),f=EffectSizes$driver)
dda = cbind(dda, setNames(qqnorm(dda$.resid, plot.it=FALSE), c("Theoretical", "Sample")))

ggplot(dda) + 
  geom_point(aes(x=Theoretical, y=Sample, colour=f))



## This code for creating a basic plot

y<-summary(mod.1)$b
ci_l<-summary(mod.1)$ci.lb
ci_h<-summary(mod.1)$ci.ub

fg1<-data.frame(cbind(y,ci_l,ci_h))
colnames(fg1)[1]<-"y"
colnames(fg1)[2]<-"ci_l"
colnames(fg1)[3]<-"ci_h"
fg1$GCD<-c("Climate change","Fragmentation", "Invasive species",
           "LUI", "Nutrient enrichment", "Pollution", rep(NA, 7))
fg1$GCD<-as.factor(fg1$GCD)


p <- ggplot(fg1, aes(x=GCD, y=y, ymin=ci_l, ymax=ci_h))+
  geom_point(aes(size = 1)) +
  geom_pointrange()+
  geom_hline(yintercept = 0, linetype=2)+
  coord_flip()+
  theme_bw() +
  theme(axis.text=element_text(size=16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab('') +
  ylab ('Effect Size')
p


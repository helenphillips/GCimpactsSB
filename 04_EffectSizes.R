## Script to calculate effect sizes


library(metafor)

# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")

## LOAD THE DATA
dataDir <- "Data/January2022"

DataCleaned <- read.csv(file.path(dataDir, "processed", "0_2_alldata.csv"), header = TRUE) # load the cleaned data (Change the name)

### Split to abundance and richness
table(DataCleaned$Measurement)

abundance <- DataCleaned[which(DataCleaned$Measurement == "Abundance"),] # 2024 # 2395

# Sometimes zeros were missing that should actually be zeros
abundance$Control_SD[which(abundance$Control_mean == 0 & is.na(abundance$Control_SD))] <- 0
abundance$Treatment_SD[which(abundance$Treatment_mean == 0 & is.na(abundance$Treatment_SD))] <- 0

# there are some climate cases that have no SDs
abundance <- abundance[which(!(is.na(abundance$Control_SD))),] 
abundance <- abundance[which(!(is.na(abundance$Treatment_SD))),] # 2378


## There's studies that rounded their SD's so that they were zero
abundance$Control_SD[which(abundance$ID == 757 & abundance$Control_SD == 0)] <- 0.01
abundance$Treatment_SD[which(abundance$ID == 757 & abundance$Treatment_SD == 0)] <- 0.01

abundance$Control_SD[which(abundance$ID == 1857 & abundance$Control_SD == 0)] <- 0.01
abundance$Treatment_SD[which(abundance$ID == 1857 & abundance$Treatment_SD == 0)] <- 0.01


EffectSizes <- escalc(measure = "ROM", # log response ratio ("ROM" in metafor)
                 m2i = Control_mean, # group 2 corresponds to the control group
                 sd2i = Control_SD,
                 n2i = Control_N,
                 
                 m1i = Treatment_mean, # group 1 is the treatment group
                 sd1i = Treatment_SD,
                 n1i = Treatment_N,
                 
                 data = abundance)


## get a warning message here that I should deal with



# we now have the log response ratio for each case
summary(EffectSizes$yi) # 221 NAs
# with their variance
summary(EffectSizes$vi) # 224 Nas

# Rename Effect sizes
EffectSizes$LRR <- EffectSizes$yi
EffectSizes$VarLRR <- EffectSizes$vi


## What have the high varLRR

EffectSizes[which(EffectSizes$VarLRR > 10),]

# TODO: check these again


## MISSING EFFECT SIZES (BECAUSE OF LEGIT ZEROS)
missing <- EffectSizes[which(is.na(EffectSizes$yi) | is.na(EffectSizes$vi)),]
table(missing$driver)


EffectSizes$LRR[which(EffectSizes$driver == "Pollution" & is.na(EffectSizes$yi))] <- min(EffectSizes$yi[EffectSizes$driver == "Pollution"], na.rm = TRUE)
EffectSizes$LRR[which(EffectSizes$driver == "Climate" & is.na(EffectSizes$yi))] <-min(EffectSizes$yi[EffectSizes$driver == "Climate"], na.rm = TRUE)
EffectSizes$LRR[which(EffectSizes$driver == "Invasives" & is.na(EffectSizes$yi))] <-min(EffectSizes$yi[EffectSizes$driver == "Invasives"], na.rm = TRUE)
EffectSizes$LRR[which(EffectSizes$driver == "LUI" & is.na(EffectSizes$yi))] <-min(EffectSizes$yi[EffectSizes$driver == "LUI"], na.rm = TRUE)
EffectSizes$LRR[which(EffectSizes$driver == "NutrientEnrichment" & is.na(EffectSizes$yi))] <-min(EffectSizes$yi[EffectSizes$driver == "NutrientEnrichment"], na.rm = TRUE)



EffectSizes$VarLRR[which(EffectSizes$driver == "Pollution" & is.na(EffectSizes$vi))] <- median(EffectSizes$vi[EffectSizes$driver == "Pollution"], na.rm = TRUE)
EffectSizes$VarLRR[which(EffectSizes$driver == "Climate" & is.na(EffectSizes$vi))] <- median(EffectSizes$vi[EffectSizes$driver == "Climate"], na.rm = TRUE)
EffectSizes$VarLRR[which(EffectSizes$driver == "Invasives" & is.na(EffectSizes$vi))] <- median(EffectSizes$vi[EffectSizes$driver == "Invasives"], na.rm = TRUE)
EffectSizes$VarLRR[which(EffectSizes$driver == "LUI" & is.na(EffectSizes$vi))] <- median(EffectSizes$vi[EffectSizes$driver == "LUI"], na.rm = TRUE)
EffectSizes$VarLRR[which(EffectSizes$driver == "NutrientEnrichment" & is.na(EffectSizes$vi))] <-median(EffectSizes$vi[EffectSizes$driver == "NutrientEnrichment"], na.rm = TRUE)



# Save data
write.csv(EffectSizes, "Data/03_Data/EffectSizes.csv")


EffectSizes$UniqueID <- paste(EffectSizes$ID, EffectSizes$Case_ID, EffectSizes$driver)

mod.1<-rma.mv(
  yi=LRR,
  V=VarLRR, 
  mods=~driver,
  random= ~1|ID,
  struct="CS",
  method="ML",
  digits=4,
  data=EffectSizes)

qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")




y<-summary(mod.1)$b
ci_l<-summary(mod.1)$ci.lb
ci_h<-summary(mod.1)$ci.ub

fg1<-data.frame(cbind(y,ci_l,ci_h))
colnames(fg1)[1]<-"y"
colnames(fg1)[2]<-"ci_l"
colnames(fg1)[3]<-"ci_h"
fg1$GCD<-c("Climate change","Fragmentation", "Invasive species",
           "LUI", "Nutrient enrichment", "Pollution")
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



test <- EffectSizes[-(which(EffectSizes$VarLRR > 10)),]

mod.test<-rma.mv(
  yi=LRR,
  V=VarLRR, 
  mods=~driver,
  random= ~1|ID,
  struct="CS",
  method="ML",
  digits=4,
  data=test)

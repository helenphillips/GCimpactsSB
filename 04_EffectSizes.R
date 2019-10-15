## Script to calculate effect sizes


library(metafor)

# DataCleaned <- read.csv("Data/02_CleanData/data.csv", h = TRUE) # load the cleaned data (Change the name)


EffectSizes <- escalc(measure = "ROM", # log response ratio ("ROM" in metafor)
                 m2i = Control_mean, # group 2 corresponds to the control group
                 sd2i = Control_SD,
                 n2i = Control_N,
                 
                 m1i = Treatment_mean, # group 1 is the treatment group
                 sd1i = Treatment_SD,
                 n1i = Treatment_N,
                 
                 data = DataCleaned)


# we now have the log response ratio for each case
summary(EffectSizes$yi)
# with their variance
summary(EffectSizes$vi)

# Rename Effect sizes
EffectSizes$LRR <- EffectSizes$yi
EffectSizes$VarLRR <- EffectSizes$vi

# Save data
write.csv(EffectSizes, "Data/03_Data/EffectSizes.csv")



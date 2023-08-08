setwd("C:/Users/helenp/WORK/GCimpactsSB")
FigoutPath <- "C:/Users/helenp/WORK/GCimpactsSB/Figures"


library("viridis") 
library("plotrix") # cornerlabels
library(metafor)


# Full dataset
hedges <- read.csv("Data/03_Data/HedgesData_cleaned_June2023.csv")


# Models (and subsets of data)
climate.mod.5 <- readRDS(file = "Models/ClimateMod_june2023.rds")
df_climate <- as.data.frame(climate.mod.5$yi) # effect sizes

lui.mod.2 <- readRDS(file = "Models/LUIMod_redo_june2023.rds")
df_lui <- as.data.frame(lui.mod.2$yi) # effect sizes

poll.mod.22 <-  readRDS(file = "Models/pollutionMod_june2023.rds")
df_poll <- as.data.frame(poll.mod.22$yi) # effect sizes

nutri.mod.2 <- readRDS(file = "Models/nutriMod_june2023.rds")
df_nutri <- as.data.frame(nutri.mod.2$yi) # effect sizes

invasivemod <- readRDS(file ="Models/invasiveMod_june2023.rds")
df_invas <- as.data.frame(invasivemod$yi) # effect sizes

fragdat <- hedges[which(hedges$driver == "HabitatLoss"),]


# Colour for bars
rbcol <- adjustcolor(viridis(3)[1], alpha.f = 0.3) # the purple colour
globalcol <- adjustcolor(viridis(3)[2], alpha.f = 0.3) # the purple colour


pdf(file = file.path(FigoutPath, "EffectSizeDistributions.pdf"),
    width = 6, height = 10)

layout(matrix(c(1,1,2,3, 4, 5, 6, 7), nrow=4, byrow = TRUE))
hist(hedges$yi, breaks = 200,
     xlim = c(-20, 20), 
     xlab = "Effect Sizes", 
     main = "", 
     col = globalcol)
corner.label(label="(a)")



#par(mfrow = c(3, 2))
hist(df_climate$`climate.mod.5$yi`, breaks = 200,
     xlim = c(-5, 5), 
     xlab = "Effect Sizes", 
     main = "", 
     col = rbcol)
corner.label(label="(b)")


hist(df_lui$`lui.mod.2$yi`, breaks = 100,
    xlim = c(-10, 5), 
     xlab = "Effect Sizes", 
     main = "",
    col = rbcol)
corner.label(label="(c)")


hist(df_poll$`poll.mod.22$yi`, breaks = 200,
     xlim = c(-20, 10), 
     xlab = "Effect Sizes", 
     main = "",
     col = rbcol)
corner.label(label="(d)")

hist(df_nutri$`nutri.mod.2$yi`, breaks = 100,
    xlim = c(-5, 10), 
     xlab = "Effect Sizes", 
     main = "",
    col = rbcol)
corner.label(label="(e)")

hist(df_invas$`invasivemod$yi`, breaks = 100,
     xlim = c(-5, 5), 
     xlab = "Effect Sizes", 
     main = "",
     col = rbcol)
corner.label(label="(f)")

hist(fragdat$yi, breaks = 100,
     xlim = c(-5, 5), 
     xlab = "Effect Sizes", 
     main = "",
     col = rbcol)
corner.label(label="(g)")
dev.off()



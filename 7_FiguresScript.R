
setwd("C:/Users/helenp/WORK/GCimpactsSB")


library(Hmisc)
library(metafor)
library(ggplot2)
library(dplyr)
library(maptools)
library(ggmap)
library(rnaturalearth)
library(countrycode)
library(viridis)
library(sf)


# install.packages("remotes")
# remotes::install_github("liamgilbey/ggwaffle")
library(ggwaffle)

# devtools::install_github("MathiasHarrer/dmetar")
library(dmetar) # For 3-level I^2

allDat <- read.csv("Data/03_Data/HedgesData_cleaned_June2023.csv")
meta <- read.csv("Data/June2023/processed/metadata.csv")


FigoutPath <- "C:/Users/helenp/WORK/GCimpactsSB/Figures"

## coordinates of studies
coords <- read.csv("Data/CoordinatesFixed.csv")
coords <- coords[!duplicated(coords[,c(1, 5:6)]),]
coords<- coords[which(!(is.na(coords$lat_dd))),]
coords<- coords[which(!(is.na(coords$long_dd))),]


coords <- st_as_sf(coords, coords = c("long_dd", "lat_dd"), 
                  crs = 4326, agr = "constant")





corner.label2 <- function(label = NULL, x = -1, y = 1, xoff = NA, yoff = NA, 
                          figcorner = FALSE, ...) {
  if (is.na(xoff)) 
    xoff <- strwidth("m")/2
  if (is.na(yoff)) 
    yoff <- strheight("m")/2
  par.usr <- par("usr")
  xpos <- par.usr[(3 + x)/2]
  ypos <- par.usr[(3 + y)/2 + 2]
  if (figcorner) {
    par.pin <- par("pin")
    xplotrange <- par.usr[2] - par.usr[1]
    yplotrange <- par.usr[4] - par.usr[3]
    par.mai <- par("mai")
    xmar <- xplotrange * par.mai[3 + x]/par.pin[1]
    ymar <- yplotrange * par.mai[2 + y]/par.pin[2]
    xpos <- xpos + x * xmar
    ypos <- ypos + y * ymar
  }
  if (!is.null(label)) {
    if (figcorner) 
      par(xpd = TRUE)
    text(xpos - x * xoff, ypos - y * yoff, label, adj = c((1 + 
                                                             x)/2, (1 + y)/2), ...)
    if (figcorner) 
      par(xpd = FALSE)
  }
  return(list(x = xpos, y = ypos))
}

## MAP




meta$Country[which(meta$Country == "Czech Rpublic" | meta$Country == "Czechia"  )] <- "Czech Republic"


meta$Country[which(meta$Country == "England" 
                   | meta$Country == "Northern Ireland"
                   | meta$Country == "Scotland"
                   | meta$Country == "the UK"
                   | meta$Country == "United Kingdom"
                   | meta$Country == "Wales")] <- "UK"

meta$Country[which(meta$Country == "United States" 
                   | meta$Country == "United States - Hawaii"
                   | meta$Country == "the US")] <- "USA"

meta$Country[which(meta$Country == "The Netherlands")] <- "Netherlands"
meta$Country[which(meta$Country == "Neherlands")] <- "Netherlands"

meta$Country[grep("Ivoire", meta$Country)] <- "Ivory Coast"

meta$Country[which(meta$Country == "CÃ´teÂ d'Ivoire")] <- "Ivory Coast"
meta$Country[which(meta$Country == "Falkland Islands/Antarctica")] <- "Falkland Islands"


meta <- meta[-which(meta$Country %in% c("Africa",
                                        "Brazil/Uruguay",
                                        "France and Italy"    ,                        
                                        "France, The Netherlands, Switzerland, Canada",
                                        "Ireland, Poland",
                                        "Netherlands, England, Portugal, Germany",
                                        "Portugal,Brazil"   ,
                                        "The Netherlands, UK" ,
                                        "UK, Netherlands, Denmark",
                                        "Sweden, UK, Czech Republic and Greece")),] # Some of these options don't exist though now

meta <- meta[which(meta$Country != ""),]



stud <- meta %>% 
  group_by(Country) %>% 
  summarise(no.stu = n())


# Number of studies in USA and China
stud[which(stud$Country == 'China'),] # 96
stud[which(stud$Country == 'USA'),] # 100

# Number of cases
china <- meta$ID[which(meta$Country == "China")]
chinadat <- allDat[which(allDat$ID %in% china),]
nrow(chinadat) # 355
355/3161 * 100 
usa <- meta$ID[which(meta$Country == "USA")]
usadat <- allDat[which(allDat$ID %in% usa),]
nrow(usadat) # 492
492/3161 * 100
# harmonize country names
stud <- stud %>%
  mutate(iso =  countrycode(stud$Country, origin = "country.name", destination = "iso3c")) %>%
  dplyr::group_by(iso) %>%
  summarise(no.stud = sum(no.stu, na.rm = TRUE))


wm <- ne_countries(scale = 110, returnclass = "sf") %>%
  left_join(stud, by = c("adm0_a3" = "iso"))

wm <- sf::st_transform(wm,"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")



## categorical map
wm$no.stud.cat <- NA
wm$no.stud.cat[which(wm$no.stud < 5)] <- "<5"
wm$no.stud.cat[which(wm$no.stud >= 5 & wm$no.stud < 10) ] <- "5 - 10"
wm$no.stud.cat[which(wm$no.stud >= 10 & wm$no.stud < 20) ] <- "10 - 20"
wm$no.stud.cat[which(wm$no.stud >= 20 & wm$no.stud < 30) ] <- "20 - 30"
wm$no.stud.cat[which(wm$no.stud >= 30 & wm$no.stud < 40) ] <- "30 - 40"
wm$no.stud.cat[which(wm$no.stud >= 90 & wm$no.stud < 100) ] <- "90 - 100"
wm$no.stud.cat[which(wm$no.stud >= 100) ] <- "100+"
wm$no.stud.cat <- as.factor(wm$no.stud.cat)

wm$no.stud.cat <- factor(wm$no.stud.cat, levels = c("<5", "5 - 10", "10 - 20", "20 - 30", "30 - 40", "90 - 100", "100+"))


# jpeg(filename = file.path(FigoutPath, "GlobalMap.jpg"),
#      width = 10, height = 6, units = "in", res = 600)

pdf(file = file.path(FigoutPath, "GlobalMap.pdf"),
     width = 10, height = 6)




map <- ggplot(wm)+
  geom_sf(aes(fill = (no.stud.cat)))+
  scale_fill_viridis_d(option = "plasma", 
                       begin = 0,
                       end = 0.9,
                       na.value = "gray",
                       name = "Number of\npublications")+
  # geom_point(data = coords, aes(x = coords$long_dd, y = coords$lat_dd),col="red", size=5)+
  geom_sf(data = coords, size = 0.8) + 
  theme_bw()+
  theme(legend.position = c(0,0),
        # legend.position = "left",
        legend.justification = c(-0.2, -0.1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=14))

map

dev.off()


#### -------------
## PLOT OF THE HARMONISED TAXANOMIC GROUPS

xx <- table(allDat$GSBA)[rev(order(table(allDat$GSBA)))]

xx_counts <- as.numeric(xx)
xx_counts_top <- xx_counts[1:13]
xx_counts_all <- c(xx_counts_top, sum(xx_counts[14:length(xx_counts)]))
counts_original<-xx_counts_all
xx_counts_all <- round(xx_counts_all/sum(xx_counts_all) * 100, 2)

xx_names <- names(xx)
xx_names_top <- xx_names[1:13]
xx_names_all <- c(xx_names_top, "Others")


counts <- xx_counts_all
counts_names<-sprintf("%s (%s)", xx_names_all, 
                      scales::percent(round(counts_original/sum(counts_original), 4)))
names(counts)<-counts_names


newcols <- c(
"#FF0000", "#FF7F00", "#FFD400", "#FFFF00",
"#BFFF00", "#6AFF00", "#00EAFF", "#0095FF",
"#AA00FF", "#FF00AA", "#EDB9B9", "#E7E9B9",
"#B9EDE0", "#B9D7ED" )#, "#DCB9ED", "#8F2323",
#"#8F6A23", "#4F8F23", "#23628F", "#6B238F",
# "#000000", "#737373", "#CCCCCC") 

#newcols2 <- c(
# "#2f4f4f", "#6b8e23", "#a0522d", "#191970",
#"#ff0000", "#ffd700", "#7fff00", "#00fa9a",
#"#00bfff", "#0000ff", "#ff00ff", "#dda0dd",
#"#ff1493", "#ffe4b5")

#  colramp <- colorRampPalette(c("blue", "red"))( length(counts) )
pdf(file = file.path(FigoutPath, "TaxanomicWaffle.pdf"),
    width = 10, height = 6)

waffle(counts, rows = 8, colors = newcols)
dev.off()


## -----------------------------------------------------
## Frequency plot of GCs and fauna size


taxadat <- allDat
taxadat$Body.Size <- as.factor(taxadat$Body.Size)
taxadat$driver <- as.factor(taxadat$driver)

taxadat$Body.Size <- factor(taxadat$Body.Size, levels = c("All sizes", "Micro-fauna", "Meso-fauna", "Macro-fauna"))
taxadat$driver <- factor(taxadat$driver, levels = c("Climate" ,"LUI" , "Pollution", "NutrientEnrichment", 
                                                    "Invasives",  "HabitatLoss"))

#jpeg(filename = file.path(FigoutPath, "BodySizexGCD.jpg"),
#     width = 10, height = 6, units = "in", res = 600)

pdf(file = file.path(FigoutPath, "BodySizexGCD.pdf"),
    width = 10, height = 6)


barplot(table(taxadat$Body.Size, taxadat$driver),
        names.arg = c("Climate change", "Land-use \nintensification", "Pollution", "Nutrient \nenrichment", "Invasive \nspecies", "Habitat \nfragmentation"),
        xlab = "",
        ylab = "Number of cases",
        axes = TRUE,
        col = rocket(4))
legend("topright", legend = c("All sizes", "Micro-fauna", "Meso-fauna", "Macro-fauna"),
       fill = rocket(4), bty = "n")


dev.off()



## ------------------------------------------------------


## MODEL FIGURES
# Main mod

# mod.2 <- readRDS(file = "Models/MainMod.rds")

mod.2 <- readRDS(file = "Models/MainMod_rerun_June2023.rds")
df_main <- as.data.frame(mod.2$X)



mainmodpred <-predict(mod.2, newmods=rbind(c(0,0,0,0,0), 
                                           c(1,0,0,0,0),
                                           c(0,1,0,0,0),
                                           c(0,0,1,0,0),
                                           c(0,0,0,1,0),
                                           c(0,0,0,0,1)
), addx=TRUE, digits=2) #

slabs <- c("Climate change", "Habitat\nfragmentation", "Invasive\nspecies", "Land-use\nintensification", "Nutrient\nenrichment", "Pollution")

order <- rev(c(1, 4, 6, 5, 3, 2)) # reverse of the actual order


# how many studies per GC
tapply(allDat$ID, allDat$driver, function(x) length(unique(x)))



pdf(file = file.path(FigoutPath, "DriversMod.pdf"),
     width = 10, height = 6)

par(mar=c(5, 10, 1, 7))
errbar(x =slabs[order], y = mainmodpred$pred[order], 
       yplus = mainmodpred$ci.ub[order],
       yminus=mainmodpred$ci.lb[order], cap=0.2,
       cex = 2)
abline(v=0, lty =2)
mtext("Effect size", side = 1, line = 3, cex = 1.5)


intce <- nrow(df_main) - (sum(df_main['driverLUI']) + 
                            sum(df_main['driverPollution']) +  
                            sum(df_main['driverNutrientEnrichment']) +  
                            sum(df_main['driverInvasives']) +  
                            sum(df_main['driverHabitatLoss'])) 



mtext(paste0("n = ", intce," (93)"), side = 4, line = 0, at = 6.1, las = 2)
mtext(paste0("n = ", sum(df_main['driverLUI'])," (238)"), side = 4, line = 0, at = 5.1, las = 2)
mtext(paste0("n = ", sum(df_main['driverPollution'])," (151)"), side = 4, line = 0, at = 4.1, las = 2)
mtext(paste0("n = ", sum(df_main['driverNutrientEnrichment'])," (157)"), side = 4, line = 0, at = 3.1, las = 2)
mtext(paste0("n = ", sum(df_main['driverInvasives'])," (36)"), side = 4, line = 0, at = 2.1, las = 2)
mtext(paste0("n = ", sum(df_main['driverHabitatLoss']), " (22)"), side = 4, line = 0, at = 1.1, las = 2)

# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(6.5,10,10,6.5), col="white", border=F)

dev.off()

## ---------------------------------------------------------------------
## Panel plot for three stressor graphs

## Climate change 

climate.mod.5 <- readRDS(file = "Models/ClimateMod_june2023.rds")
nrow(climate.mod.5$X) # this is the dataframe (ish)
df_climate <- as.data.frame(climate.mod.5$X)

climatedat <-predict(climate.mod.5, newmods=rbind(c(0,0,0),
                                                  c(1,0,0),
                                                  c(0,1,0),
                                                  c(0,0,1)
), addx=TRUE, digits=2) #


climate_slabs <- c("Carbon dioxide increase", "Temperature change", "Drought", "Flooding")


## lui

lui.mod.2 <- readRDS(file = "Models/LUIMod_redo_june2023.rds")
df_lui <- as.data.frame(lui.mod.2$X)


luidat <-predict(lui.mod.2, newmods=rbind(c(0,0,0,0,0,0),  # intercept
                                          c(0,0,0,0,1,0),
                                          c(0,0,0,0,0,1),
                                          
                                          c(1,0,0,0,0,0),# fire
                                          c(1,0,0,0,1,0),
                                          c(1,0,0,0,0,1),
                                          
                                          c(0,1,0,0,0,0), 
                                          c(0,1,0,0,1,0),
                                          c(0,1,0,0,0,1),# Harvesting   
                                          
                                          c(0,0,1,0,0,0),
                                          c(0,0,1,0,1,0),
                                          c(0,0,1,0,0,1),# organic
                                          
                                          c(0,0,0,1,0,0),
                                          c(0,0,0,1,1,0),
                                          c(0,0,0,1,0,1)# tillage
), addx=TRUE, digits=2) #

macro <- rev(c(1, 4, 7, 10, 13)) #3 just the macrofauna
micro <- macro + 2

lui_slabs <- rev(c("Grazing", "Fire", "Harvesting", "Conventional farming", "Tillage"))


# just micro
df_luimicro <- df_lui[which(df_lui$`Body.SizeMicro-fauna` == 1),]

# not micro
df_luinotmicro <- df_lui[which(df_lui$`Body.SizeMicro-fauna` == 0),]



## pollution

poll.mod.22 <-  readRDS(file = "Models/pollutionMod_june2023.rds")
df_poll <- as.data.frame(poll.mod.22$X)


polldat <-predict(poll.mod.22, newmods=rbind(c(0), c(1)), 
                  addx=TRUE, digits=2) #

poll_slabs <- c("Metals", "Pesticides")



## nutrient
nutri.mod.2 <- readRDS(file = "Models/nutriMod_june2023.rds")

df_nut <- as.data.frame(nutri.mod.2$X)


nutridat <-predict(nutri.mod.2, newmods=rbind(c(0,0,0,0,0,0,0),
                                              # intercept (synthetic Fertilizers)
                                              c(1,0,0,0,0,0,0),
                                              #  "Ca-liming + Wood ash"                                         
                                              c(0,1,0,0,0,0,0),
                                              # Compost
                                              c(0,0,1,0,0,0,0),
                                              # Manure + Slurry 
                                              c(0,0,0,1,0,0,0),
                                              # OMixture
                                              c(0,0,0,0,1,0,0),
                                              # rOther Organic fertilisers
                                              c(0,0,0,0,0,1,0),
                                              # Residue 
                                              c(0,0,0,0,0,0,1)# Sludge 
), addx=TRUE, digits=2) #

# Big gap is because I can't get the left margin to be bigger...annoying

nut_slabs <- c("Synthetic Fertilizers", "Ca-liming + Wood ash", "Compost", "Manure + Slurry", 
               "                        Multiple fertilizer types", 
               "Other Organic fertilisers", "Residue + Mulch", "Sludge")


nut_ord <- rev(c(1, 2, 3,8, 4,  7, 6, 5))




## create one big dataframe

panel_dat <- rbind(as.data.frame(nutridat)[nut_ord,1:6], 
                   rep(NA, 6), # sneaky way to get a gap in the plot
                   as.data.frame(polldat)[,1:6],
                   rep(NA, 6),
                   as.data.frame(luidat)[macro,1:6],
                   rep(NA, 6),
                   as.data.frame(climatedat)[,1:6] )


panel_dat$driver <- c(rep("nutri", length(nut_slabs)), NA, 
                      rep("poll", length(poll_slabs)), NA,
                      rep("lui", length(lui_slabs)), NA, 
                      rep("climate", times = 4) )


all_slabs <- c(nut_slabs[nut_ord], "", poll_slabs, "",  lui_slabs, "", climate_slabs)




#jpeg(filename = file.path(FigoutPath, "stressorsPanel.jpg"),
#     width = 12, height = 8, units = "in", res = 600)
 pdf(file = file.path(FigoutPath, "stressorsPanel.pdf"),
    width = 12, height = 8)


par(mar=c(5, 10, 2, 6))
errbar(x =all_slabs, y = panel_dat$pred, 
       yplus = panel_dat$ci.ub,
       yminus=panel_dat$ci.lb, cap=0.2,
       cex = 2)
abline(v=0, lty =2)

# tapply(allDat$ID, allDat$GCDType, function(x) length(unique(x)))


# Climate n
mtext("Effect size", side = 1, line = 3, cex = 1.5)
mtext(paste0("n = ", sum(df_climate['GCDTypeWaterAvailability-Flood']), " (15)"), side = 4, line = 0, at = 22.1, las = 2)
mtext(paste0("n = ", sum(df_climate['GCDTypeWaterAvailability-Drought']), " (28)"), side = 4, line = 0, at = 21.1, las = 2)
mtext(paste0("n = ", sum(df_climate$GCDTypeTemperature), " (42)"), side = 4, line = 0, at = 20.1, las = 2)
intce <- nrow(df_climate) - (sum(df_climate['GCDTypeWaterAvailability-Flood']) + sum(df_climate['GCDTypeWaterAvailability-Drought']) +  sum(df_climate$GCDTypeTemperature ))
mtext(paste0("n = ", intce, " (28)"), side = 4, line = 0, at = 19.1, las = 2)

# add in the co2 axis label

# LUI n's
mtext(paste0("n = ", sum(df_lui['GCDTypeFire']), " (27)"), side = 4, line = 0, at = 16.1, las = 2)
mtext(paste0("n = ", sum(df_lui['GCDTypeHarvesting']), " (36)"), side = 4, line = 0, at = 15.1, las = 2)
mtext(paste0("n = ", sum(df_lui['GCDTypeOrganic versus Inorganic']), " (38)"), side = 4, line = 0, at = 14.1, las = 2)
mtext(paste0("n = ", sum(df_lui['GCDTypeTillage']), " (73)"), side = 4, line = 0, at = 13.1, las = 2)

intce <- nrow(df_lui) - (sum(df_lui['GCDTypeFire']) + sum(df_lui['GCDTypeHarvesting']) + 
                           sum(df_lui['GCDTypeOrganic versus Inorganic']) +
                           sum(df_lui['GCDTypeTillage']))
mtext(paste0("n = ", intce, " (42)"), side = 4, line = 0, at = 17.1, las = 2)



# pollution n's
mtext(paste0("n = ", sum(df_poll['GCDTypePesticides']), " (51)"), side = 4, line = 0, at = 11.1, las = 2)
intce <- nrow(df_poll) - (sum(df_poll['GCDTypePesticides']))
mtext(paste0("n = ", intce, " (77)"), side = 4, line = 0, at = 10.1, las = 2)



# nutrient n's
intce <- nrow(df_nut) - (sum(df_nut['GCDTypeCa-liming + Wood ash']) + sum(df_nut['GCDTypeCompost']) + sum(df_nut['GCDTypeSludge (including Biosolids)'])+
                           sum(df_nut['GCDTypeManure + Slurry']) + +sum(df_nut['GCDTypeResidue + Mulch']) +
                           sum(df_nut['GCDTypeOther Organic fertilisers (NOT including compost and Urea)']) +
                           sum(df_nut['GCDTypeMixture']))

mtext(paste0("n = ", intce, " (88)"), side = 4, line = 0, at = 8.1, las = 2)

mtext(paste0("n = ", sum(df_nut['GCDTypeCa-liming + Wood ash']), " (14)"), side = 4, line = 0, at = 7.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeCompost']), " (11)"), side = 4, line = 0, at = 6.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeSludge (including Biosolids)']), " (10)"), side = 4, line = 0, at = 5.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeManure + Slurry']), " (46)"), side = 4, line = 0, at = 4.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeResidue + Mulch']), " (32)"), side = 4, line = 0, at = 3.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeOther Organic fertilisers (NOT including compost and Urea)']), " (15)"), side = 4, line = 0, at = 2.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeMixture']), " (8)"), side = 4, line = 0, at = 1.1, las = 2)



## after I've made the plot, repeat the plot over the grey shading
polygon(x=c(-1.5,-1.5,1.4,1.4), y=c(18.5,22.5,22.5,18.5), col="#DCDCDC7D", border=F)
polygon(x=c(-1.5,-1.5,1.4,1.4), y=c(9.5,11.5,11.5,9.5), col="#DCDCDC7D", border=F)

par(new = TRUE)

errbar(x =all_slabs, y = panel_dat$pred, 
       yplus = panel_dat$ci.ub,
       yminus=panel_dat$ci.lb, cap=0.2,
       cex = 2)
abline(v=0, lty =2)


## Add in microfauna for lui

points(y = 13:17+0.1, x = luidat$pred[micro], col = "darkgrey", pch = 19, cex = 2)
segments(luidat$ci.lb[micro], 13:17+0.1, luidat$ci.ub[micro], 13:17+0.1, col = "darkgrey")

legend(0.55,14.5, legend=c("Macro/Meso-fauna", "Micro-fauna"), pch = 19,
       col = c("black", "darkgrey"), cex = 1, bty = "n", pt.cex = 1.5)



## labelling
mtext("Land-use\nintensification", side = 2, line = 5.5, at = 15, cex = 1.2)
mtext("Climate\nchange", side = 2, line = 5.5, at = 20.5, cex = 1.2)
mtext("Pollution", side = 2, line = 5.5, at = 10.5, cex = 1.2)
mtext("Nutrient\nenrichment", side = 2, line = 5.5, at = 4.5, cex = 1.2)

mtext("(a)", side = 2, line = 8, at = 22.2, cex = 1.2, las = 2)
mtext("(b)", side = 2, line = 8, at = 17.2, cex = 1.2, las = 2)
mtext("(c)", side = 2, line = 8, at = 12.2, cex = 1.2, las = 2)
mtext("(d)", side = 2, line = 8, at = 8.2, cex = 1.2, las = 2)

# remove not needed ticks on y
segments(-2, 9, -1.588, 9, col = "white", lwd  = 1.3)
segments(-2, 12, -1.588, 12, col = "white", lwd  = 1.3)
segments(-2, 18, -1.588, 18, col = "white", lwd = 1.3)

# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(22.5,25,25,22.5), col="white", border=F)



dev.off()



## -------------------------------------------------

## Taxa graph

mod.gsba.taxa <- readRDS(file = "Models/GSBAMod_june2023.rds")

taxa_dat <- allDat[allDat$GSBA %in% c("Acari", "Collembola",  "Earthworms", "Nematodes"),]



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
                                             c(0,0,1,0,0, 0,1,0,  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
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

gcds <- rep(c("Climate change", "habitat", "invasives", "LUI", "Nutrient", "Pollution"), each = 4)
taxa <- rep(c("Acari","Collembola", "Earthworms", "Nematodes"), times = 6)

# This is because I can't get the left margin to be bigger...annoying
taxaspace <- paste("                          ", taxa)

# exclude habitat frag and invasives (little data)
subs <- rev(c(1:4, 13:16, 21:24, 17:20))

# for number of cases
table(taxa_dat$driver, taxa_dat$GSBA)
## for number of articles
with(taxa_dat, tapply(ID, list(GSBA, driver), function(x) length(unique(x))))
# I regret not soft-coding this




#jpeg(filename = file.path(FigoutPath, "GSBAMod.jpg"),
#     width = 10, height = 6, units = "in", res = 600)
pdf(file = file.path(FigoutPath, "GSBAMod_2023.pdf"),
       width = 10, height = 6)



par(mar=c(5, 10, 2, 8))

errbar(x =taxaspace[subs], y = t_dat$pred[subs], 
       yplus = t_dat$ci.ub[subs],
       yminus=t_dat$ci.lb[subs], 
       cex = 2,
       xaxt="n")

mtext("Effect size", side = 1, line = 3, cex = 1.5)

mtext("Climate\nchange", side = 2, line = 4.5, at = 14.5, cex = 1.2)
mtext("Land-use\nintensification", side = 2, line = 4.5, at = 10.5, cex = 1.2)
mtext("Pollution\n", side = 2, line = 4.5, at = 6.5, cex = 1.2)
mtext("Nutrient\nenrichment", side = 2, line = 4.5, at = 2.5, cex = 1.2)




polygon(x=c(-1.5,-1.5,1.4,1.4), y=c(0.5,4.5,4.5,0.5), col="#DCDCDC7D", border=F)
polygon(x=c(-1.5,-1.5,1.4,1.4), y=c(8.5,12.5,12.5,8.5), col="#DCDCDC7D", border=F)

par(new = TRUE)
errbar(x =taxaspace[subs], y = t_dat$pred[subs], 
       yplus = t_dat$ci.ub[subs],
       yminus=t_dat$ci.lb[subs], 
       cex = 2,
       xaxt="n")

abline(v=0, lty =2)


# climate chate
mtext("n = 146 (33)", side = 4, line = 0, at = 16.1, las = 2)
mtext("n = 81 (38)", side = 4, line = 0, at = 15.1, las = 2)
mtext("n = 22 (8)", side = 4, line = 0, at = 14.1, las = 2)
mtext("n = 109 (41)", side = 4, line = 0, at = 13.1, las = 2)
# LUI
mtext("n = 136 (51)", side = 4, line = 0, at = 12.1, las = 2)
mtext("n = 82 (49)", side = 4, line = 0, at = 11.1, las = 2)
mtext("n = 174 (76)", side = 4, line = 0, at = 10.1, las = 2)
mtext("n = 187 (84)", side = 4, line = 0, at = 9.1, las = 2)
# pollution
mtext("n = 105 (33)", side = 4, line = 0, at = 8.1, las = 2)
mtext("n = 99 (41)", side = 4, line = 0, at = 7.1, las = 2)
mtext("n = 98 (33)", side = 4, line = 0, at = 6.1, las = 2)
mtext("n = 218 (63)", side = 4, line = 0, at = 5.1, las = 2)
# nutrient enrichment
mtext("n = 181 (40)", side = 4, line = 0, at = 4.1, las = 2)
mtext("n = 90 (40)", side = 4, line = 0, at = 3.1, las = 2)
mtext("n = 136 (41)", side = 4, line = 0, at = 2.1, las = 2)
mtext("n = 152 (66)", side = 4, line = 0, at = 1.1, las = 2)


# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(16.5,22,22,16.5), col="white", border=F)




dev.off()


##### ----------------------------------------------------------

## PUBLICATION BIAS PLOTS ---------------

invas.mod.3 <- readRDS(file = "Models/invasiveMod.rds")





jpeg(filename = file.path(FigoutPath, "FunnelPlots.jpg"),
     width = 13, height = 8, units = "in", res = 600)


par(mfrow = c(2,3), mar = c(2, 5, 0, 0))
funnel(mod.2, main="", xlab = "", ylab = "")
corner.label2(label = "(a)", cex = 2)


par(mar = c(2, 2, 0, 0))
funnel(climate.mod.5,  main = "Climate change", xlab = "", ylab = "")
corner.label2(label = "(b) Climate change", cex = 2)


funnel(lui.mod.2,  main = "", xlab = "", ylab = "")
corner.label2(label = "(c) LUI", cex = 2)


par(mar = c(5, 5, 0, 0))
funnel(poll.mod.22,  main = "",xlab = "", ylab = "")
corner.label2(label = "(d) Pollution", cex = 2)

par(mar = c(5, 2, 0, 0))
funnel(nutri.mod.2, main = "", xlab = "", ylab = "")
corner.label2(label = "(e) Nutrient enrichment", cex = 2)


funnel(invas.mod.3,  main = "", xlab = "", ylab = "")
corner.label2(label = "(f) Invasive species", cex = 2)

mtext("Standard Error", side = 2,  outer = TRUE, line = -2, cex =2)
mtext("Residual Value", side = 1,  outer = TRUE, line = -2, cex =2)

dev.off()




#### --------------------------------------------------
## MEASUREMENT TYPE
measurement.mod.1 <- readRDS("Models/MeasurementMod_June2023.rds")
df1 <- as.data.frame(measurement.mod.1$X)

testdf <- as.data.frame(measurement.mod.1$mf.r)


measurement_dat <- allDat[which(allDat$UniqueID %in% testdf$UniqueID),]


measurementdat <-predict(measurement.mod.1, newmods=rbind(c(0,0,0),
                                                          c(1,0,0),
                                                          c(0,1,0),
                                                          c(0,0,1)
), addx=TRUE, digits=2) #

slabs <- c("Abundance", "Biomass", "Richness", "Shannon")

tapply(measurement_dat$ID, measurement_dat$Measurement, function(x) length(unique(x)))


pdf(file = file.path(FigoutPath, "MeasurementMod.pdf"),
    width = 10, height = 4)

par(mar=c(4, 8, 1.5, 7))
errbar(x =slabs, y = measurementdat$pred, 
       yplus = measurementdat$ci.ub,
       yminus=measurementdat$ci.lb, cap=0.2,
       cex = 2, ylim = c(-1,0))
abline(v=0, lty =2)
mtext("Effect size", side = 1, line = 3, cex = 1.5)
mtext(paste0("n = ", sum(df1['MeasurementShannon']), " (123)"), side = 4, line = 0, at = 4.1, las = 2)
mtext(paste0("n = ", sum(df1['MeasurementRichness']), " (153)"), side = 4, line = 0, at = 3.1, las = 2)
mtext(paste0("n = ", sum(df1$MeasurementBiomass), " (124)"), side = 4, line = 0, at = 2.1, las = 2)
intce <- nrow(df1) - (sum(df1['MeasurementBiomass']) + sum(df1['MeasurementRichness']) +  sum(df1$MeasurementRichness ))
mtext(paste0("n = ", intce, " (574)"), side = 4, line = 0, at = 1.1, las = 2)

dev.off()


#
#### --------------------------------------------------


## Looking at the raw values of pollution
poll <- readRDS("Models/pollutionDataFrame.rds")


poll_cut <- poll[which(poll$effect > -20),]
poll_cut <- poll_cut[which(poll_cut$effect < 20),]


poll_cut$GCDTypeSource <- as.factor(poll_cut$GCDTypeSource)
poll_cut$GCDTypeSource <- factor(poll_cut$GCDTypeSource, 
                                 levels = c("Metals - Mining/Smelting", "Metals - Waste/sewage", "Metals - Urban/transport",
                                            "Metals - Industrial", "Metals - Others", "Pesticides - Farming", "Pesticides - Others"))


jpeg(filename = file.path(FigoutPath, "PollutionRaw.jpg"),
     width = 12, height = 6, units = "in", res = 600)

par(mar=c(5, 5, 1, 1))

poll_labels <- c("Metals\nMining/Smelting", "Metals\nWaste/sewage", "Metals\nUrban/transport",
                 "Metals\nIndustrial", "Metals\nOthers", "Pesticides\nFarming", "Pesticides\nOthers")


boxplot(poll_cut$effect ~ poll_cut$GCDTypeSource,
        ylab = "", xlab = "", 
        axes=FALSE, outline=FALSE, ylim = c(-16,15))
stripchart(poll_cut$effect ~ poll_cut$GCDTypeSource,
           vertical=TRUE,
           method = "jitter",
           pch = 19,
           col = 1:7,
           add = TRUE)



axis(1, at = 1:7, labels = poll_labels, padj = 1)
axis(2)
mtext(side = 2, text = "Effect size", cex = 1.5, line =3)

dev.off()

#### --------------------------------------------------



## STUFF FOR TABLES

# Missing stressors
table(allDat$GCDType, allDat$driver)

## Taxanomic harmonisation table


## Small clean up to make my life easier
# put all first words with capitals
allDat$TaxaGroup2 <- paste(toupper(substr(allDat$TaxaGroup, 1, 1)), 
              substr(allDat$TaxaGroup, 2, nchar(allDat$TaxaGroup)), sep="")
# replace & with 'and'
allDat$TaxaGroup2 <- gsub("&", "and", allDat$TaxaGroup2)

taxaharmonisation <- as.data.frame(table(allDat$TaxaGroup2, allDat$GSBA))
taxaharmonisation <-taxaharmonisation[which(taxaharmonisation$Freq > 0),]


write.csv(taxaharmonisation, "TaxaHarmonisation Table_June2023Revision.csv", row.names = FALSE)


table(allDat$GSBA, allDat$Body.Size)

table(allDat$TaxaGroup[which(allDat$GSBA == "Insects")], allDat$Body.Size[which(allDat$GSBA == "Insects")])
table(allDat$TaxaGroup[which(allDat$GSBA == "Poly-groups: Arachnida")], allDat$Body.Size[which(allDat$GSBA == "Poly-groups: Arachnida")])

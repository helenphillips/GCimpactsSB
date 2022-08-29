
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

# devtools::install_github("MathiasHarrer/dmetar")
library(dmetar) # For 3-level I^2

allDat <- read.csv("Data/03_Data/HedgesData_cleaned.csv")
meta <- read.csv("Data/February2022/processed/metadata.csv")


FigoutPath <- "C:/Users/helenp/WORK/GCimpactsSB/Figures"




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


jpeg(filename = file.path(FigoutPath, "GlobalMap.jpg"),
     width = 10, height = 6, units = "in", res = 600)

# pdf(file = file.path(FigoutPath, "GlobalMap.pdf"),
#     width = 10, height = 6)




map <- ggplot(wm)+
  geom_sf(aes(fill = (no.stud.cat)))+
  scale_fill_viridis_d(option = "plasma", 
                       begin = 0,
                       end = 0.9,
                       na.value = "gray",
                       name = "Number of\npublications")+
  theme_bw()+
  theme(legend.position = c(0,0),
        # legend.position = "left",
        legend.justification = c(-0.2, -0.1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=14))

map

dev.off()



## -----------------------------------------------------
## Frequency plot of GCs and fauna size


taxadat <- allDat
taxadat$Body.Size <- as.factor(taxadat$Body.Size)
taxadat$driver <- as.factor(taxadat$driver)

taxadat$Body.Size <- factor(taxadat$Body.Size, levels = c("All sizes", "Micro-fauna", "Meso-fauna", "Macro-fauna"))
taxadat$driver <- factor(taxadat$driver, levels = c("Climate" ,"LUI" , "Pollution", "NutrientEnrichment", 
                                                    "Invasives",  "HabitatLoss"))

jpeg(filename = file.path(FigoutPath, "BodySizexGCD.jpg"),
     width = 10, height = 6, units = "in", res = 600)

barplot(table(taxadat$Body.Size, taxadat$driver),
        names.arg = c("Climate change", "Land-use \nintensification", "Pollution", "Nutrient \nenrichment", "Invasive \nspecies", "Habitat \nfragmentation"),
        xlab = "",
        ylab = "Number of cases",
        axes = TRUE)
legend("topright", legend = c("All sizes", "Micro-fauna", "Meso-fauna", "Macro-fauna"),
       fill = gray.colors(4), bty = "n")


dev.off()

## ------------------------------------------------------


## MODEL FIGURES
# Main mod

mod.2 <- readRDS(file = "Models/MainMod.rds")

mainmodpred <-predict(mod.2, newmods=rbind(c(0,0,0,0,0), 
                                           c(1,0,0,0,0),
                                           c(0,1,0,0,0),
                                           c(0,0,1,0,0),
                                           c(0,0,0,1,0),
                                           c(0,0,0,0,1)
), addx=TRUE, digits=2) #

slabs <- c("Climate Change", "Habitat\nfragmentation", "Invasive\nspecies", "Land-use\nintensification", "Nutrient\nenrichment", "Pollution")

order <- rev(c(1, 4, 6, 5, 3, 2)) # reverse of the actual order


# how many studies per GC
tapply(allDat$ID, allDat$driver, function(x) length(unique(x)))



jpeg(filename = file.path(FigoutPath, "DriversMod.jpg"),
     width = 10, height = 6, units = "in", res = 600)
#pdf(file = file.path(FigoutPath, "DriversMod.pdf"),
#     width = 10, height = 6)

par(mar=c(5, 10, 1, 7))
errbar(x =slabs[order], y = mainmodpred$pred[order], 
       yplus = mainmodpred$ci.ub[order],
       yminus=mainmodpred$ci.lb[order], cap=0.2,
       cex = 2)
abline(v=0, lty =2)
mtext("Effect size", side = 1, line = 3, cex = 1.5)

mtext("n = 443 (92)", side = 4, line = 0, at = 6.1, las = 2)
mtext("n = 877 (240)", side = 4, line = 0, at = 5.1, las = 2)
mtext("n = 797 (151)", side = 4, line = 0, at = 4.1, las = 2)
mtext("n = 786 (158)", side = 4, line = 0, at = 3.1, las = 2)
mtext("n = 186 (36)", side = 4, line = 0, at = 2.1, las = 2)
mtext("n = 100 (22)", side = 4, line = 0, at = 1.1, las = 2)

# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(6.5,10,10,6.5), col="white", border=F)

dev.off()


## ---------------------------------------------------------------------
## Panel plot for three stressor graphs

## Climate change 

climate.mod.5 <- readRDS(file = "Models/ClimateMod.rds")
nrow(climate.mod.5$X) # this is the dataframe (ish)
df_climate <- as.data.frame(climate.mod.5$X)

climatedat <-predict(climate.mod.5, newmods=rbind(c(0,0,0),
                                                  c(1,0,0),
                                                  c(0,1,0),
                                                  c(0,0,1)
), addx=TRUE, digits=2) #

climate_slabs <- c("Gas change", "Temperature change", "Drought", "Flooding")


## lui

lui.mod.2 <- readRDS(file = "Models/LUIMod.rds")
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

lui_slabs <- rev(c("Grazing", "Fire", "Harvesting", "Inorganic", "Tillage"))


# just micro
df_luimicro <- df_lui[which(df_lui$`Body.SizeMicro-fauna` == 1),]

# not micro
df_luinotmicro <- df_lui[which(df_lui$`Body.SizeMicro-fauna` == 0),]






## nutrient
nutri.mod.2 <- readRDS(file = "Models/nutriMod.rds")

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


nut_slabs <- c("Synthetic Fertilizers", "Ca-liming + Wood ash", "Compost", "Manure + Slurry", "                        Multiple fertilizer types", 
               "Other Organic fertilisers", "Residue + Mulch", "Sludge")


nut_ord <- rev(c(1, 2, 3,8, 4,  7, 6, 5))




## create one big dataframe

panel_dat <- rbind(as.data.frame(nutridat)[nut_ord,1:6], 
                   rep(NA, 6), # sneaky way to get a gap in the plot
                   as.data.frame(luidat)[macro,1:6],
                   rep(NA, 6),
                   as.data.frame(climatedat)[,1:6] )


panel_dat$driver <- c(rep("nutri", length(nut_slabs)), NA, rep("lui", length(lui_slabs)), NA, rep("climate", times = 4) )


all_slabs <- c(nut_slabs[nut_ord], "", lui_slabs, "", climate_slabs)




jpeg(filename = file.path(FigoutPath, "stressorsPanel.jpg"),
     width = 12, height = 8, units = "in", res = 600)
# pdf(file = file.path(FigoutPath, "stressorsPanel.pdf"),
#    width = 12, height = 8)


par(mar=c(5, 10, 2, 6))
errbar(x =all_slabs, y = panel_dat$pred, 
       yplus = panel_dat$ci.ub,
       yminus=panel_dat$ci.lb, cap=0.2,
       cex = 2)
abline(v=0, lty =2)


mtext("Effect size", side = 1, line = 3, cex = 1.5)
mtext(paste0("n = ", sum(df_climate['GCDTypeWaterAvailability-Flood'])), side = 4, line = 0, at = 19.1, las = 2)
mtext(paste0("n = ", sum(df_climate['GCDTypeWaterAvailability-Drought'])), side = 4, line = 0, at = 18.1, las = 2)
mtext(paste0("n = ", sum(df_climate$GCDTypeTemperature )), side = 4, line = 0, at = 17.1, las = 2)
intce <- nrow(df_climate) - (sum(df_climate['GCDTypeWaterAvailability-Flood']) + sum(df_climate['GCDTypeWaterAvailability-Drought']) +  sum(df_climate$GCDTypeTemperature ))
mtext(paste0("n = ", intce), side = 4, line = 0, at = 16.1, las = 2)



mtext(paste0("n = ", sum(df_lui['GCDTypeFire'])), side = 4, line = 0, at = 13.1, las = 2)
mtext(paste0("n = ", sum(df_lui['GCDTypeHarvesting'])), side = 4, line = 0, at = 12.1, las = 2)
mtext(paste0("n = ", sum(df_lui['GCDTypeOrganic versus Inorganic'])), side = 4, line = 0, at = 11.1, las = 2)
mtext(paste0("n = ", sum(df_lui['GCDTypeTillage'])), side = 4, line = 0, at = 10.1, las = 2)

intce <- nrow(df_lui) - (sum(df_lui['GCDTypeFire']) + sum(df_lui['GCDTypeHarvesting']) + 
                           sum(df_lui['GCDTypeOrganic versus Inorganic']) +
                           sum(df_lui['GCDTypeTillage']))
mtext(paste0("n = ", intce), side = 4, line = 0, at = 14.1, las = 2)





intce <- nrow(df_nut) - (sum(df_nut['GCDTypeCa-liming + Wood ash']) + sum(df_nut['GCDTypeCompost']) + sum(df_nut['GCDTypeSludge (including Biosolids)'])+
                           sum(df_nut['GCDTypeManure + Slurry']) + +sum(df_nut['GCDTypeResidue + Mulch']) +
                           sum(df_nut['GCDTypeOther Organic fertilisers (NOT including compost and Urea)']) +
                           sum(df_nut['GCDTypeMixture']))

mtext(paste0("n = ", intce), side = 4, line = 0, at = 8.1, las = 2)

mtext(paste0("n = ", sum(df_nut['GCDTypeCa-liming + Wood ash'])), side = 4, line = 0, at = 7.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeCompost'])), side = 4, line = 0, at = 6.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeSludge (including Biosolids)'])), side = 4, line = 0, at = 5.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeManure + Slurry'])), side = 4, line = 0, at = 4.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeResidue + Mulch'])), side = 4, line = 0, at = 3.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeOther Organic fertilisers (NOT including compost and Urea)'])), side = 4, line = 0, at = 2.1, las = 2)
mtext(paste0("n = ", sum(df_nut['GCDTypeMixture'])), side = 4, line = 0, at = 1.1, las = 2)



## after I've made the plot, repeat the plot over the grey shading

polygon(x=c(-1.5,-1.5,1.4,1.4), y=c(0.5,8.5,8.5,0.5), col="#DCDCDC7D", border=F)
polygon(x=c(-1.5,-1.5,1.4,1.4), y=c(15.5,19.5,19.5,15.5), col="#DCDCDC7D", border=F)

par(new = TRUE)

errbar(x =all_slabs, y = panel_dat$pred, 
       yplus = panel_dat$ci.ub,
       yminus=panel_dat$ci.lb, cap=0.2,
       cex = 2)
abline(v=0, lty =2)


## Add in microfauna for lui

points(y = 10:14+0.1, x = luidat$pred[micro], col = "darkgrey", pch = 19, cex = 2)
segments(luidat$ci.lb[micro], 10:14+0.1, luidat$ci.ub[micro], 10:14+0.1, col = "darkgrey")

legend(0.55,14.5, legend=c("Macro/Meso-fauna", "Micro-fauna"), pch = 19,
       col = c("black", "darkgrey"), cex = 1, bty = "n", pt.cex = 1.5)



## labelling
mtext("Land-use\nintensification", side = 2, line = 5.5, at = 12, cex = 1.2)
mtext("Climate\nchange", side = 2, line = 5.5, at = 17.5, cex = 1.2)
mtext("Nutrient\nenrichment", side = 2, line = 5.5, at = 4.5, cex = 1.2)

mtext("(a)", side = 2, line = 8, at = 19.2, cex = 1.2, las = 2)
mtext("(b)", side = 2, line = 8, at = 14.2, cex = 1.2, las = 2)
mtext("(c)", side = 2, line = 8, at = 8.2, cex = 1.2, las = 2)

# remove not needed ticks on y
segments(-2, 9, -1.565, 9, col = "white")
segments(-2, 15, -1.565, 15, col = "white")

# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(19.5,22,22,19.5), col="white", border=F)



dev.off()



## -------------------------------------------------

## Taxa graph

mod.gsba.taxa <- readRDS(file = "Models/GSBAMod.rds")

taxa_dat <- allDat[allDat$GSBA %in% c("Acari", "Collembola",  "Earthworms", "Nematodes"),]
table(taxa_dat$driver, taxa_dat$GSBA)

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

gcds <- rep(c("Climate change", "habitat", "invasives", "LUI", "Nutrient", "Pollution"), each = 4)
taxa <- rep(c("Acari","Collembola", "Earthworms", "Nematodes"), times = 6)

taxaspace <- paste("                          ", taxa)

# exclude habitat frag and invasives (little data)
subs <- rev(c(1:4, 13:16, 21:24, 17:20))



jpeg(filename = file.path(FigoutPath, "GSBAMod.jpg"),
     width = 10, height = 6, units = "in", res = 600)
# pdf(file = file.path(FigoutPath, "GSBAMod.pdf"),
#       width = 10, height = 6)



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



mtext("n = 141", side = 4, line = 0, at = 16.1, las = 2)
mtext("n = 78", side = 4, line = 0, at = 15.1, las = 2)
mtext("n = 22", side = 4, line = 0, at = 14.1, las = 2)
mtext("n = 112", side = 4, line = 0, at = 13.1, las = 2)
mtext("n = 138", side = 4, line = 0, at = 12.1, las = 2)
mtext("n = 84", side = 4, line = 0, at = 11.1, las = 2)
mtext("n = 174", side = 4, line = 0, at = 10.1, las = 2)
mtext("n = 188", side = 4, line = 0, at = 9.1, las = 2)
mtext("n = 105", side = 4, line = 0, at = 8.1, las = 2)
mtext("n = 99", side = 4, line = 0, at = 7.1, las = 2)
mtext("n = 98", side = 4, line = 0, at = 6.1, las = 2)
mtext("n = 218", side = 4, line = 0, at = 5.1, las = 2)
mtext("n = 180", side = 4, line = 0, at = 4.1, las = 2)
mtext("n = 90", side = 4, line = 0, at = 3.1, las = 2)
mtext("n = 136", side = 4, line = 0, at = 2.1, las = 2)
mtext("n = 160", side = 4, line = 0, at = 1.1, las = 2)


# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(16.5,22,22,16.5), col="white", border=F)




dev.off()


##### ----------------------------------------------------------

## PUBLICATION BIAS PLOTS ---------------

invas.mod.3 <- readRDS(file = "Models/invasiveMod.rds")
poll.mod.1d <- readRDS(file = "Models/pollutionMod.rds")





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
funnel(poll.mod.1d,  main = "",xlab = "", ylab = "")
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
measurement.mod.1 <- readRDS("Models/MeasurementMod.rds")
df1 <- as.data.frame(measurement.mod.1$X)

measurementdat <-predict(measurement.mod.1, newmods=rbind(c(0,0,0),
                                                          c(1,0,0),
                                                          c(0,1,0),
                                                          c(0,0,1)
), addx=TRUE, digits=2) #

slabs <- c("Abundance", "Biomass", "Richness", "Shannon")



jpeg(filename = file.path(FigoutPath, "MeasurementMod.jpg"),
     width = 10, height = 4, units = "in", res = 600)


par(mar=c(4, 8, 1.5, 7))
errbar(x =slabs, y = measurementdat$pred, 
       yplus = measurementdat$ci.ub,
       yminus=measurementdat$ci.lb, cap=0.2,
       cex = 2, ylim = c(-1,0))
abline(v=0, lty =2)
mtext("Effect size", side = 1, line = 3, cex = 1.5)
mtext(paste0("n = ", sum(df1['MeasurementBiomass'])), side = 4, line = 0, at = 4.1, las = 2)
mtext(paste0("n = ", sum(df1['MeasurementRichness'])), side = 4, line = 0, at = 3.1, las = 2)
mtext(paste0("n = ", sum(df1$MeasurementRichness)), side = 4, line = 0, at = 2.1, las = 2)
intce <- nrow(df1) - (sum(df1['MeasurementBiomass']) + sum(df1['MeasurementRichness']) +  sum(df1$MeasurementRichness ))
mtext(paste0("n = ", intce), side = 4, line = 0, at = 1.1, las = 2)

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

boxplot(poll_cut$effect ~ poll_cut$GCDTypeSource,
        ylab = "Effect size", xlab = "", 
        names = "")
stripchart(poll_cut$effect ~ poll_cut$GCDTypeSource,
           vertical=TRUE,
           method = "jitter",
           pch = 19,
           col = 1:7,
           add = TRUE)

poll_labels <- c("Metals\nMining/Smelting", "Metals\nWaste/sewage", "Metals\nUrban/transport",
                 "Metals\nIndustrial", "Metals\nOthers", "Pesticides\nFarming", "Pesticides\nOthers")

axis(1, at = 1:7, labels = poll_labels, padj = 1)
axis(2)
mtext(side = 2, text = "Effect size", cex = 1.5, line =3)

dev.off()

#### --------------------------------------------------



## STUFF FOR TABLES

# Missing stressors
table(allDat$GCDType, allDat$driver)


taxaharmonisation <- as.data.frame(table(allDat$TaxaGroup, allDat$GSBA))
taxaharmonisation <-taxaharmonisation[which(taxaharmonisation$Freq > 0),]
write.csv(taxaharmonisation, "TaxaHarmonisation Table.csv", row.names = FALSE)


table(allDat$GSBA, allDat$Body.Size)

table(allDat$TaxaGroup[which(allDat$GSBA == "Insects")], allDat$Body.Size[which(allDat$GSBA == "Insects")])
table(allDat$TaxaGroup[which(allDat$GSBA == "Poly-groups: Arachnida")], allDat$Body.Size[which(allDat$GSBA == "Poly-groups: Arachnida")])

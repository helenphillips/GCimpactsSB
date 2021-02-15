
setwd("~/WORK/GCimpactsSB")

### LOAD LIBRARIES ------------------------

library(ggplot2)
library(metafor)

### LOAD DATA -----------------------------

dat_frag <- read.csv("Data/DataExtractionTable - HabitatLoss.csv")
dat_cc <- read.csv("Data/DataExtractionTable - ClimateChange.csv")
dat_inv <- read.csv("Data/DataExtractionTable - Invasives.csv")
dat_lui <- read.csv("Data/DataExtractionTable - LandUseIntensification.csv")
dat_lui_r <- read.csv("Data/DataExtractionTable_Rowan - LandUseIntensification.csv")
dat_ne <- read.csv("Data/DataExtractionTable - NutrientEnrichment.csv")
dat_ne_y <- read.csv("Data/DataExtractionTable - NutrientEnrichment_yiming.csv")
dat_poll <- read.csv("Data/DataExtractionTable - Pollution.csv")


### REMOVE COLS ------------------------------

neededCols <- c("ID","Case_ID","CaseDescription", "System" ,  
                "Control_mean",
                "Control_SD","Control_N","Treatment_mean","Treatment_SD"  ,         
                "Treatment_N","Measurement","MeasurementUnits","Error",
               "TaxaGroup","TaxaBodysize" )


dat_frag <- dat_frag[,which(names(dat_frag) %in% neededCols)]
dat_cc <- dat_cc[,which(names(dat_cc) %in% neededCols)]
dat_inv <- dat_inv[,which(names(dat_inv) %in% neededCols)]
dat_lui <- dat_lui[,which(names(dat_lui) %in% neededCols)]
dat_lui_r <- dat_lui_r[,which(names(dat_lui_r) %in% neededCols)]
dat_ne <- dat_ne[,which(names(dat_ne) %in% neededCols)]
dat_ne_y <- dat_ne_y[,which(names(dat_ne_y) %in% neededCols)]
dat_poll <- dat_poll[,which(names(dat_poll) %in% neededCols)]

## ADD GCD COL ---------------


dat_frag$GCD <- 'fragmentation'
dat_cc$GCD <- 'climatechange'
dat_inv$GCD <- 'invasives'
dat_lui$GCD <- 'lui'
dat_lui_r$GCD <- 'lui'
dat_ne$GCD <- 'nutrient'
dat_ne_y$GCD <- 'nutrient'

dat_poll$GCD <- 'pollution'


##JOIN TOGETHER -------------------

dat <- rbind(dat_frag, dat_cc, dat_inv,dat_lui,dat_lui_r,
             dat_ne, dat_ne_y, dat_poll)


dat$Measurement <- as.factor(dat$Measurement)
dat$GCD <- as.factor(dat$GCD)


### REMOVE ROWS WITH NO ERRORS -----------------

dat <- dat[!(is.na(dat$Control_SD)),]
dat <- droplevels(dat)


### REMOVE ROWS WITH NO PAPER ID -----------------
dat <- dat[!(is.na(dat$ID)),]
dat <- droplevels(dat)

###### CONVERSION OF ERROR -----------------------
# Should all be SDs

table(dat$Error)

## SEs
ses <- which(dat$Error == "SE")

dat$Control_SD[ses] <- dat$Control_SD[ses] * sqrt(dat$Control_N[ses])
dat$Treatment_SD[ses] <- dat$Treatment_SD[ses] * sqrt(dat$Treatment_N[ses])
dat$Error[ses] <- "SD"

## CI
## https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm#:~:text=The%20standard%20deviation%20for%20each,should%20be%20replaced%20by%205.15.
ci <- which(dat$Error == "CI")

dat$Control_SD[ci] <- sqrt(dat$Control_N[ci]) * ((2 * dat$Control_SD[ci])/3.92)
dat$Treatment_SD[ci] <- sqrt(dat$Treatment_N[ci]) * ((2 * dat$Treatment_SD[ci])/3.92)
dat$Error[ci] <- "SD"


ci <- which(dat$Error == "CI95")

dat$Control_SD[ci] <- sqrt(dat$Control_N[ci]) * ((2 * dat$Control_SD[ci])/3.92)
dat$Treatment_SD[ci] <- sqrt(dat$Treatment_N[ci]) * ((2 * dat$Treatment_SD[ci])/3.92)
dat$Error[ci] <- "SD"

## Get rid of CV and DT
dat <- dat[-which(dat$Error == "CV"),]
dat <- dat[-which(dat$Error == "DT"),]
dat <- droplevels(dat)
## unknowns marked as SD (the largest error)

others <- which(dat$Error != "SD")
dat$Error[others] <- "SD"


#### CALCULATE LOG-RESPONSE --------------------------

lrr_dat<-
  escalc(measure="ROM", ## log-response ratio
         m1i=Treatment_mean, ## treatment mean
         m2i=Control_mean, ## control mean
         sd1i=Treatment_SD, ## standard deviations control
         sd2i=Control_SD,
         n1i=Treatment_N,
         n2i=Control_N, 
         data=dat,
         var.names=c("LRR","LRR_var"),
         digits=4)




hist(lrr_dat$LRR)
hist(lrr_dat$LRR_var)



## ONLY ABUNDANCE DATA FOR NOW --------


lrr_dat <- lrr_dat[which(lrr_dat$Measurement == "Abundance"),]

## remove NAs
lrr_dat <- lrr_dat[!(is.na(lrr_dat$LRR_var)),]

## THREE RIDICULOUS VARIANCES

lrr_dat <- lrr_dat[which(lrr_dat$LRR_var < 50000),]
plot(lrr_dat$LRR_var)

lrr_dat <- lrr_dat[which(lrr_dat$LRR_var < 1000),]
plot(lrr_dat$LRR_var)

lrr_dat <- lrr_dat[which(lrr_dat$LRR_var < 400),]
plot(lrr_dat$LRR_var)

## add small amount to all vars

lrr_dat$LRR_var <-lrr_dat$LRR_var + 0.5

## MODEL ------
# mod.int<-rma.mv(
#   yi=LRR,
#   V=LRR_var, 
#   mods=~1,
#   random= ~1|ID,
#   struct="CS",
#   method="ML",
#   digits=4,
#   data=lrr_dat)



mod.1<-rma.mv(
  yi=LRR,
  V=LRR_var, 
  mods=~GCD - 1,
  random= ~1|ID,
  struct="CS",
  method="ML",
  digits=4,
  data=lrr_dat)



qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
## Could be better


summary(mod.1)



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

jpeg(filename = "GCD_carleton.jpeg", quality = 100, res = 300, width = 3500, height = 2000)

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

dev.off()


## FRAGMENTATION ANALYSIS

frag_dat <- lrr_dat[lrr_dat$GCD == "fragmentation",]
frag.mod<-rma.mv(
  yi=LRR,
  V=LRR_var, 
  mods=~GCD - 1,
  random= ~1|ID,
  struct="CS",
  method="ML",
  digits=4,
  data=lrr_dat)







library(dplyr)
# load packages for maps
library(maptools)
library(ggmap)
library(rnaturalearth)
library(countrycode)
library(viridis)
library(ggplot2)
library(patchwork)
library(reshape2)


meta <- read.csv("Data/DataExtractionTable - Metadata.csv")
meta_r <- read.csv("Data/DataExtractionTable_Rowan - Metadata.csv")
meta_y <- read.csv("Data/DataExtractionTable_yiming - Metadata.csv")


metaCols <- c("ID", "Year", "Country")

meta <- meta[,which(names(meta) %in% metaCols)]
meta_r <- meta_r[,which(names(meta_r) %in% metaCols)]
meta_y <- meta_y[,which(names(meta_y) %in% metaCols)]


allmeta <- rbind(meta, meta_r, meta_y)
unique(allmeta$Country)[order(unique(allmeta$Country))]


allmeta$Country[which(allmeta$Country == "Czech Rpublic" | allmeta$Country == "Czechia"  )] <- "Czech Republic"


allmeta$Country[which(allmeta$Country == "England" 
                      | allmeta$Country == "Northern Ireland"
                      | allmeta$Country == "Scotland"
                      | allmeta$Country == "the UK"
                      | allmeta$Country == "United Kingdom"
                      | allmeta$Country == "Wales")] <- "UK"

allmeta$Country[which(allmeta$Country == "United States" 
                      | allmeta$Country == "the US")] <- "USA"

allmeta$Country[which(allmeta$Country == "The Netherlands")] <- "Netherlands"
allmeta$Country[which(allmeta$Country == "CÃ´te d'Ivoire")] <- "Ivory Coast"
allmeta$Country[which(allmeta$Country == "Falkland Islands/Antarctica")] <- "Falkland Islands"


allmeta <- allmeta[-which(allmeta$Country %in% c("Africa",
"Brazil/Uruguay",
"France and Italy"    ,                        
"France, The Netherlands, Switzerland, Canada",
"Ireland, Poland",
"Netherlands, England, Portugal, Germany",
"Portugal,Brazil"   ,
"The Netherlands, UK" ,
"UK, Netherlands, Denmark" )),]

allmeta <- allmeta[-which(is.na(allmeta$ID)),]

stud <- allmeta %>% 
  group_by(Country) %>% 
  summarize(no.stu = n())


# harmonize country names
stud <- stud %>%
  mutate(iso =  countrycode(stud$Country, origin = "country.name", destination = "iso3c")) %>%
  dplyr::group_by(iso) %>%
  summarize(no.stud = sum(no.stu, na.rm = TRUE))


wm <- ne_countries(scale = 110, returnclass = "sf") %>%
  left_join(stud, by = c("adm0_a3" = "iso"))

wm <- sf::st_transform(wm,"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

mybreaks <- c(5, 25, 50, 75, 100)

map <- ggplot(wm)+
  geom_sf(aes(fill = (no.stud)))+
  scale_fill_viridis(option = "viridis", 
                     begin = 0,
                     end = 1,
                     na.value = "gray", breaks = mybreaks,
                     name = "Number of\nstudies")+
  theme_bw()+
  theme(legend.position = c(0,0),
        # legend.position = "left",
        legend.justification = c(-0.2, -0.1),
        legend.title = element_text(size=15),
        legend.text = element_text(size=14))

map


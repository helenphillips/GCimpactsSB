## Fixing the coordinates in the file

setwd("C:/Users/helenp/WORK/GCimpactsSB")
## LOAD THE DATA
dataDir <- "Data/September2022"
dat <- read.csv(file = file.path(dataDir, "processed", "alldata.csv"))


dat$UniqueID <- paste(dat$ID, dat$Case_ID, dat$driver)


# One row has latitude but no longitude was available
dat$Latitude[which(dat$ID == 169)] <- NA


# remember this will overwrite, so need to start with longitude
dat$Longitude[which(dat$Latitude == "Devecser, Hungary")] <-  "17.450361"
dat$Latitude[which(dat$Latitude == "Devecser, Hungary")] <-  "47.106098"

dat$Longitude[which(dat$Latitude == "Kosogorsky near Tula, Russia")] <-  "37.628746"
dat$Latitude[which(dat$Latitude == "Kosogorsky near Tula, Russia")] <-  "54.188831"

dat$Longitude[which(dat$Latitude == "Shenyang Experimental Station of Ecology")] <-  "123.374861"
dat$Latitude[which(dat$Latitude == "Shenyang Experimental Station of Ecology")] <-  "41.518705"

dat$Longitude[which(dat$Latitude == "Derio, Spain")] <-  "-2.888318"
dat$Latitude[which(dat$Latitude == "Derio, Spain")] <-  "43.295085"
                            
dat$Longitude[which(dat$Latitude == "Ghent (nearby)")] <-  "3.748882"
dat$Latitude[which(dat$Latitude == "Ghent (nearby)")] <-  "51.040514"

dat$Longitude[which(dat$Latitude == "Biesbosch")] <-  "4.763134"
dat$Latitude[which(dat$Latitude == "Biesbosch")] <-  "51.750642"

dat$Longitude[which(dat$Latitude == "Atcala de Henares")] <-  "-3.371655"
dat$Latitude[which(dat$Latitude == "Atcala de Henares")] <-  "40.488516"

dat$Longitude[which(dat$Latitude == "Naples")] <-  "14.406152"
dat$Latitude[which(dat$Latitude == "Naples")] <-  "40.795076"

dat$Longitude[which(dat$Latitude == "Krompachy, Slovakia")] <-  "20.879764"
dat$Latitude[which(dat$Latitude == "Krompachy, Slovakia")] <-  "48.905454"

dat$Longitude[which(dat$Latitude == "Donana Ntional Park, Spain")] <-  "-6.434663"
dat$Latitude[which(dat$Latitude == "Donana Ntional Park, Spain")] <-  "37.042276"

dat$Longitude[which(dat$Latitude == "Mount Paddusstieva")] <-  "18.816727"
dat$Latitude[which(dat$Latitude == "Mount Paddusstieva")] <-  "68.354176"
# Abisko Scientific Research Station, based on paper info
                                                                                                     
dat$Longitude[which(dat$Latitude == "Mount Slattatjakka")] <-  "18.668255"
dat$Latitude[which(dat$Latitude == "Mount Slattatjakka")] <-  "68.332905"

dat$Longitude[which(dat$Latitude == "Tianjin Normal University campus.")] <-  "117.183297"
dat$Latitude[which(dat$Latitude == "Tianjin Normal University campus.")] <-  "39.095519"

dat$Longitude[which(dat$Latitude == "Mortagne-du-Nord, France")] <-  "3.454121"
dat$Latitude[which(dat$Latitude == "Mortagne-du-Nord, France")] <-  "50.498901"

dat$Longitude[which(dat$Latitude == "Research Institute for Soil Science and Agricultural Chemistry (RISSAC) of the Hungarian Academy of Sciences at Nagyhorcsok, Hungary")] <-  "19.059678"
dat$Latitude[which(dat$Latitude == "Research Institute for Soil Science and Agricultural Chemistry (RISSAC) of the Hungarian Academy of Sciences at Nagyhorcsok, Hungary")] <-  "47.474162"

dat$Longitude[which(dat$Latitude == "Copper smelting works, Glogow, poland")] <-  "16.046837"
dat$Latitude[which(dat$Latitude == "Copper smelting works, Glogow, poland")] <-  "51.668480"

dat$Longitude[which(dat$Latitude == "Hygum site, Jutland")] <-  "9.303022"
dat$Latitude[which(dat$Latitude == "Hygum site, Jutland")] <-  "56.021977"

dat$Longitude[which(dat$Latitude == "KoilaKozanis")] <- "21.790661"
dat$Latitude[which(dat$Latitude == "KoilaKozanis")] <- "40.330928"
 
# "germany" # not doing this one, because the sites are all over germany
dat$Latitude[which(dat$Latitude == "germany")] <- ""




dat$Longitude[which(dat$Latitude == "48.43922897879104, 18.917915872750914")] <- "18.917915872750914" 
dat$Latitude[which(dat$Latitude == "48.43922897879104, 18.917915872750914")] <- "48.43922897879104"


## Is it Brazil, so is south and west
dat$Latitude[which(dat$Latitude == "-22d10m13.53s")] <- "22d10m13.53sS"
dat$Longitude[which(dat$Latitude == "-47d53m58.12s")] <- "47d53m58.12sW"



dat$Latitude[which(dat$Latitude == "40d-43mN")] <- "40d43mN"

# Coordinates make no sense (in paper). But gave location as Guianga, Tugbok District, Davao City
dat$Longitude[which(dat$Latitude == "7006m29.83sN")] <- "125.507742"
dat$Latitude[which(dat$Latitude == "7006m29.83sN")] <- "7.151770"


#UTM? Named as Micheville, France
dat$Longitude[which(dat$Latitude == "6897475")] <- "5.880795"
dat$Latitude[which(dat$Latitude == "6897475")] <- "49.460592"


dat$Latitude[which(dat$Latitude == "50.3000dN")] <- "50.3000" 
dat$Longitude[which(dat$Longitude == "19.9833dE")] <- "19.9833" 

dat$Latitude[which(dat$Latitude == "9d34mE")] <- "56d29mN"
dat$Longitude[which(dat$Longitude == "56d29mN")] <- "9d34mE"


# Not sure what I/they were going for. Longitude was correct
dat$Latitude[which(dat$Latitude == "22\024d36mN")] <- "22d36mN"

# unique(dat$Longitude)

dat$Latitude[which(dat$Latitude == "40d14m46.5066s")] <- "40d14m46.5066sN"
dat$Longitude[which(dat$Longitude == "-8d20m23.9964s")] <- "8d20m23.9964sW"

dat$Longitude[which(dat$Longitude == "-47d53m58.12s")] <- "47d53m58.12sW"

dat$Longitude[which(dat$Longitude == "48d35m49sN")] <- "20d13m40sE"


dat$Latitude[which(dat$Latitude == "43d20l13sN")] <- "43d20m13sN"

## Fixing based on the conversion from DMS to decimal degress ----

dat$Latitude[which(dat$Latitude == "37.7833S")] <- "-37.7833"
dat$Longitude[which(dat$Longitude == "144.9667E")] <- "144.9667"

dat$Latitude[which(dat$Latitude == "49.5N")] <- "49.5"
dat$Longitude[which(dat$Longitude == "119.31W")] <- "-119.31"


dat$Latitude[which(dat$Latitude == "22s14m00sS")] <- "22d14m00sS"


# incorrect format, but when I checked the paper, it would be better to round to 9d10
dat$Longitude[which(dat$Longitude == "834mE")] <- "9d10mE"


dat$Longitude[which(dat$Longitude == "16d19m26.18'sE")] <- "16d19m26.18sE"

dat$Latitude[which(dat$Latitude == "42d22N")] <- "42d22mN"


## fixing based on matching coordinates to countries

dat$Longitude[which(dat$ID == "163")] <- "76d07m00sW" #authors got it wrong
dat$Longitude[which(dat$ID == "305")] <- "123d30m3.384sE" #authors got it wrong

dat$Latitude[which(dat$ID == "330")] <-"39.051467"
dat$Longitude[which(dat$ID == "330")] <-"22.961876" # coordinate was just off coast, this is coordinate for region

dat$Latitude[which(dat$ID == "357")] <-"-77.633333" #who ever did the original data entry missed the minus

dat$Latitude[which(dat$ID == "464")] <- "51d04mN"
dat$Longitude[which(dat$ID == "464")] <- "02d06mW" # there were two different coordinates, and one was definitely in the sea

dat$Latitude[which(dat$ID == "509")] <- "49.884490"
dat$Longitude[which(dat$ID == "509")] <- "3.009181" #unclear/incorrect format in paper, this is the lcoation of the town

dat$Longitude[which(dat$ID == "525")] <- "84.33" #who ever did the original data entry missed the minus

dat$Latitude[which(dat$ID == "525")] <- "52.3" 
dat$Longitude[which(dat$ID == "525")] <- "10.43" # data entry mixed up lat and long 

dat$Latitude[which(dat$ID == "525")] <- "-45.02" 
dat$Longitude[which(dat$ID == "525")] <- "170.76" # wasn't in dms 

dat$Latitude[which(dat$ID == "570")] <- "49d14m48sN" # wrong in paper
dat$Latitude[which(dat$ID == "571")] <- "48d55m20sN" # wrong in paper

dat$Latitude[which(dat$ID == "628")] <- "68.337168"
dat$Longitude[which(dat$ID == "628")] <- "18.848387" #paper was off

dat$Latitude[which(dat$ID == "913")] <- "53d28m43sN"
dat$Longitude[which(dat$ID == "913")] <- "6d14m06sE" # changing to what was in the paper for ease

dat$Latitude[which(dat$ID == "977")] <- "42d24mN"
dat$Longitude[which(dat$ID == "977")] <- "128d06mE" # using the plot level coordinates rather than general area

dat$Latitude[which(dat$ID %in% c("982", "983", "984"))] <- "52.3"
dat$Longitude[which(dat$ID %in% c("982", "983", "984"))] <- "10.433" # wrong way round

dat$Latitude[which(dat$ID == "986")] <- "56d23mN"
dat$Longitude[which(dat$ID == "986")] <- "10d57mE"

dat$Latitude[which(dat$ID == "1041")] <- "71d52m40sS"
dat$Longitude[which(dat$ID == "1041")] <- "68d15m57sW" #missed the minus. Just put back into presented coordinates

dat$Latitude[which(dat$ID == "1049")] <- "47.732275"
dat$Longitude[which(dat$ID == "1049")] <- "16.773987" # coordinates didn't match the name given (or country)


dat$Latitude[which(dat$ID == "1050")] <- "48d21mN"
dat$Longitude[which(dat$ID == "1050")] <- "85d21mW" # missedminus, but also put back into presented coordiantes

dat$Latitude[which(dat$ID == "1062")] <- "54d47mS"
dat$Longitude[which(dat$ID == "1062")] <- "68d16mW" # missedminus, but also put back into presented coordiantes

dat$Latitude[which(dat$ID == "1107")] <- "31d28mN"
dat$Longitude[which(dat$ID == "1107")] <- "121d56mE" # missedminus, but also put back into presented coordiantes

dat$Latitude[which(dat$ID == "1119")] <- "33d08mN"
dat$Longitude[which(dat$ID == "1119")] <- "117d21mE" # latitude put into long, but also put back into presented coordiantes

dat$Latitude[which(dat$ID == "1199")] <- "54.662927"
dat$Longitude[which(dat$ID == "1199")] <- "-2.253358" # coordinates didn't match location name in paper

dat$Latitude[which(dat$ID == "1265")] <- "52.140"
dat$Longitude[which(dat$ID == "1265")] <- "4.838" # put into wrong format during extraction

dat$Latitude[which(dat$ID == "1267")] <- "60.989078"
dat$Longitude[which(dat$ID == "1267")] <- "25.513673" # coordinates didn't match location name in paper


dat$Latitude[which(dat$ID == "1282")] <- "54d51mS"
dat$Longitude[which(dat$ID == "1282")] <- "68d36mW" # pmissedminus, but also put back into presented coordiantes

dat$Latitude[which(dat$UniqueID  %in% c("1708 a Invasives", "1708 c Invasives"))] <- "47d16m0.00sN"
dat$Longitude[which(dat$UniqueID  %in% c("1708 a Invasives", "1708 c Invasives"))] <- "94d23m48.60sW"
dat$Latitude[which(dat$UniqueID  %in% c("1708 b Invasives", "1708 d Invasives"))] <- "46d26m3.06sN"
dat$Longitude[which(dat$UniqueID  %in% c("1708 b Invasives", "1708 d Invasives"))] <- "91d19m36.00sW"

dat$Latitude[which(dat$ID == "1711")] <- "39.076854" # 
dat$Longitude[which(dat$ID == "1711")] <- "-86.555186" # incorrect sign (extracted incorrectly), but coordinates didn't match location given in paper

dat$Latitude[which(dat$ID == "1766")] <- "37d24mN" # 
dat$Longitude[which(dat$ID == "1766")] <- "122d13mW"  # pmissedminus, but also put back into presented coordiantes

dat$Latitude[which(dat$ID == "1822")] <- "42d34m20sN" 
dat$Longitude[which(dat$ID == "1822")] <- "9d2m6sW"  # seconds missed off (put back into presented coordinated format)


dat$Latitude[which(dat$ID == "2005")] <- "42d34m20sN" 
dat$Longitude[which(dat$ID == "2005")] <- "9d2m6sW"  # seconds missed off (put back into presented coordinated format)


dat$Latitude[which(dat$UniqueID 
                   %in% c("2005 e Climate", "2005 f Climate", "2005 g Climate", "2005 h Climate"))] <- "-60.71"
dat$Longitude[which(dat$UniqueID
                    %in% c("2005 e Climate", "2005 f Climate", "2005 g Climate", "2005 h Climate"))] <- "-45.59"

dat$Latitude[which(dat$UniqueID 
                   %in% c("2005 i Climate", "2005 j Climate", "2005 k Climate", "2005 l Climate"))] <- "-67.61"
dat$Longitude[which(dat$UniqueID
                    %in% c("2005 i Climate", "2005 j Climate", "2005 k Climate", "2005 l Climate"))] <- "-68.22"

dat$Latitude[which(dat$UniqueID 
                   %in% c("2005 a Climate", "2005 b Climate", "2005 c Climate", "2005 d Climate"))] <- "-51.76"
dat$Longitude[which(dat$UniqueID
                    %in% c("2005 a Climate", "2005 b Climate", "2005 c Climate", "2005 d Climate"))] <- "-59.03" # missed signs


dat$Latitude[which(dat$ID == "2195")] <- "49d25m5sN" # 
dat$Longitude[which(dat$ID == "2195")] <- "1d7m41sE"  # wrong way round


dat$Latitude[which(dat$ID == "2286")] <- "24d32mN" # 
dat$Longitude[which(dat$ID == "2286")] <- "91d15mE"  # their figure didn't match coordinates given. Changed coordinates based on figure


dat$Latitude[which(dat$ID == "2320")] <- "51d24m44sN" # 
dat$Longitude[which(dat$ID == "2320")] <- "05d02m45sE"  # their figure didn't match coordinates given. Changed coordinates based on figure



dat$Latitude[which(dat$UniqueID  == "2375 a Climate")] <- "52d20m3sN"
dat$Longitude[which(dat$UniqueID  == "2375 a Climate")] <- "6d27m27sW"  # just switching back to original
dat$Latitude[which(dat$UniqueID  == "2375 b Climate")] <- "51d15m0sN"
dat$Longitude[which(dat$UniqueID  == "2375 b Climate")] <- "22d34m0sE"  # just switching back to original

dat$Longitude[which(dat$ID == "2448")] <- "128d16.3mE" # typo

dat$Latitude[which(dat$ID == "2541")] <- "46d50mS" # wrong direction
dat$Longitude[which(dat$ID == "2541")] <- "37d50mE"  # also maybe a typo in both

dat$Latitude[which(dat$ID == "2734")] <- "14d04mN" # wrong direction

dat$Latitude[which(dat$ID == "2779")] <- "60d43mS" # wrong direction
dat$Longitude[which(dat$ID == "2779")] <- "45d36mW" 

dat$Latitude[which(dat$ID == "2858")] <- "35d54mN" # wrong direction
dat$Longitude[which(dat$ID == "2858")] <- "84d21mW" 

dat$Latitude[which(dat$ID == "2918")] <- "61.010101" # coordiantes didbn't match locaiton givne
dat$Longitude[which(dat$ID == "2918")] <- "25.385887" 

dat$Latitude[which(dat$ID == "2934")] <- "32d35m5sN" # switched lat and long (also put back into orgiinal)
dat$Longitude[which(dat$ID == "2934")] <- "119d42m0sE" 

dat$Latitude[which(dat$ID == "2960")] <- "45d34mN" # mised minus
dat$Longitude[which(dat$ID == "2960")] <- "84d40mW" 


dat$Latitude[which(dat$ID == "2972")] <- "37d30m50sN" # authors must have switched them
dat$Longitude[which(dat$ID == "2972")] <- "6d13m22sE" 

dat$Latitude[which(dat$ID == "3020")] <- "64d00m01sN" # missed minus
dat$Longitude[which(dat$ID == "3020")] <- "21d11m09sW" 

dat$Latitude[which(dat$ID == "3031")] <- "46d50m50sN" # transcribed wrong way
dat$Longitude[which(dat$ID == "3031")] <- "6d34m30sE" 

dat$Longitude[which(dat$ID == "3089")] <- "-96.811" # missed minus sign 

dat$Latitude[which(dat$ID == "3098")] <- "60.344895" # coordinates didn't match location given
dat$Longitude[which(dat$ID == "3098")] <- "17.549498" 

dat$Longitude[which(dat$ID == "3132")] <- "119d28mE" # typo

dat$Latitude[which(dat$ID == "3165")] <- "32d46mS" # coordinates didn't match location given
dat$Longitude[which(dat$ID == "3165")] <- "152d05mE" 

dat$Latitude[which(dat$ID == "3172")] <- "56d33mN" # typo?
dat$Longitude[which(dat$ID == "3172")] <- "13d13mE" 

dat$Latitude[which(dat$ID == "3199")] <- "32d35mN" # typo?
dat$Longitude[which(dat$ID == "3199")] <- "119d42mE" 


dat$Longitude[which(dat$ID == "1814")] <- "175d50mE"
dat$Longitude[which(dat$ID == "1729")] <- "175d25mE"


dat$Latitude[which(dat$ID == "532")] <- "69" # making into decimal degrees fully
dat$Longitude[which(dat$ID == "532")] <- "88" 

dat$Latitude[which(dat$ID == "690")] <- "44.65138" # making into decimal degrees fully
dat$Longitude[which(dat$ID == "690")] <- "10.08330" 

dat$Latitude[which(dat$ID == "803")] <- "45.6" # making into decimal degrees fully
dat$Longitude[which(dat$ID == "803")] <- "-89.4" 

dat$Latitude[which(dat$ID == "64")] <- NA # coordinates wrong in paper. not able to figure out
dat$Longitude[which(dat$ID == "64")] <- NA 

dat$Latitude[which(dat$ID == "2916")] <- "-30" # putting into dd fully
dat$Longitude[which(dat$ID == "2916")] <- "150" 

dat$Latitude[which(dat$ID == "2993")] <- "59" # putting into dd fully
dat$Longitude[which(dat$ID == "2993")] <- "18" 

dat$Latitude[which(dat$ID == "2217")] <- "16" # putting into dd fully
dat$Longitude[which(dat$ID == "2217")] <- "-61" 


dat$Latitude[which(dat$ID == "782")] <- "35.80" # putting into dd fully
dat$Longitude[which(dat$ID == "782")] <- "-93.50" 

dat$Latitude[which(dat$ID == "744")] <- "40d49mN" # putting into dd fully
dat$Longitude[which(dat$ID == "744")] <- "107d46mW" 

dat$Latitude[which(dat$ID == "2535")] <- "51.966135" # putting into dd fully
dat$Longitude[which(dat$ID == "2535")] <- "20.163874" 

dat$Latitude[which(dat$ID == "44")] <- "49" # putting into dd fully
dat$Longitude[which(dat$ID == "44")] <- "8" 

dat$Latitude[which(dat$ID == "432")] <- "51.813" # putting into dd fully
dat$Longitude[which(dat$ID == "432")] <- "0.381" 

dat$Latitude[which(dat$ID == "835")] <- "53d10mN" # missing directions
dat$Longitude[which(dat$ID == "835")] <- "5d51mE" 

dat$Latitude[which(dat$ID == "544")] <- "52d18mN" # missing directions
dat$Longitude[which(dat$ID == "544")] <- "10d26mE" 




### CREATE SIMPLE DATAFRAME


datgps <- dat[,c("ID", "UniqueID", "Latitude" , "Longitude")]

#### CLEANING -----
# first remove white space

datgps$Latitude <- gsub(" ", "", datgps$Latitude)
datgps$Longitude <- gsub(" ", "", datgps$Longitude)


# create new columns
datgps$lat_dd <- NA
datgps$long_dd <- NA




# if already in decimal degrees, can copy over as a number
# first get vector in dms

indms_lat <- grep("[A-z]+", datgps$Latitude, value = FALSE, perl = TRUE)
indms_long <- grep("[A-z]+", datgps$Longitude, value = FALSE, perl = TRUE)


# then 'not' those
datgps$lat_dd[-(indms_lat)] <- as.numeric(datgps$Latitude[-(indms_lat)])
datgps$long_dd[-(indms_long)] <- as.numeric(datgps$Longitude[-(indms_long)])

# Then hopefully convert the rest
library(sp)

  
datgps$lat_dd[indms_lat] <- as.numeric(char2dms(datgps$Latitude[indms_lat], chd='d', chm='m', chs='s'))
datgps$long_dd[indms_long] <- as.numeric(char2dms(datgps$Longitude[indms_long], chd='d', chm='m', chs='s'))


# row for each unique combo
coord <- datgps[!duplicated(datgps[,c(1, 5:6)]),]


# need to add in NAs to blanks
coord$Latitude[which(coord$Latitude == "")] <- NA
coord$Longitude[which(coord$Longitude == "")] <- NA

## check NAs. 
which(is.na(coord$lat_dd) & !(is.na(coord$Latitude)))
which(is.na(coord$long_dd) & !(is.na(coord$Longitude))) # same rows

coord[which(is.na(coord$lat_dd) & !(is.na(coord$Latitude))),]




coord <- coord[,c(1, 6, 5)]

coord<- coord[which(!(is.na(coord$lat_dd))),]
coord<- coord[which(!(is.na(coord$long_dd))),]


#coord$X<-coord$Group.1
# coord<-coord[1:3]
names(coord)<-c("ID", "Long", "Lat")



dsSPDF<-SpatialPointsDataFrame(coord[,2:3], data.frame(coord[,2:3]))
proj4string(dsSPDF)<-CRS("+proj=longlat")

library(maps)

map("world",border="gray87",fill=TRUE, col="gray87")

points(dsSPDF, col="black", bg="black", cex= 1, pch=19)

## All seem fine now


##  (the code below was used earlier in the process as a check)

## SAVING -------------


write.csv(datgps, file = "Data/CoordinatesFixed.csv", row.names = FALSE)






#### CHECKING ----------------



library(sp)
library(rworldmap)
library(maptools)
# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
coords2country = function(points)
{
  # prepare a SpatialPolygons object with one poly per country
  countries = map('world', fill=TRUE, col="transparent", plot=FALSE)
  names = sapply(strsplit(countries$names, ":"), function(x) x[1])
  
  # clean up polygons that are out of bounds
  filter = countries$x < -180 & !is.na(countries$x)
  countries$x[filter] = -180
  
  filter = countries$x > 180 & !is.na(countries$x)
  countries$x[filter] = 180
  
  IDs <- sapply(strsplit(countries$names, ":"), function(x) x[1])
  countriesSP = map2SpatialPolygons(countries, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  
  # convert our list of points to a SpatialPoints object
  pointsSP = SpatialPoints(points, proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  
  # Return the state names of the Polygons object containing each point
  myNames = sapply(countriesSP@polygons, function(x) x@ID)
  myNames[indices]
}

coord2 <- coord
coord2$country  <- coords2country(coord[,2:3])



### BRING IN META DATA TO CHECK AGAINST -----

meta <- read.csv("Data/February2022/processed/metadata.csv")


meta_short <- meta[,c("ID", "Country")]

coord2 <- merge(coord2, meta_short, by.x = "X", by.y = "ID", all.x = TRUE)

write.csv(coord2, file = "Data/checkingCoordinates.csv", row.names = FALSE)




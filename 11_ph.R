## Fixing the ph values

## when a study gave ph values for control and treatment sites, we recorded both
## with a semi-colon between the two numbers

## need to process them (maybe take the mean?)


setwd("C:/Users/helenp/WORK/GCimpactsSB")
## LOAD THE DATA
dataDir <- "Data/September2022"
dat <- read.csv(file = file.path(dataDir, "processed", "alldata.csv"))


dat$UniqueID <- paste(dat$ID, dat$Case_ID, dat$driver)


## also seems at least one as a - 

grep(";", dat$SoilPH)
grep("-", dat$SoilPH)

## could also check by looking for ones which have more that 3 digits
checking1 <- nchar(gsub("[^0-9]+", "", dat$SoilPH)) ## remove non-digits then count
checking1 <- which(checking1 > 3)
checking2 <- grep(";", dat$SoilPH)

setdiff(checking1, checking2) # which in checking1 is not in checking 2
# i.e., which has more than 3 digits, but doesn't have a semi-colon

dat[c(371, 1651, 1652, 2531),] # the one -, two have 3 digits after decimal point, and one with a comma between the two values



#### SUBSET TO MAKE EASIER -----------------

phs <- dat[,c('ID', 'UniqueID', 'SoilPH')]

### MAKING ALL SEMI-COLONS ---------------
phs$SoilPH <- gsub("-", ";", phs$SoilPH)
phs$SoilPH <- gsub(",", ";", phs$SoilPH)

## ID THE ONES TO FIX ---------------

tofix <- grep(";", phs$SoilPH)

phs$ph_fixed <- phs$SoilPH
## FIX THEM ---------------------------

foo <- data.frame(do.call('rbind', strsplit(as.character(phs$SoilPH[tofix]),';',fixed=TRUE)))

foo$X2[which(foo$X2 == "NA")] <- NA # to prevent the warning showing later
foo$X1 <- as.numeric(foo$X1)
foo$X2 <- as.numeric(foo$X2)

newphs <- rowMeans(foo[,c('X1', 'X2')], na.rm=TRUE)



phs$ph_fixed[tofix] <- newphs

### SAVE ----------------

phs$SoilPH <- NULL


write.csv(phs, file = "Data/Phs_Fixed.csv", row.names = FALSE)






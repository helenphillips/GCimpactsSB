## SetWD

## LIBRARIES -----------------------------------------------

# setwd("~/WORK/GCimpactsSB")
setwd("C:/Users/helenp/WORK/GCimpactsSB")


## data dir ------
dataDir <- "Data/January2022"




## Metadata files first
meta <- read.csv(file.path(dataDir, "metadata.csv"))


## filter out the fully factorial ones
ff_dat <- meta[which(meta$FullyFactorial == TRUE),]
summary(ff_dat)
ff_dat$totalDrivers <- rowSums(ff_dat[ , 14:19], na.rm=TRUE)


## STORAGE -------------



res <- data.frame(matrix(NA, nrow = 6, ncol = 6))
colnames(res) <- names(ff_dat[,14:19])
rownames(res) <- names(ff_dat[,14:19])



drivers <-  names(ff_dat[,14:19])

for(d in 1:length(drivers)){

  drive <- drivers[d]
  tmp_drivers <- drivers[-d]
## cases --------------------------
  temp_dat  <- ff_dat[which(ff_dat[which(names(ff_dat) == drive)] == TRUE),] # 34
# how many are only lui
  res[d, which(names(res) == drive)] <- length(which(temp_dat$totalDrivers == 1))

   # length(which(lui$totalDrivers > 2)) # 3

  for(dd in tmp_drivers){
    res[d, which(names(res) == dd)] <- 
      length(which(temp_dat[which(names(ff_dat) == dd)]  == TRUE))
    
    
  }
  
}







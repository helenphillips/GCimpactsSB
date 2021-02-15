
metafor_PlotAllCovs <- function(mod, dat, covariate1, covariate2, col1 ='black', pch1 = 19){

  library(Hmisc)
  
  print("Make sure all covariates are either 'factors' or 'continuous' in the dataframe given!")
  print("Covariates will be plotted in order provided")
  
  print("This function currently assumes that there are two covariates AND that both covariates are factors")
  print("if you need other variations of this, tell Helen to stop being lazy")
  
  # Check that the covariates are in the dataframe
if(!(all(c(covariate1, covariate2) %in% names(dat)))){
  stop("All covariates must be columns in the dataframe given")}

  col_cov1 <- which(names(dat) == covariate1)
  col_cov2 <- which(names(dat) == covariate2)
  
  newdat <- expand.grid(cov1=levels(dat[,col_cov1]),
                        cov2=levels(dat[,col_cov2]))
  
  names(newdat) <- c(parse(text=covariate1), parse(text=covariate2))
  
  newdat <- model.matrix(~ get(covariate1) + get(covariate2), data=newdat)
  
  
  attributes(newdat)$dimnames[[2]] <- gsub('get(covariate1)', covariate1, attributes(newdat)$dimnames[[2]], fixed = TRUE)
  attributes(newdat)$dimnames[[2]] <- gsub('get(covariate2)', covariate2, attributes(newdat)$dimnames[[2]], fixed = TRUE)
  

  df <- predict(mod, newmods=newdat[,-1], addx=TRUE)
  



nlevels1 <- length(grep(covariate1, dimnames(df$X)[[2]])) + 1
colscov1 <- c(1, grep(covariate1, dimnames(df$X)[[2]])) # 1 for the intercept

nlevels2 <- length(grep(covariate2, dimnames(df$X)[[2]])) + 1

y.label <- levels(dat[,col_cov1])




df <- as.data.frame(df)


chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
secondfactor <- chunk(1:nrow(df),nlevels2) # x = vector, n = number in chunk

## Make the dataframe nicer to work with
df$V1 <- NA
df$V2 <- NA

for(t in 1:length(secondfactor)){
  df$V2[secondfactor[[t]]] <- levels(dat[,col_cov2])[t]
  df$V1[secondfactor[[t]]] <- levels(dat[,col_cov1])
}


df[,grep('X.', names(df), fixed = TRUE)] <- NULL

# Order by the first covariate, then the second
df <- df[
  with(df, order(V1, V2)),
]

df$fullname <- paste(df$V1, df$V2)


len_cols <- length(col1)
if(len_cols != 1 && len_cols != nrow(df)){
  stop("Number of colours given needs to either be 1 OR the length of the variable combinations")
}

len_pchs <- length(pch1)
if(len_pchs != 1 && len_pchs != nrow(df)){
  stop("Number of plotting symbols given needs to either be 1 OR the length of the variable combinations")
}


par(mar=c(3, 20, 1, 1))


errbar(df$fullname, df$pred, 
       (df$ci.ub), 
       (df$ci.lb), 
       xlab = "Effect",  
       col = col1, pch = pch1)

ys <-seq(1, nrow(df), by= 1)
points(x = df$pred, y = ys, col=col1, pch = pch1)


hlines <- seq(1, nrow(df), by = nlevels2) - 0.5
hlines <- hlines[-1]

abline(v = 0, lty=2)

 
#mtext(levels(dat$FragmentationDesign)[1], 2, line = 3) 

abline(h=hlines,  col = "grey", lty = 3)
# axis(2, at = hlines, labels = NA, tck = -0.25, lty = 2, outer = TRUE)

 
}
 # offset <- 0

# for(t in 2:length(secondfactor)){
#   
#   offset <- offset + 0.2
#   
#   ys <-seq(1, 5, by= 1) - offset
#   
#   points(x = df$pred[secondfactor[[t]]], y = ys, col='grey', pch = 19)
#   segments(x0=df$ci.lb[secondfactor[[t]]], y0=ys, x1 = df$ci.ub[secondfactor[[t]]], y1 = ys,
#            col = 'grey')
#   
#   
# }
# 






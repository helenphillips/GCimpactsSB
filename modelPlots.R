
metafor_PlotAllCovs <- function(mod, dat, covariate1, covariate2, col1 ='black', pch1 = 19, pchcex = 1,
                                givenLabels = TRUE, newlabels, GreyMissing = TRUE, removeMissing = FALSE){
  
  
  # mod = the model you want to plot
  # dat = the dataframe as specified in the model
  # covariate1 = the first FACTOR covariate (string, so in quotes). The plot will be grouped by this by this one
  # covariate2 = the second FACTOR covariate (string)
  # col1 = either a single colour or a vector of colours. If a vector, has to be as long as covariate1 multiplied by covariate 2
  # pch1 = either a single number for the plotting point type, or a vector. If a vector, has to be as long as covariate1 multiplied by covariate 2
  # pchcex = the size adjustment for the plotting points
  # givenLabels = (TRUE/FALSE). Use the y axis labels (TRUE) as given by the function itself. Recommended on first attempt
  # newlabels = if givenLabels is FALSE, a new vector of labels for the y axis
  # GreyMissing = (TRUE/FALSE) For any points that have no underlying data, make teh points grey
  # removeMissing = (TRUE/FALSE) For any points that have no underlying data, remove them
  
  
  
  library(Hmisc)
  
  print("Make sure all covariates are either 'factors' or 'continuous' in the dataframe given!")
  print("Covariates will be plotted in order provided")
  
  print("This function currently assumes that there are two covariates AND that both covariates are factors")
  print("if you have a continuous variable (AND two factors) use metafor_PlotContCovs function")
  
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
  
  
  ## Descriptors
  df$n <- as.vector(table(dat[,col_cov1], dat[,col_cov2]))
  
  studies <- aggregate(dat$ID, list(dat[,col_cov2], dat[,col_cov1]), function(x) {c(Min = length(unique(x)))})
  studies$fullname <- paste(studies$Group.2, studies$Group.1)
  
  
  
  

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
  df <- merge(df, studies, by = "fullname", all.x = TRUE)
  df$Group.1 <- NULL
  df$Group.2 <- NULL
  names(df)[names(df) == 'x'] <- "nStudies"
  df$nStudies[is.na(df$nStudies)] <- 0
  
  ## A rough way to change the ylabels
  if(!(givenLabels)){
    if(length(newlabels) != nrow(df)){stop("newlabels must be the same length as the dataframe i.e., all possible combinations")}
    df$fullname <- newlabels}
  
  
  len_cols <- length(col1)
  if(len_cols != 1 && len_cols != nrow(df)){
    stop("Number of colours given needs to either be 1 OR the length of the variable combinations")
  }
  
  len_pchs <- length(pch1)
  if(len_pchs != 1 && len_pchs != nrow(df)){
    stop("Number of plotting symbols given needs to either be 1 OR the length of the variable combinations")
  }
  
  
  if(GreyMissing){ # Grey out any factor levels that are missing from underlying data
    
    if(len_cols == 1){ # if colour is only one
      col1 <- rep(col1, times = nrow(df))
    }
    
    zero <- which(df$n == 0)
    col1[zero] <- "grey"
    
  }
  
  
  if(removeMissing){ # remove missing factor levels
    
    zero <- which(df$n == 0)
    df <- df[-zero,]
    
    if(len_cols != 1){ col1 <- col1[-zero]}
    if(len_pchs != 1){ pch1 <- pch1[-zero]}
    
    
  }
  
  
  
  par(mar=c(3, 20, 2, 5))
  
  
  errbar(df$fullname, df$pred, 
         (df$ci.ub), 
         (df$ci.lb), 
         xlab = "Effect",  
         col = col1, pch = pch1)
  
  
  df$allN <- paste0(df$n, "(", df$nStudies, ")")
  
  axis(4, at = 1:length(df$fullname), labels = df$allN, las = 2, tick = FALSE)
  
  ys <-seq(1, nrow(df), by= 1)
  points(x = df$pred, y = ys, col=col1, pch = pch1, cex = pchcex)
  
  
  hlines <- seq(1, nrow(df), by = nlevels2) - 0.5
  hlines <- hlines[-1]
  
  abline(v = 0, lty=2)
  
  
  #mtext(levels(dat$FragmentationDesign)[1], 2, line = 3) 
  if(!(removeMissing)){
    abline(h=hlines,  col = "grey", lty = 3)
    # axis(2, at = hlines, labels = NA, tck = -0.25, lty = 2, outer = TRUE)
  }
  
  ## Add significance star
  Signi <- which(!(df$ci.lb < 0 & df$ci.ub > 0))
  text( df$ci.ub[Signi] + 0.1, Signi, "*", cex = pchcex)
 
  
   
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



metafor_PlotContCovs <- function(mod, dat, contCov, covariate1, covariate2, 
                                 whichCov1lev, whichCov2lev,
                                 col1 ='black', pch1 = 19, pchcex = 1,
                                 ylabel = "Effect Size", xlabel = "Continuous Effect"){
  
  contcov <- which(names(dat) == contCov) # continuous variable column
  col_cov1 <- which(names(dat) == covariate1)
  col_cov2 <- which(names(dat) == covariate2)
  
  newdat <- expand.grid(contcov=seq(min(dat[,contcov]), max(dat[,contcov]), length.out = 100),
                        cov1=levels(dat[,col_cov1]),
                        cov2=levels(dat[,col_cov2]))
  
  names(newdat) <- c(parse(text=contCov), parse(text=covariate1), parse(text=covariate2))
  
  newdat <- model.matrix(~ get(contCov) + get(covariate1) + get(covariate2), data=newdat)
  
  
  attributes(newdat)$dimnames[[2]] <- gsub('get(covariate1)', covariate1, attributes(newdat)$dimnames[[2]], fixed = TRUE)
  attributes(newdat)$dimnames[[2]] <- gsub('get(covariate2)', covariate2, attributes(newdat)$dimnames[[2]], fixed = TRUE)
  attributes(newdat)$dimnames[[2]] <- gsub('get(contCov)', contCov, attributes(newdat)$dimnames[[2]], fixed = TRUE)
  
  
  df <- predict(mod, newmods=newdat[,-1], addx=TRUE)
  
  
  # Col for continuous var
  
  dfcont <- grep(contCov, attributes(df$X)$dimnames[[2]])
  
  
  # Jst the reference for covariate 1
  dfcov1 <- grep(whichCov1lev, attributes(df$X)$dimnames[[2]])
  
  # Just the reference for covariate 2
  dfcov2 <- grep(whichCov2lev, attributes(df$X)$dimnames[[2]])
  
  ncoldf <- 1:ncol(df$X)
  ncoldf <- ncoldf[!(ncoldf %in% c(dfcont, dfcov1, dfcov2))]
  
  df <- df[df$X[,dfcov1] == 1 & df$X[,dfcov2] == 1]
  
  
  for(i in ncoldf){
    print(i)
    df <- df[df$X[,i] == 0]
    
  }
  
  
  
  
  pt_cols <- col1
  
  ci_cols <- adjustcolor(pt_cols, alpha.f = 0.3)
  
  plot(-1e+05, -1e+05, ylim = c(min(df$ci.lb,na.rm = TRUE), max(df$ci.ub, na.rm = TRUE)),
       xlim = c(min(df$X[,dfcont],na.rm = TRUE), 
                max(df$X[,dfcont], na.rm = TRUE)),  ylab = ylabel, xlab = xlabel)
  
  X.Vec <- c(df$X[,dfcont], max(df$X[,dfcont]), 
             rev(df$X[,dfcont]), min(df$X[,dfcont]))
  Y.Vec <- c(df$ci.lb, tail(df$ci.ub, 1), rev(df$ci.ub), df$ci.lb[1])
  
  polygon(X.Vec, Y.Vec, col = ci_cols, border = NA)
  
  points(df$X[,dfcont], df$pred, 
         col = pt_cols,  type = "l", lwd = 5)
  
  
  rug(dat[,contcov])
  
}



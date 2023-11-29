setwd("C:/Users/helenp/WORK/GCimpactsSB")

hedges <- read.csv("Data/03_Data/HedgesData_cleaned_June2023.csv")

library(metafor)
library(tidyverse)


mod2 <- readRDS("Models/MainMod_rerun_June2023.rds")
luimod <- readRDS("Models/LUIMod_redo_june2023.rds")

### INITIAL HELPER FUNCTIONS -----------------


firstup <- function(x, upper = TRUE) {
  if(upper){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  } else{ x }
}

print.orchard <- function(x, ...){
  return(print(x$mod_table))
}

weighted_var <- function(x, weights){
  weight_var <- sum(x * weights) / sum(weights)
  return(weight_var)
}

num_studies <- function(data, mod, group){
  
  # Summarize the number of studies within each level of moderator
  table <- data        %>%
    dplyr::group_by({{mod}}) %>%
    dplyr::summarise(stdy = length(unique({{group}})))
  
  table <- table[!is.na(table$moderator),]
  # Rename, and return
  colnames(table) <- c("Parameter", "Num_Studies")
  return(data.frame(table))
  
}

pred_interval_esmeans <- function(model, mm, mod, ...){
  
  tmp <- summary(mm)
  tmp <- tmp[ , ]
  test.stat <- stats::qt(0.975, tmp$df[[1]])
  
  if(length(model$tau2) <= 1){ # including gamma2
    sigmas <- sum(model$sigma2)
    PI <- test.stat * base::sqrt(tmp$SE^2 + sigmas)
  } else {
    sigmas <- sum(model$sigma2)
    taus   <- model$tau2
    gammas <- model$gamma2
    w_tau <- model$g.levels.k
    w_gamma <- model$g.levels.k
    
    if(mod == "1"){
      tau <- weighted_var(taus, weights = w_tau)
      gamma <- weighted_var(gamma, weights = w_gamma)
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau + gamma)
      
    } else {
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus + gammas)
    }
  }
  
  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI
  
  # renaming "overall" to ""
  if(tmp[1,1] == "overall"){tmp[,1] <- "intrcpt"}
  
  return(tmp)
}
###############################3

##







### modelres code --------------

# I have also added in the raw_data code, to go through that at the point in model res that
# it is called

# I have deleted any code that is not needed (i.e., checks that my code passes, or alternatives
# that my data doesn't need)


# model = mod2
# mod = "driver"
model = lui
mod = "GCDType"
group = "ID"


df_mod = 1.0e6 # almost identical to z value


# Extract the data from the model object
# data <- model$data # this doesn't exist
# data <- hedges # so substituted this in



lui <- hedges[which(hedges$driver == "LUI"),] 
lui <- droplevels(lui[which(lui$Body.Size != "All sizes"),])
lui$GCDType[which(lui$GCDType %in% c("Defoliation"))] <- "Grazing"
lui$GCDType[which(lui$GCDType %in% c("Intensity", "Landscape", "Management", "Weeding", "Planting", "Intensification"))] <- "Management"
lui$GCDType[which(lui$GCDType %in% c("Logging"))] <- "Harvesting"
lui$GCDType[which(lui$GCDType %in% c("Water"))] <- "Irrigation"
lui$GCDType[which(lui$GCDType %in% c("degradation", "Disturbance"))] <- "Degradation"
lui <- lui[which(lui$GCDType != "Human population"),]
lui <- lui[which(lui$GCDType != "Mono- versus poly-culture"),] # too little
lui <- lui[which(lui$GCDType != "Irrigation"),] # too little
lui <- lui[which(lui$GCDType != "Management"),] # too little
lui <- lui[which(lui$GCDType != "Degradation"),] # too little

lui$GCDType <- relevel(as.factor(lui$GCDType), ref = "Grazing")

data <- lui

# Check if missing values exist and use complete case data


at = NULL
by = NULL
weights = NULL
upper = TRUE


# if(is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
  grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = data, coef = model$b,
                        vcov = stats::vcov(model), df = model$k-1) ## NOTE: Added data argument emmeans >vers 1.7.4. Object is unstable so feeding in the relevant arguments from model object directly. Note, we should think about df!
  mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights)
  
  # getting prediction intervals
  mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
  

    mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                            estimate = mm_pi[,"emmean"],
                            lowerCL = mm_pi[,"lower.CL"],
                            upperCL = mm_pi[,"upper.CL"],
                            lowerPR = mm_pi[,"lower.PI"],
                            upperPR = mm_pi[,"upper.PI"])
    
 
  
  # Extract data
  # data2 <- get_data_raw(model, mod, group)
  
  # model, mod, group, N = NULL, at = NULL, subset = TRUE
    N = NULL
    at = NULL
    subset = FALSE
  

  # Extract the data from the model object
  # data <- model$data 


    # Extract effect sizes
    #yi <- model$yi
    #vi <- model$vi
    
    yi <- data$yi
    vi <- data$vi
    
    
    # type <- attr(model$yi, "measure") # not sure what this does, but only length = 1

    # Get moderator
    moderator <- as.character(data[[mod]]) # Could default to base instead of tidy
    moderator <- firstup(moderator) 
  
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  # data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  
  data_reorg <- data.frame(yi, vi, moderator, stdy)
  
  
  #names(data_reorg)[4] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  





  data2 <- data_reorg
  
  


mod_table$name <- factor(mod_table$name,
                         levels = mod_table$name,
                         labels = mod_table$name)


output <- list(mod_table = mod_table,
               data = data2)

class(output) <- c("orchard", "data.frame")




## orchard plot code
# not pasted here

# pasted below so I can edit it
# so run the function below this one first

orchard_plot_ME(output, mod = "driver", group = "ID", xlab = "Effect Size",
            
                                                  transfm = c("none"), 
                                                   legend.pos = c("bottom.right") # "none" - no legends
                                                 )



orchard_plot_ME(output, mod = "GCDType", group = "ID", xlab = "Effect Size",
                
                transfm = c("none"), 
                legend.pos = c("bottom.right") # "none" - no legends
)



#' @title orchard_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, it creates an orchard plot from mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals.
#' @param object model object of class \code{rma.mv}, \code{rma}, or \code{orchard} table of model results.
#' @param mod the name of a moderator. Defaults to \code{"1"} for an intercept-only model. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding varaibles in 'by'. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param xlab The effect size measure label.
#' @param N The name of the column in the data specifying the sample size so that each effect size estimate is scaled to the sample size, N. Defaults to \code{NULL}, so that precision is used for scaling each raw effect size estimate instead of sample size.
#' @param alpha The level of transparency for effect sizes represented in the orchard plot.
#' @param angle The angle of y labels. The default is 90 degrees.
#' @param cb If \code{TRUE}, it uses 20 colour blind friendly colors.
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD). Defaults to \code{"none"}.
#' @param condition.lab Label for the condition being marginalized over.
#' @param trunk.size Size of the mean, or central point.
#' @param branch.size Size of the confidence intervals.
#' @param twig.size Size of the prediction intervals.
#' @param legend.pos Where to place the legend. To remove the legend, use \code{legend.pos = "none"}.
#' @param k.pos Where to put k (number of effect sizes) on the plot. Users can specify the exact position or they can use specify \code{"right"}, \code{"left"},  or \code{"none"}. Note that numeric values (0, 0.5, 1) can also be specified and this would give greater precision.
#' @param colour Colour of effect size shapes. By default, effect sizes are colored according to the \code{mod} argument. If \code{TRUE}, they are colored according to the grouping variable
#' @param fill If \code{TRUE}, effect sizes will be filled with colours. If \code{FALSE}, they will not be filled with colours.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @param flip Logical, defaults to \code{TRUE}, indicating whether the plot should be flipped.
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accounting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1,
#' random=list(~1|ExptID, ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")
#' orchard_plot(results, mod = "Grazer.type",
#' group = "ExptID", xlab = "log(Response ratio) (lnRR)")
#' # or
#' orchard_plot(eklof_MR, mod = "Grazer.type", group = "ExptID",
#' xlab = "log(Response ratio) (lnRR)")
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' orchard_plot(lim_MR, mod = "Phylum", group = "Article",
#' xlab = "Correlation coefficient", transfm = "tanh", N = lim$N)
#' }
#' @export

orchard_plot_ME <- function(object, mod = "1", group, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         trunk.size = 3, branch.size = 1.2, twig.size = 0.5,
                         transfm = c("none", "tanh"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE, flip = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  orchaRd::mod_results(object, mod, group,  N,
                                       by = by, at = at, weights = weights, upper = upper)
    } else {
      results <-  orchaRd::mod_results(object, mod = "1", group,  N,
                                       by = by, at = at, weights = weights, upper = upper)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  data_trim <- results$data
  # making sure factor names match
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision (1/SE)"
  
  if(is.null(N) == FALSE){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
    #latex2exp::TeX()
  }
  
  if(transfm == "tanh"){
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- Zr_to_r(mod_table[,cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }else{
    label <- xlab
  }
  
  # Add in total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))
  
  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])
  
  # the number of groups in a moderator & data points
  group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  # colour blind friendly colours with grey
  cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  # setting fruit colour
  if(colour == TRUE){
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  }else{
    color <- data_trim$mod
    color2 <- mod_table$name
  }
  
  # whether we fill fruit or not
  if(fill == TRUE){
    fill <- color
  }else{
    fill <- NULL
  }
  
  # whether marginal
  if(names(mod_table)[2] == "condition"){
    
    # the number of levels in the condition
    condition_no <- length(unique(mod_table[, "condition"]))
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
                              size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = color2), size = twig.size, position = ggplot2::position_dodge2(width = 0.3), fatten = trunk.size) +
      # this will only work for up to 5 different conditions
      # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
      ggplot2::scale_shape_manual(values =  20 + (1:condition_no))  +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::labs(shape = condition.lab) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle)) 
       
    if(flip){
      plot <- plot + ggplot2::coord_flip()
    }
    
  } else {
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +
      coord_cartesian(xlim = c(-5000, 5000)) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name,  ymin = lowerPR, ymax = upperPR, fill = color2), size = twig.size, fatten = trunk.size, shape = 21) +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle))  
      
    #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
    if(flip){
      plot <- plot + ggplot2::coord_flip()
    }
  }
  
  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting colors in
  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values = cbpl) +
      ggplot2::scale_colour_manual(values = cbpl)
  }
  
  # putting k and g in
  if(k == TRUE && g == FALSE && k.pos == "right"){
    plot <- plot +
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)
  } else if(k == TRUE && g == FALSE && k.pos == "left") {
    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos == "right"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "right", size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos == "left"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "left", size = 3.5)
  } else if (k == TRUE && g == FALSE && k.pos%in%c('right','left','none')==FALSE) {
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
                                                                 1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no]),
                                     parse = TRUE, size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos%in%c('right','left','none')==FALSE) {
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
                                                                 1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no],
                                                                                          "~", "(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, size = 3.5)
  }
  return(plot)
}


#' @title Zr_to_r
#' @description Converts Zr back to r (Pearson's correlation coefficient)
#' @param df data frame of results of class 'orchard'
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

Zr_to_r <- function(df){
  return(sapply(df, tanh))
}




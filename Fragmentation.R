
setwd("~/WORK/GCimpactsSB")

### LOAD LIBRARIES ------------------------

library(ggplot2)
library(metafor)

### LOAD DATA -----------------------------

dat <- read.csv("Data/DataExtractionTable - HabitatLoss.csv")


### REMOVE COLS ------------------------------

neededCols <- c("ID","Case_ID","CaseDescription", "System" ,  
                "FragmentationDesign","Control_mean",
                "Control_SD","Control_N","Treatment_mean","Treatment_SD"  ,         
                "Treatment_N","Measurement","MeasurementUnits","Error",
                "Control_absValue","Control_absValueUnit","Control_description","Treatment_absValue",  
                "Treatment_absValueUnit","Treatment_description","TaxaGroup","TaxaBodysize" )


dat <- dat[,which(names(dat) %in% neededCols)]

### REMOVE ROWS WITH NO ERRORS -----------------

dat <- dat[!(is.na(dat$Control_SD)),]

###### CONVERSION OF ERROR -----------------------
# Should all be SDs

table(dat$Error)

## SEs

ses <- which(dat$Error == "SE")

dat$Control_SD[ses] <- dat$Control_SD[ses] * sqrt(dat$Control_N[ses])
dat$Treatment_SD[ses] <- dat$Treatment_SD[ses] * sqrt(dat$Treatment_N[ses])
dat$Error[ses] <- "SD"

## unknowns marked as SD (the largest error)

others <- which(dat$Error != "SD")
dat$Error[others] <- "SD"

#### CALCULATE LOG-RESPONSE --------------------------

frag_dat<-
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




hist(frag_dat$LRR)
hist(frag_dat$LRR_var)



###### A MODEL ---------------------------

mod.1<-rma.mv(
  yi=LRR,
  V=LRR_var, 
  mods=~FragmentationDesign + Measurement,
  random= ~1|ID,
  struct="CS",
  method="ML",
  digits=4,
  data=frag_dat)



qqnorm(residuals(mod.1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.1,type="pearson"),col="red")
## Could be better


summary(mod.1)


mod.2<-rma.mv(
  yi=LRR,
  V=LRR_var, 
  mods=~FragmentationDesign,
  random= ~1|ID,
  struct="CS",
  method="ML",
  digits=4,
  data=frag_dat)


anova(mod.1, mod.2)


out_int<-rma.mv(
  yi=LRR,
  V=LRR_var, 
  mods=~FragmentationDesign-1,
  random= ~1|ID,
  struct="CS",
  method="REML",
  digits=4,
  data=frag_dat)

y<-summary(out_int)$b
ci_l<-summary(out_int)$ci.lb
ci_h<-summary(out_int)$ci.ub

fg1<-data.frame(cbind(y,ci_l,ci_h))
colnames(fg1)[1]<-"y"
colnames(fg1)[2]<-"ci_l"
colnames(fg1)[3]<-"ci_h"
fg1$Design<-c("Corridors/Connectivity","Edge effects", "Fragmentation per se",
              "Habitat amount", "Isolation")
fg1$Design<-as.factor(fg1$Design)

fg1


pdf('Fragmentation_mod.pdf')

p <- ggplot(fg1, aes(x=Design, y=y, ymin=ci_l, ymax=ci_h))+
  geom_pointrange()+
  geom_hline(yintercept = 0, linetype=2)+
  coord_flip()+
  xlab('Variable')
p
dev.off()


## MCMCGLMM approach?? # https://ourcodingclub.github.io/tutorials/mcmcglmm/
# This paper suggested bayesian for data with outliers
# https://onlinelibrary.wiley.com/doi/10.1002/bimj.201800071
mcmc.mod1 <- MCMCglmm(LRR ~ FragmentationDesign, 
                      random = ~ID, 
                      family = "gaussian", data=frag_dat, 
                      verbose=FALSE, nitt = 100000, thin = 40) 


hist(mcmc(mcmc.mod1$VCV)[,"ID"]) # Random effect seems ok

plot(mcmc.mod1$Sol) # intercepts looks good

plot(mcmc.mod1$VCV) # random effect

a <- 1000
prior <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a), # G for every random effect
                        G1 = list(V = diag(1), fix = 1))) # for the idh() #  Studies with higher standard error have been given lower statistical weight.

mcmc.mod2 <- MCMCglmm(LRR ~ FragmentationDesign, 
                      random = ~ID + idh(LRR_var):units, 
                      family = "gaussian", data=frag_dat, 
                      prior = prior,
                      verbose=FALSE, nitt = 100000, thin = 40) 
  
  
plot(mcmc.mod2$VCV)


# Model checking
xsim <- simulate(mcmc.mod2) # reruns 100 new models, based around the same variance/covariance structures but with simulated data.

plot(frag_dat$LRR, I(1/frag_dat$LRR_var))
points(xsim, I(1/frag_dat$LRR_var), col = "red") # here you can plot the data from both your simulated and real datasets and compare them
# Not bad?
summary(mcmc.mod2)
traceplot(mcmc(mcmc.mod2$Sol))

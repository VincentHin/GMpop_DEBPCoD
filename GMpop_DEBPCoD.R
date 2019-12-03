## Specify your working directory
getwd()
setwd("~/research/git/GMpop_DEBPCoD/")

## Source functions into separate environment
detach(mmf)
mmf <- new.env()
source('GMpop_DEBPCoD_funcs.R', local = mmf)
attach(mmf); rm(mmf)

## Import timeseries output (GM_pop):
GMpop_default <- import.pop.output(folder = "~/research/MarineMammals/GM_output/competition/timeseries/July2019_commit_e7688e0/default/", baseflnm = "GM_pop", bifurcation = FALSE, confine.R0 = FALSE)
lapply(GMpop_default, tail)

source("Figure1.R")

## Import bifurcation output (DIST*_RESAMP*)
popbif.DISTs_RESAMP0 <-  import.pop.output(baseflnm = 'DISTs_RESAMP0.0', bifurcation = TRUE)
popbif.DISTs_RESAMP0.25 <-  import.pop.output(baseflnm = 'DISTs_RESAMP0.25', bifurcation = TRUE)
popbif.DISTw_RESAMP0.25 <-  import.pop.output(baseflnm = 'DISTw_RESAMP0.25', bifurcation = TRUE)

source("Figure2.R")


## --- popbif R0 output : ->
popbif.DIST_RESAMP0.R0 <- days.to.years(popbif.DISTs_RESAMP0$R0)
popbif.DISTs_RESAMP0.25.R0 <- days.to.years(popbif.DISTs_RESAMP0.25$R0)
popbif.DISTw_RESAMP0.25.R0 <- days.to.years(popbif.DISTw_RESAMP0.25$R0)

## Create popbif.DIST_RESAMP.R0 data.frame
popbif.DIST_RESAMP0.R0$dist_resamp <- 'Summer_0'
popbif.DISTs_RESAMP0.25.R0$dist_resamp <- 'Summer_0.25'
popbif.DISTw_RESAMP0.25.R0$dist_resamp <- 'Winter_0.25'
popbif.DIST_RESAMP.R0 <- rbind(popbif.DIST_RESAMP0.R0, popbif.DISTs_RESAMP0.25.R0, popbif.DISTw_RESAMP0.25.R0)
rm(popbif.DIST_RESAMP0.R0, popbif.DISTs_RESAMP0.25.R0, popbif.DISTw_RESAMP0.25.R0)

## Calculate additional columns
popbif.DIST_RESAMP.R0 <- size.from.age(popbif.DIST_RESAMP.R0)
popbif.DIST_RESAMP.R0 <- set.ageclass(popbif.DIST_RESAMP.R0)
popbif.DIST_RESAMP.R0 <- set.age3class(popbif.DIST_RESAMP.R0)
popbif.DIST_RESAMP.R0$yearclass <- floor(popbif.DIST_RESAMP.R0$age_yrs)
popbif.DIST_RESAMP.R0$ageclass[which(popbif.DIST_RESAMP.R0$ageclass == 'scenescent')] <- 'mature'
popbif.DIST_RESAMP.R0 <- calc.expected.age(popbif.DIST_RESAMP.R0) # Slooooowww lookup table.....

## Assign a cause of death to each individual 
## (1 = No starvation, 2 = Reduced life exp. but not starving, 3 = Reduced life exp. & starving)
popbif.DIST_RESAMP.R0$death_cause <- NA
popbif.DIST_RESAMP.R0$death_cause[popbif.DIST_RESAMP.R0$days_reduced <= 5] <- 1
popbif.DIST_RESAMP.R0$death_cause[popbif.DIST_RESAMP.R0$days_reduced > 5 & popbif.DIST_RESAMP.R0$fatratio > 0.15] <- 2
popbif.DIST_RESAMP.R0$death_cause[popbif.DIST_RESAMP.R0$days_reduced > 5 & popbif.DIST_RESAMP.R0$fatratio < 0.15] <- 3
popbif.DIST_RESAMP.R0$death_cause <- as.factor(popbif.DIST_RESAMP.R0$death_cause)

## Remove Winter_0.25 at 34 days disturbance, which goes extinct
idx <- which(popbif.DIST_RESAMP.R0$dist_resamp == "Winter_0.25" & popbif.DIST_RESAMP.R0$par1 == 34)
popbif.DIST_RESAMP.R0 <- popbif.DIST_RESAMP.R0[-idx,]

## Split into a list based on disturbance period and dist_resamp
popbif.DIST_RESAMP.R0.lst <- dlply(popbif.DIST_RESAMP.R0, .variables = c('dist_resamp','par1'))

## Calculate cause of death per age class
popbif.DIST_RESAMP.R0.deathcause <- ldply(popbif.DIST_RESAMP.R0.lst, function(x){
  
  ## Calculate some statistics in data.frame
  data.frame('dist_resamp' = head(x$dist_resamp,1),
    'par1' = head(x$par1,1),
    'value' = c(length(which(x$death_cause == 1)) / nrow(x),
      length(which(x$death_cause == 2)) / nrow(x),
      length(which(x$death_cause == 3)) / nrow(x),
      length(which(x$death_cause == 1 & x$age_yrs > 10)) / nrow(x),
      length(which(x$death_cause == 2 & x$age_yrs > 10)) / nrow(x),
      length(which(x$death_cause == 3 & x$age_yrs > 10)) / nrow(x),
      length(which(x$death_cause == 1 & x$age < 1223)) / nrow(x),
      length(which(x$death_cause == 2 & x$age < 1223)) / nrow(x),
      length(which(x$death_cause == 3 & x$age < 1223)) / nrow(x),
      length(which(x$death_cause == 1 & x$age > 1223 & x$age_yrs < 10)) / nrow(x),
      length(which(x$death_cause == 2 & x$age > 1223 & x$age_yrs < 10)) / nrow(x),
      length(which(x$death_cause == 3 & x$age > 1223 & x$age_yrs < 10)) / nrow(x)),
    'deathcause' = rep(factor(c(1:3), levels = c(3,2,1)),4),
    'ageclass' = rep(factor(c('all', 'ageplus10', 'calf', 'agemin10'), levels = c('all','calf','agemin10','ageplus10')), each = 3)
  )
}
)
popbif.DIST_RESAMP.R0.deathcause$.id <- NULL
head(popbif.DIST_RESAMP.R0.deathcause)

source("Figure3.R")
source("Figure5.R")

## Calculate bootstrapped means +/- 95% CL
popbif.DIST_RESAMP.R0.stats <- ldply(popbif.DIST_RESAMP.R0.lst, function(x){
  
  mean.ageatdeath <- smean.cl.boot(x$age_yrs, na.rm = TRUE)
  mean.femcalves.oldfems <- smean.cl.boot(x$IDtotalfemcalves[x$age_yrs > 10], na.rm = TRUE)
  mean.femcalves <- smean.cl.boot(x$IDtotalfemcalves, na.rm = TRUE)
  mean.agefirstrepro <- smean.cl.boot(x$IDagefirstrepro_zero, na.rm = TRUE)
  mean.agefirstweaning <- smean.cl.boot(x$IDagefirstweaning_zero, na.rm = TRUE)
  
  ## Calculate some statistics in data.frame
  data.frame('dist_resamp' = head(x$dist_resamp,1),
    'par1' = head(x$par1,1),
    'no.fems' = nrow(x),
    'mean.ageatdeath' = as.numeric(mean.ageatdeath[1]),
    'min.ageatdeath' = as.numeric(mean.ageatdeath[2]),
    'max.ageatdeath' = as.numeric(mean.ageatdeath[3]),
    'mean.totalfemcalves' = as.numeric(mean.femcalves[1]),
    'min.totalfemcalves' = as.numeric(mean.femcalves[2]),
    'max.totalfemcalves' = as.numeric(mean.femcalves[3]),
    'mean.totalfemcalves_oldfems' = as.numeric(mean.femcalves.oldfems[1]),
    'min.totalfemcalves_oldfems' = as.numeric(mean.femcalves.oldfems[2]),
    'max.totalfemcalves_oldfems' = as.numeric(mean.femcalves.oldfems[3]),
    'mean.agefirstrepro' = as.numeric(mean.agefirstrepro[1]),
    'min.agefirstrepro' = as.numeric(mean.agefirstrepro[2]),
    'max.agefirstrepro' = as.numeric(mean.agefirstrepro[3]),
    'mean.agefirstweaning' = as.numeric(mean.agefirstweaning[1]),
    'min.agefirstweaning' = as.numeric(mean.agefirstweaning[2]),
    'max.agefirstweaning' = as.numeric(mean.agefirstweaning[3]))
})
##
popbif.DIST_RESAMP.R0.stats$.id <- NULL
head(popbif.DIST_RESAMP.R0.stats)

source("Figure6.R")
##
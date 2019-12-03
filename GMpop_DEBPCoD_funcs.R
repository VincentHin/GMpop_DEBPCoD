if (!exists("par.defaults")) par.defaults <- par()
par.defaults$cin <- par.defaults$cra <- par.defaults$csi <- par.defaults$cxy <- par.defaults$din <-  par.defaults$page <- NULL

# Load some packages
library('stringr'); library("RColorBrewer"); library('PSPManalysis'); library('readr'); library('tibble'); library('gridExtra')
library('reshape2');  library('ggplot2'); library('plyr'); library('tidyr'); library('Hmisc'); library('ggExtra')
library('scales')

## : -> GGOPTS ##
ggopts <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.key.height = unit(9, "pt"),
    legend.position = 'bottom',
    legend.box = "horizontal",
    legend.margin = margin(0,0,0,0), 
    legend.box.margin = margin(0,0,0,0),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(size = 9, hjust = 0.5))

ggopts_lollipop <- theme(plot.margin = unit(c(0.25,0.5,0,0.5), 'cm'),
                         strip.background = element_blank(),
                         strip.text.x = element_text(hjust = 0),
                         panel.spacing.y = unit(0.1, 'in'),
                         panel.border = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         panel.grid.major.y = element_blank(),
                         axis.ticks = element_blank())

## : -> PL.OPTS ##
pl.opts <- list(
  lwidth = 3,
  cex.text = 0.65,
  cex.form = 0.65,
  cex.labs = 0.75,
  cex.axs = 0.65)

names.pop.R0 <- c('par1', 'par2', 'time', 'resource', 
  'IDindividualID', 'age', 'reserves', 'survival', 'life_expectancy', 
  'IDtotalfemcalves', 'IDtotalcalves', 'IDtotalweaned',
  'IDagefirstreceptive', 'IDagefirstrepro', 'IDagefirstweaning', 'IDstarvdays')
names.pop.state <- c('survival', 'age', 'reserves', 'totalmort',
  'IDmale', 'IDfemale', 'IDwaiting', 'IDpregnant', 'IDlactating', 'IDindividualID', 'IDmotherID', 'IDmotherIndex', 
  'IDminloglifespan', 'IDpregwait', 'IDinseminated', 'IDcalved', 'IDweaned', 
  'IDtotalfemcalves', 'IDtotalcalves', 'IDtotalweaned', 
  'IDagefirstreceptive', 'IDagefirstbirth', 'IDagefirstweaning',
  'IDlength', 'IDbones', 'IDweight', 'IDweightM', 'IDfatratio', 
  'IDingest', 'IDringest', 'IDmingest', 'IDmaint', 'IDgrowth', 'IDpregcosts', 'IDlactcosts', 'IDnetenergy', 
  'IDmortality', 'IDbackground', 'IDstarvation', 'IDstarvdays', 'IDfeedinglevel')
names.pop <- c('time', 'time_yrs', 'total', 'male', 'female', 'resting', 'pregnant', 'lactating', 'preglact', 'waiting', 'waitlact', 
  'yearling', 'male_yearling', 'female_yearling', 'suckling', 'male_suckling', 'female_suckling', 'juvenile', 'male_juvenile', 'female_juvenile', 
  'maturegrowing', 'male_maturegrowing', 'female_maturegrowing', 'mature', 'male_mature', 'female_mature', 'scenescent', 'male_scenescent', 'female_scenescent',
  'meanR0', 'sdR0', 'minR0', 'maxR0', 'meanTotalCalves', 'sdTotalCalves', 'minTotalCalves', 'maxTotalCalves', 'meanTotalWeaned', 'sdTotalWeaned', 'minTotalWeaned', 'maxTotalWeaned',
  'meanAgeFirstReceptive', 'sdAgeFirstReceptive', 'minAgeFirstReceptive', 'maxAgeFirstReceptive', 'meanAgeFirstBirth', 'sdAgeFirstBirth', 'minAgeFirstBirth', 'maxAgeFirstBirth',
  'meanAgeFirstWeaning', 'sdAgeFirstWeaning', 'minAgeFirstWeaning', 'maxAgeFirstWeaning', 
  'meanIBI', 'sdIBI', 'minIBI', 'maxIBI',
  'meanIWI', 'sdIWI', 'minIWI', 'maxIWI',
  'meanStarvDays', 'sdStarvDays', 'minStarvDays', 'maxStarvDays', 
  'meanReproFemales', 'sdReproFemales', 'meanStarvingFemales', 'sdStarvingFemales', 
  'Observations', 'ObsAgeFirstReceptive', 'ObsAgeFirstBirth', 'ObsAgeFirstWeaning', 'ObsStarvDays', 'ObsIBI', 'ObsIWI', 'resource')
names.deaths <- c('time_yrs', 'male_yearling', 'male_suckling', 'male_juvenile', 'male_growing', 'male_mature', 'male_scenescent',
  'female_yearling', 'female_suckling', 'female_juvenile', 'female_growing', 'female_mature', 'female_scenescent',
  'resting_yearling', 'resting_suckling', 'resting_juvenile', 'resting_growing', 'resting_mature', 'resting_scenescent',
  'pregnant_yearling', 'pregnant_suckling', 'pregnant_juvenile', 'pregnant_growing', 'pregnant_mature', 'pregnant_scenescent',
  'lactating_yearling', 'lactating_suckling', 'lactating_juvenile', 'lactating_growing', 'lactating_mature', 'lactating_scenescent',
  'preglact_yearling', 'preglact_suckling', 'preglact_juvenile', 'preglact_growing', 'preglact_mature', 'preglact_scenescent',
  'waiting_yearling', 'waiting_suckling', 'waiting_juvenile', 'waiting_growing', 'waiting_mature', 'waiting_scenescent',
  'waitlact_yearling', 'waitlact_suckling', 'waitlact_juvenile', 'waitlact_growing', 'waitlact_mature', 'waitlact_scenescent')
names.pop_bif <- c(names.pop, 'bifpar')

## Define colours
GMcolours <- 
  list('total' = adjustcolor('black', alpha.f=0.75), 
    'male' = adjustcolor("#4292C6", alpha.f=0.75),
    'female' = adjustcolor("#D95F02", alpha.f=0.75),
    'resource'=adjustcolor('forestgreen', alpha.f = 0.75),
    'calf' = "#56B4E9",
    'calf_trans' = adjustcolor("#56B4E9", alpha.f = 0.75),
    'resting' = "#009E73",
    'resting_trans' = adjustcolor("#009E73", alpha.f = 0.75),
    'waiting' = "#0072B2",
    'waiting_trans' = adjustcolor("#0072B2", alpha.f = 0.75),
    'pregnant' = 'darkgoldenrod2',
    'pregnant_trans' = adjustcolor('darkgoldenrod2', alpha.f = 0.75),
    'lactating' = "#E7298A",
    'lactating_trans' = adjustcolor("#E7298A", alpha.f = 0.75),
    'waitlact' = '#7570B3',
    'waitlact_trans' = adjustcolor('#7570B3', alpha.f = 0.75),
    'preglact' = "#A6761D",
    'preglact_trans' = adjustcolor("#A6761D", alpha.f = 0.75),
    'netenergy' = '#F8766D',
    'netenergy_trans' = adjustcolor('#F8766D', alpha.f = 0.75),
    'class1' = "#00AFBB",
    'class2' = "#E7B800",
    'class3' = "#FC4E07"
)
##

## Read population output from bifurcation or timeseries run. Searches within 'folder' for folder or files that match 'name'.
## Creates one list per name run. List elements are: .out, _deaths.out, (.minmax.out, .avg.out, .gavg.out, var.out_), rep, esf, R0
import.pop.output <- function(baseflnm = NULL, folder = getwd(), bifurcation = TRUE, confine.R0 = TRUE){
  if(is.null(baseflnm)) return('Enter a name or a component of a name')
  ## strip .cvf part of the name (if any)
  baseflnm <- gsub(pattern = '.cvf$', replacement = '', x = baseflnm)
  
  cat(paste('Check for files/folders named', baseflnm, 'in directory', folder, '\n'))
  if(bifurcation == TRUE) nms <- names.pop_bif else nms <- names.pop
  
  outfl <- file.path(folder, paste0(baseflnm, ".out"))
  minmaxfl <- file.path(folder, paste0(baseflnm, ".minmax.out"))
  gavgfl <- file.path(folder, paste0(baseflnm, ".gavg.out"))
  avgfl <- file.path(folder, paste0(baseflnm, ".avg.out"))
  varfl <- file.path(folder, paste0(baseflnm, ".var.out"))
  deathfl <- file.path(folder, paste0(baseflnm, "_deaths.out"))
  R0fl <- file.path(folder, paste0(baseflnm, "_R0.dat"))
  repfl <- file.path(folder, paste0(baseflnm, ".rep"))
  esffl <- file.path(folder, paste0(baseflnm, ".esf"))
  
  ## Import .out file
  if(length(outfl) != 0){
    cat(paste('\n Importing file:', outfl, '\n'))
    out.lst <- list(read.table(file = outfl, header = FALSE, sep = '\t', col.names = nms))
    names(out.lst) <- 'out'
  } else {
    out.lst <- NULL
  }
  
  ## Import deathfl
  if(length(deathfl) != 0){
    cat(paste('\n Importing file:', deathfl, '\n'))
    death.lst <- list(read.table(deathfl, col.names = names.deaths))
    names(death.lst) <- 'deaths.out'
  } else {
    death.lst <- NULL
  }
  
  ## Import R0 file
  if(length(R0fl) != 0){
    cat(paste('\n Importing file:', R0fl, '\n'))
    R0.lst <- lapply(R0fl, read.table, col.names =  names.pop.R0)
    names(R0.lst) <- 'R0'
  } else {
    R0.lst <- NULL
  }
  
  ## Import esf file
  if(length(esffl) != 0){
    cat(paste('\n Importing file:', esffl, '\n'))
    esf.lst <- lapply(esffl, read.table, col.names = names.pop.state, skip = 2)
    names(esf.lst) <- 'esf'
  } else {
    esf.lst <- NULL
  }
  
  ## Import rep file
  if(length(repfl) != 0){
    cat(paste('\n Importing file:', repfl, '\n'))
    rep.lst <- lapply(repfl, read_lines)
    names(rep.lst) <- 'rep'
  } else {
    rep.lst <- NULL
  }
  
  lst <- c(out.lst, death.lst, R0.lst, esf.lst, rep.lst)
  
  if(bifurcation == TRUE){
    
    ## Import .minmax.out file
    if(length(minmaxfl) != 0){
      cat(paste('\n Importing file:', minmaxfl, '\n'))
      minmax.lst <- list(read_delim(file = minmaxfl, delim = '\t', col_names = c(nms, 'period1', 'period2')))
      names(minmax.lst) <- 'minmax.out'
      minmax.lst$minmax.out$stat <- c('min','max')
    } else {
      minmax.lst <- NULL
    }
    
    ## Import .gavg.out file
    if(length(gavgfl) != 0){
      cat(paste('\n Importing file:', gavgfl, '\n'))
      gavg.lst <- list(read_delim(file = gavgfl, delim = '\t', col_names = c(nms, 'period1', 'period2')))
      names(gavg.lst) <- 'gavg.out'
      gavg.lst$gavg.out$stat <- 'gavg'
    } else {
      gavg.lst <- NULL
    }
    
    ## Import .avg.out file
    if(length(avgfl) != 0){
      cat(paste('\n Importing file:', avgfl, '\n'))
      avg.lst <- list(read_delim(file = avgfl, delim = '\t', col_names = c(nms, 'period1', 'period2')))
      names(avg.lst) <- 'avg.out'
      avg.lst$avg.out$stat <- 'avg'
    } else {
      avg.lst <- NULL
    }
    
    ## Import .var.out file
    if(length(varfl) != 0){
      cat(paste('\n Importing file:', varfl, '\n'))
      var.lst <- list(read_delim(file = varfl, delim = '\t', col_names = c(nms, 'period1', 'period2')))
      names(var.lst) <- 'var.out'
      var.lst$var.out$stat <- 'var'
    } else {
      var.lst <- NULL
    }
    
    lst <- c(lst, minmax.lst, gavg.lst, avg.lst, var.lst)
  }
  
  if(length(R0fl) != 0 & confine.R0 == TRUE) lst <- confine.R0.output(input.lst = lst)
  
  return(lst)
}
##

## Takes a list of pop bif output and confines R0 output to output bifperiod. Is used within import.pop.output()
confine.R0.output <- function(input.lst = NULL){
  
  ## Get unique bifurcation parameters
  bifpars <- unique(input.lst$out$bifpar)
  
  ## Get list of time ranges from .out data
  idx <-
    llply(as.list(bifpars), function(x) { 
      idx <- range( input.lst$out$time[input.lst$out$bifpar == x] )
      c(idx, x)
    })
  
  ## Add bifpar column to R0 data
  input.lst$R0$bifpar <- NA
  
  input.lst$R0 <- 
    ldply(idx, function(x) {
      set <- subset(input.lst$R0, time >= x[1] & time <= x[2])
      if(nrow(set) != 0) set$bifpar <- x[3]
      return(set)
    })
  
  cat('\n Confined R0 file to time range of bifurcation output generation.\n\n')
  
  return(input.lst)
}
##

## Transposes columns IDagefirstrepro, IDagefirstweaning, IDagefirstreceptive, IBI and IWI from days to years
## Calculate some additional columns that exclude zero. Calculates PercWeaned from IDearlydeaths.
days.to.years <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify input dataframe')
  
  input.df$IDagefirstrepro <- input.df$IDagefirstrepro / 365
  input.df$IDagefirstrepro_zero <- input.df$IDagefirstrepro
  input.df$IDagefirstrepro_zero[input.df$IDagefirstrepro_zero == 0] <- NA
  input.df$Reproductive <- as.numeric(input.df$IDagefirstrepro > 0)
  if(length(input.df$age) > 0) input.df$age_yrs <- input.df$age / 365
  if(!is.null(input.df$IDearlydeaths)) input.df$PercWeaned <- (1 - input.df$IDearlydeaths / (2 * input.df$IDreprofemale + input.df$IDearlydeaths))
  
  if(length(input.df$IDagefirstweaning) > 0){
    input.df$IDagefirstweaning <- input.df$IDagefirstweaning / 365
    input.df$IDagefirstweaning_zero <- input.df$IDagefirstweaning
    input.df$IDagefirstweaning_zero[input.df$IDagefirstweaning_zero == 0] <- NA
  }
  
  if(length(input.df$IDagefirstreceptive) > 0){
    input.df$IDagefirstreceptive <- input.df$IDagefirstreceptive / 365
    input.df$IDagefirstreceptive_zero <- input.df$IDagefirstreceptive
    input.df$IDagefirstreceptive_zero[input.df$IDagefirstreceptive_zero == 0] <- NA
  }
  
  idx.zero <- match(c('meanIBI', 'sdIBI', 'minIBI', 'maxIBI', 'meanIWI', 'sdIWI', 'minIWI', 'maxIWI'), names(input.df))
  if(!is.na(idx.zero)[1]){
    input.df[,idx.zero] <- input.df[,idx.zero] / 365
    input.df$meanIBI_zero <- input.df$meanIBI
    input.df$minIBI_zero <- input.df$minIBI
    input.df$maxIBI_zero <- input.df$maxIBI
    input.df$meanIWI_zero <- input.df$meanIWI
    input.df$minIWI_zero <- input.df$minIWI
    input.df$maxIWI_zero <- input.df$maxIWI
    
    input.df$meanIBI_zero[input.df$meanIBI_zero == 0] <- NA
    input.df$minIBI_zero[input.df$minIBI_zero == 0] <- NA
    input.df$maxIBI_zero[input.df$maxIBI_zero == 0] <- NA
    input.df$meanIWI_zero[input.df$meanIWI_zero == 0] <- NA
    input.df$minIWI_zero[input.df$minIWI_zero == 0] <- NA
    input.df$maxIWI_zero[input.df$maxIWI_zero == 0] <- NA
  }
  
  return(input.df)
}
##

## Calculate the structure length / mass from age
size.from.age <- function(input.df = NULL, pars = list(linf = 450, lb = 177, k = 0.00045, omega1 = 8.5E-5, omega2 = 2.6)){
  if(is.null(input.df)) return('Please specify input dataframe')
  if(is.null(pars)) return('Please specfic parameters a list names pars')
  
  lengthage = expression(linf - (linf - lb)*exp(-k*a))
  structage = expression(omega1 * (linf - (linf - lb)*exp(-k*a))^omega2)
  
  input.df$length <- eval(lengthage, c(pars, list(a = input.df$age)))
  input.df$structure <- eval(structage, c(pars, list(a = input.df$age)))
  
  if(!is.null(input.df$reserves)){
    input.df$weight <- input.df$structure + input.df$reserves
    input.df$fatratio <- input.df$reserves / input.df$weight
  }
  return(input.df)
}
##

## Set ageclass (5 ageclasses)
set.ageclass <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify an input dataframe with input.df')
  
  # Female is initiated as juveniles (age at weaning)
  input.df$ageclass <- 'scenescent'
  input.df$ageclass[input.df$age_yrs < 25] <- 'mature'
  input.df$ageclass[input.df$age_yrs < 15] <- 'maturegrowing'
  input.df$ageclass[input.df$age_yrs < 8] <- 'juvenile'
  input.df$ageclass[input.df$age_yrs < (1223/365)] <- 'calve'
  
  return(input.df)
}
##

## Set ageclass (3 ageclasses)
set.age3class <- function(input.df = NULL){
  if(is.null(input.df)) return('Please specify an input dataframe with input.df')
  
  # Female is initiated as juveniles (age at weaning)
  input.df$age3class <- '8plus'
  input.df$age3class[input.df$age_yrs < 8] <- '8min'
  input.df$age3class[input.df$age_yrs < (1223/365)] <- 'calve'
  
  return(input.df)
}
##

## Calculate the expected age at death from a life_expectancy column of input.df. Requires a lookup table and is slooooowwww....
calc.expected.age <- function(input.df = NULL, pars = list(alpha1 = 4.01e-4, alpha2 = 6.04e-6, beta1 = 5.82e-4, beta2 = 3.01E-4)){
  if(is.null(input.df)) return('Please specify a data.frame')
  if(is.null(input.df$life_expectancy)) input.df$life_expectancy <- exp(-input.df$IDminloglifespan)
  
  surv <- expression(exp((alpha1*exp(-beta1*a)*beta2-alpha2*exp(beta2*a)*beta1+alpha2*beta1-beta2*alpha1)/(beta1*beta2)))
  ages <- seq(0,60*365,by=0.5)
  survage <- eval(surv, c(pars, list(a = ages)))
  
  ## Using lapply is slighly faster then for loop
  life_exp <- as.list(input.df$life_expectancy)
  age_exp <- lapply(life_exp, function(x) { ages[which.min(abs(survage - x))] })
  
  ## For loop approach (slower)
  # age_exp <- vector(mode = "list", nrow(input.df))
  # life_exp <- as.list(input.df$life_expectancy)
  # for(i in 1:length(life_exp)){
  #   idx <- which.min(abs(survage - life_exp[[i]]))
  #   age_exp[[i]] <- ages[idx]
  # }
  
  input.df$age_expectancy <- unlist(age_exp)
  input.df$days_reduced <- (input.df$age_expectancy - input.df$age)
  
  ##
  return(input.df)
}
##

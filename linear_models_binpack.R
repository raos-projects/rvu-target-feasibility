#!/usr/bin/env Rscript

##################################################
## Project: RVU Target Feasibility Study
## Script purpose: Model outcome metrics across simulations
##                 as functions of simulation parameters
## Date: 12/25/2020
## Author: Saieesh Rao
##################################################

####  LIBRARIES   ####
library(sandwich)
library(lmtest)
library(dplyr)
library(magrittr)
library(purrr)
library(future)
library(furrr)
library(data.table)

####  CONSTANTS  ####

on.hpc = TRUE

if(on.hpc){
  root_results = "/scratch/t.sur.rsaieesh/results-binpack/"
} else {
  root_results = "X:/results-binpack/"
}

options(verbose = TRUE)

##    RAND SEED   ##

#The base seed is "PHEMISTER" converted into phone numbers.
#Dallas B. Phemister was a former professor and chairman of
#the Department of Surgery at the University of Chicago. He
#was a former president of the American Surgical Association
#and the American College of Surgeons, as well as a member of
#the editorial board of Annals of Surgery.

SEED <- 743647837

##  CASEMIX COMPLEXITY  ##

#This is a metric defined by collecting cases with identical
#principal procedures and splitting them into 3 quantiles
#based on operating time. It's a crude measure that satisfies
#an intuitive assumption that more complex cases take more
#time but also generate more RVUs on average. 

Casemix.LOW = 1 #"Fast"
Casemix.MEDIUM = 2 #"Medium"
Casemix.HIGH = 3 #"Slow"
Casemix.ALL = 0 #"All"

casemixlist <- c(Casemix.ALL,Casemix.LOW,Casemix.MEDIUM,Casemix.HIGH)

##  SCHEDULE  ##

#no longer using scheduling patterns; instead, cluster of
#cases to be grouped at one time (5,10,20,30,99999) by the
#bin-packing algorithm. Sched.ANNUAL is an arbitrarily large
#constant meant to encompass all cases in an annual caselog.

Sched.ANNUAL = 99999

##  SPECIALTIES  ##

cardiac <- "Cardiac Surgery"
general <- "General Surgery"
gyne <- "Gynecology"
ir <- "Interventional Radiologist"
neuro <- "Neurosurgery"
ortho <- "Orthopedics"
ent <- "Otolaryngology (ENT)"
plastics <- "Plastics"
thoracic <- "Thoracic"
urology <- "Urology"
vascular <- "Vascular"

# list of specialties with ir removed
speclist <- c(cardiac, general, gyne, neuro, ortho,
              ent, plastics, thoracic, urology, vascular)

####   SUBR: IMPORT SINGLE SIMULATION DATA   ####

# defines a new import function, which is used to read in a single csv
# file from the scratch folder. The file specified contains block-level
# data on all 100 or 1000 trials performed under a single set of 
# simulation parameters. The filename is generated algorithmically from
# the simulation parameters passed as arguments to the import function.

import.sim <- function(.specialty, .rvu_goal, .block_hours,
                       .cluster_size, .casemix, .turnover_time) {
  
  #systematically generate filename and path based on trial settings
  
  if(.cluster_size == Sched.ANNUAL){
    middle_path = 'binpack-global/bp_'
  }
  else{
    middle_path = 'binpack/bp_'
  }
  
  if(on.hpc){
    path <- paste0("/scratch/t.sur.rsaieesh/",middle_path,.specialty,
                   "_RVUs-",.rvu_goal,
                   "_blocksize-",.block_hours,
                   "_casemix-",.casemix,
                   "_turnover-",.turnover_time,
                   "-cluster-",.cluster_size,
                   ".csv")
  } else {
    path <- paste0("X:/",middle_path,.specialty,
                   "_RVUs-",.rvu_goal,
                   "_blocksize-",.block_hours,
                   "_casemix-",.casemix,
                   "_turnover-",.turnover_time,
                   "-cluster-",.cluster_size,
                   ".csv")
  }
  
  print(paste("importing",path))
  data.single <- fread(path)  #read in trial data
  
  #mark sim specs associated with each simulation
  data.single %<>% 
    mutate(specialty_spec = .specialty) %>% 
    mutate(cluster.size_spec = .cluster_size) %>%
    mutate(casemix_spec = .casemix) %>%
    mutate(rvu.target_spec = .rvu_goal) %>%
    mutate(block.size_spec = .block_hours*60) %>%
    mutate(turnover.time_spec = .turnover_time)
  
  #return data.table containing requested trial data as output from function
  return(data.single)            
}


####  SUBR: GENERIC CALCULATION OF MULTIVARIATE LINEAR MODELS   ####

# Create table containing coefficients from multivariate linear regression.
# Dependent response variable in regression is specified as '.var'. 
# Examples include hourly RVUs, net overtime, net turnover time, number
# of cases, etc. The independent variables are turnover time, block length,
# and optionally RVU target for dependent vars which are countable (not rates).
# Options to invert turnover time and block size allow for regressions against
# the multiplicative inverse (1/x) of those values rather than the raw values.
# This is especially helpful when certain metrics decrease with increasing
# values of turnover time or block size (such as n.blocks ~ block.size_spec).

coeff.table <- function(.data, .var, .regress_by_rvu=FALSE,
                        .invert_turnover=FALSE, .invert_blocksize=FALSE) {
  var <- substitute(.var)
  #convert mins to hours for regression
  .data %<>% mutate(block.size_spec = block.size_spec / 60)

  # algorithmically generate regression formulas from coeff.table arguments.
  # formulas returned as strings and reconstituted as formulas in appropriate environment
  # order of formula creation is RVU, Block Size, Turnover time (order is important later)
  synth.formula <- function(x.regress_by_rvu, x.invert_turnover, x.invert_blocksize) {
    
    #formula base
    fx <- c('eval(var) ~ ', 'eval(var) ~ ')
    
    if(is.na(x.regress_by_rvu)) {
      #only regress by RVU; can return formula now
      #invert arguments are irrelevant in this case
      fx_rvu <- c('rvu.target_spec','rvu.total')
      return(paste0(fx,fx_rvu))
      
    } else {
      if(x.regress_by_rvu) {
        #only inlcude rvu var if TRUE, otherwise default blank
        fx_rvu <- c('rvu.target_spec + ','rvu.total + ')
      } else {
        fx_rvu <- c()
      }
      #blocksize and turnover time are always included; variation is
      #whether to invert values or use raw values directly
      if(x.invert_blocksize) {
        fx_bls <- c('I(block.size_spec^-1) + ','I(block.size_spec^-1) + ')
      } else {
        fx_bls <- c('block.size_spec + ','block.size_spec + ')
      }
      if(x.invert_turnover) {
        fx_tot <- c('I(turnover.time_spec^-1)','I(turnover.time_spec^-1)')
      } else {
        fx_tot <- c('turnover.time_spec','turnover.time_spec')
      }
      
      #combine formula terms as string
      fx <- paste0(fx, paste0(fx_rvu,fx_bls,fx_tot))
      
      return(fx)
    }
  }
  
  formulas <- synth.formula(x.regress_by_rvu = .regress_by_rvu,
                            x.invert_turnover = .invert_turnover,
                            x.invert_blocksize = .invert_blocksize)
  formula_target <- formulas[1]
  formula_total <- formulas[2]
  
  # remove turnover time == 0 if inverting turnover.time_spec (can't divide by zero)
  if(.invert_turnover){
    .data %<>% filter(turnover.time_spec > 0)
    
    #exit if no data because turnover_time == 0 and was filtered for all data
    if(.data %>% dim() %>% .[1] == 0){
      return(NULL) #will result in invisible row during row binding
    }
  }
  
  # calculate ordinary least-squares (OLS) linear models with given formula
  model_target <- lm(as.formula(formula_target), .data)
  model_total <- lm(as.formula(formula_total), .data)
  
  # obtain model parameters and related values, depending on form of formula
  if(is.na(.regress_by_rvu) || .regress_by_rvu){
    
    #additional weighted least-squares (WLS) regression, weighting by RVU metric
    wls_target <- lm(as.formula(formula_target),
                     data = .data, weights = .data$rvu.target_spec^-1)
    wls_total <- lm(as.formula(formula_total),
                    data = .data, weights = .data$rvu.total^-1)
    
    # OLS model using RVU target
    coeffs_target <- model_target$coefficients
    summary_target <- summary(model_target)
    regvars_target <- summary_target$coefficients %>% row.names()
    #heteroscedastic consistent estimators
    hce_target <- model_target %>% coeftest(vcovHC=vcovHC(.,type='HC3'))
    hce_target_tv <- hce_target[,'t value']
    hce_target_pv <- hce_target[,'Pr(>|t|)']
    hce_target_se <- hce_target[,'Std. Error']
    AIC_target <- AIC(model_target)
    
    # WLS model using RVU target
    summary_wls_target <- summary(wls_target)
    coeffs_wls_target <- summary_wls_target %>% coefficients()
    wls_target_coeff <- coeffs_wls_target[,'Estimate']
    wls_target_tv <- coeffs_wls_target[,'t value']
    wls_target_pv <- coeffs_wls_target[,'Pr(>|t|)']
    wls_target_se <- coeffs_wls_target[,'Std. Error']
    AIC_wls_target <- AIC(wls_target)
    
    # OLS model using Total RVUs
    coeffs_total <- model_total$coefficients
    summary_total <- summary(model_total)
    regvars_total <- summary_total$coefficients %>% row.names()
    #heteroscedastic consistent estimators
    hce_total <- model_total %>% coeftest(vcov=vcovHC(.,type='HC3')) 
    hce_total_tv <- hce_total[,'t value']
    hce_total_pv <- hce_total[,'Pr(>|t|)']
    hce_total_se <- hce_total[,'Std. Error']
    AIC_total <- AIC(model_total)
    
    # WLS model using Total RVUs
    summary_wls_total <- summary(wls_total)
    coeffs_wls_total <- summary_wls_total %>% coefficients()
    wls_total_coeff <- coeffs_wls_total[,'Estimate']
    wls_total_tv <- coeffs_wls_total[,'t value']
    wls_total_pv <- coeffs_wls_total[,'Pr(>|t|)']
    wls_total_se <- coeffs_wls_total[,'Std. Error']
    AIC_wls_total <- AIC(wls_total)

  }
  else{
    #formula_target and formula_total are identical when .regress_by_rvu == FALSE
    #(there's no RVU term in the formula) so just pick one for model calculations.
    #additionally no need for weighting by RVU metric if not regressing by RVUs.
    
    coeffs <- model_target$coefficients
    summary <- summary(model_target)
    regvars <- summary$coefficients %>% row.names()
    #heteroscedastic consistent estimators
    hce <- model_target %>% coeftest(vcov=vcovHC(.,type='HC3'))
    #for each of below, [1='(Intercept)', 2=rvu.target_spec,
    #3=block.size_spec OR I(block.size_spec^-1),
    #4=turnover.time_spec OR I(turnover.time_spec^-1)]
    hce_tv <- hce[,'t value']
    hce_pv <- hce[,'Pr(>|t|)']
    hce_se <- hce[,'Std. Error']
    AIC_model <- AIC(model_target)
  }
  
  
  .specialty <- .data$specialty_spec %>% unique()         #should be one value
  .cluster_size <- .data$cluster.size_spec %>% unique()   #should be one value
  .casemix <- .data$casemix_spec %>% unique()             #should be one value
  .turnover_time <- .data$turnover.time_spec %>% unique()
  .block_hours <- .data$block.size_spec %>% unique()
  .rvu_target <- .data$rvu.target_spec %>% unique()
  #RVU total expected to be multiple values; mean just used for table
  .rvu_total = mean(.data$rvu.total)
  
  #if list OR null, need to convert to single value in output table
  if(length(.specialty) != 1) .specialty <- NA         
  if(length(.cluster_size) != 1) .cluster_size <- NA
  if(length(.casemix) != 1) .casemix <- NA
  if(length(.block_hours) != 1) .block_hours <- NA
  if(length(.turnover_time) != 1) .turnover_time <- NA 
  if(length(.rvu_target) != 1) .rvu_target <- NA
  
  #if regression contains RVU as independent var
  if(is.na(.regress_by_rvu) || .regress_by_rvu){
    
    output_target <- tribble(
      ~surgspec, ~cluster_size, ~casemix, ~turnover_time, ~block_hours,
      ~rvu_target, ~rvu_total, ~invert.turnover_time, ~invert.block_hours,
      ~coeff.rvu_target, ~coeff.rvu_total, ~coeff.block_hours,
      ~coeff.turnover_time, ~intercept,
      #heteroscedasticity consistent standard errors
      ~HCSE.rvu_target, ~HCSE.rvu_total, ~HCSE.block_hours,
      ~HCSE.turnover_time, ~HCSE.intercept,
      ~tvalue.rvu_target, ~tvalue.rvu_total, ~tvalue.block_hours,
      ~tvalue.turnover_time, ~tvalue.intercept,
      ~pvalue.rvu_target, ~pvalue.rvu_total, ~pvalue.block_hours,
      ~pvalue.turnover_time, ~pvalue.intercept,
      ~rsquared.mult, ~rsquared.adjust,
      ~stderr.regression, ~AIC_target, ~AIC_total,
      ~coeff.WLS.rvu_target, ~coeff.WLS.rvu_total, ~coeff.WLS.block_hours,
      ~coeff.WLS.turnover_time, ~intercept.WLS,
      ~SE.WLS.rvu_target, ~SE.WLS.rvu_total, ~SE.WLS.block_hours,
      ~SE.WLS.turnover_time, ~SE.WLS.intercept,
      ~tvalue.WLS.rvu_target, ~tvalue.WLS.rvu_total, ~tvalue.WLS.block_hours,
      ~tvalue.WLS.turnover_time, ~tvalue.WLS.intercept,
      ~pvalue.WLS.rvu_target, ~pvalue.WLS.rvu_total, ~pvalue.WLS.block_hours,
      ~pvalue.WLS.turnover_time, ~pvalue.WLS.intercept,
      ~rsquared.mult.WLS, ~rsquared.adjust.WLS,
      ~stderr.regression.WLS,~AIC_target.WLS,~AIC_total.WLS,
      
      .specialty, .cluster_size, .casemix, .turnover_time, .block_hours,
      .rvu_target, .rvu_total, .invert_turnover, .invert_blocksize,
      coeffs_target["rvu.target_spec"], NA, NA,
      NA, coeffs_target["(Intercept)"],
      NA,NA,NA,NA, hce_target_se["(Intercept)"],
      NA,NA,NA,NA, hce_target_tv["(Intercept)"],
      NA,NA,NA,NA, hce_target_pv["(Intercept)"],
      # NA,NA,NA,NA,summary_target$coefficients["(Intercept)",4],
      summary_target$r.squared, summary_target$adj.r.squared,
      summary_target$sigma, AIC_target, NA,
      wls_target_coeff["rvu.target_spec"], NA, NA,
      NA, wls_target_coeff["(Intercept)"],
      NA,NA,NA,
      NA, wls_target_se["(Intercept)"],
      NA,NA,NA,
      NA, wls_target_tv["(Intercept)"],
      NA,NA,NA,
      NA, wls_target_pv["(Intercept)"],
      summary_wls_target$r.squared, summary_wls_target$adj.r.squared,
      summary_wls_target$sigma, AIC_wls_target, NA
    )
    output_total <- tribble(
      ~surgspec, ~cluster_size, ~casemix, ~turnover_time, ~block_hours,
      ~rvu_target, ~rvu_total, ~invert.turnover_time, ~invert.block_hours,
      ~coeff.rvu_target, ~coeff.rvu_total, ~coeff.block_hours,
      ~coeff.turnover_time, ~intercept,
      #heteroscedasticity consistent errors
      ~HCSE.rvu_target, ~HCSE.rvu_total, ~HCSE.block_hours,
      ~HCSE.turnover_time, ~HCSE.intercept,
      ~tvalue.rvu_target, ~tvalue.rvu_total, ~tvalue.block_hours,
      ~tvalue.turnover_time, ~tvalue.intercept,
      ~pvalue.rvu_target, ~pvalue.rvu_total, ~pvalue.block_hours,
      ~pvalue.turnover_time, ~pvalue.intercept,
      ~rsquared.mult, ~rsquared.adjust,
      ~stderr.regression,~AIC_target,~AIC_total,
      ~coeff.WLS.rvu_target, ~coeff.WLS.rvu_total, ~coeff.WLS.block_hours,
      ~coeff.WLS.turnover_time, ~intercept.WLS,
      ~SE.WLS.rvu_target, ~SE.WLS.rvu_total, ~SE.WLS.block_hours,
      ~SE.WLS.turnover_time, ~SE.WLS.intercept,
      ~tvalue.WLS.rvu_target, ~tvalue.WLS.rvu_total, ~tvalue.WLS.block_hours,
      ~tvalue.WLS.turnover_time, ~tvalue.WLS.intercept,
      ~pvalue.WLS.rvu_target, ~pvalue.WLS.rvu_total, ~pvalue.WLS.block_hours,
      ~pvalue.WLS.turnover_time, ~pvalue.WLS.intercept,
      ~rsquared.mult.WLS, ~rsquared.adjust.WLS,
      ~stderr.regression.WLS,~AIC_target.WLS,~AIC_total.WLS,
      
      .specialty, .cluster_size, .casemix, .turnover_time, .block_hours,
      .rvu_target, .rvu_total, .invert_turnover, .invert_blocksize,
      NA, coeffs_total["rvu.total"], NA,
      NA, coeffs_total["(Intercept)"],
      NA,NA,NA,
      NA, hce_total_se["(Intercept)"],
      NA,NA,NA,
      NA, hce_total_tv["(Intercept)"],
      NA,NA,NA,
      NA, hce_total_pv["(Intercept)"],
      # NA,NA,NA,NA,summary_target$coefficients["(Intercept)",4],
      summary_total$r.squared, summary_total$adj.r.squared,
      summary_total$sigma, NA, AIC_total,
      NA, wls_total_coeff["rvu.total"], NA,
      NA, wls_total_coeff["(Intercept)"],
      NA,NA,NA,
      NA, wls_total_se["(Intercept)"],
      NA,NA,NA,
      NA, wls_total_tv["(Intercept)"],
      NA,NA,NA,
      NA, wls_total_pv["(Intercept)"],
      summary_wls_total$r.squared, summary_wls_total$adj.r.squared,
      summary_wls_total$sigma, NA, AIC_wls_total
    )
    
    #add rvu associated coeffs, since these regressions have rvus as independent variables
    output_target$HCSE.rvu_target <- hce_target_se["rvu.target_spec"]
    output_total$HCSE.rvu_total <- hce_total_se["rvu.total"]
    output_target$tvalue.rvu_target <- hce_target_tv["rvu.target_spec"]
    output_total$tvalue.rvu_total <- hce_total_tv["rvu.total"]
    output_target$pvalue.rvu_target <- hce_target_pv["rvu.target_spec"]
    output_total$pvalue.rvu_total <- hce_total_pv["rvu.total"]
    
    output_target$SE.WLS.rvu_target <- wls_target_se["rvu.target_spec"]
    output_total$SE.WLS.rvu_total <- wls_total_se["rvu.total"]
    output_target$tvalue.WLS.rvu_target <- wls_target_tv["rvu.target_spec"]
    output_total$tvalue.WLS.rvu_total <- wls_total_tv["rvu.total"]
    output_target$pvalue.WLS.rvu_target <- wls_target_pv["rvu.target_spec"]
    output_total$pvalue.WLS.rvu_total <- wls_total_pv["rvu.total"]
    
    #if regressed by things besides RVU
    if(!is.na(.regress_by_rvu)){
      
      #if block.size_spec was a regression variable:
      if(TRUE %in%
         (c('I(block.size_spec^-1)','block.size_spec') %in% regvars_target)){

        # missing name returns zero, so pick bigger number (non-zero) as index
        bs_ind <- max(regvars_target %>% {which(.=='block.size_spec')},
                      regvars_target %>% {which(.=='I(block.size_spec^-1)')})
        
        output_target$coeff.block_hours <- coeffs_target[bs_ind]
        output_total$coeff.block_hours <- coeffs_total[bs_ind]
        output_target$HCSE.block_hours <- hce_target_se[bs_ind]
        output_total$HCSE.block_hours <- hce_total_se[bs_ind]
        output_target$tvalue.block_hours <- hce_target_tv[bs_ind]
        output_total$tvalue.block_hours <- hce_total_tv[bs_ind]
        output_target$pvalue.block_hours <- hce_target_pv[bs_ind]
        output_total$pvalue.block_hours <- hce_total_pv[bs_ind]
        
        output_target$coeff.WLS.block_hours <- wls_target_coeff[bs_ind]
        output_total$coeff.WLS.block_hours <- wls_total_coeff[bs_ind]
        output_target$SE.WLS.block_hours <- wls_target_se[bs_ind]
        output_total$SE.WLS.block_hours <- wls_total_se[bs_ind]
        output_target$tvalue.WLS.block_hours <- wls_target_tv[bs_ind]
        output_total$tvalue.WLS.block_hours <- wls_total_tv[bs_ind]
        output_target$pvalue.WLS.block_hours <- wls_target_pv[bs_ind]
        output_total$pvalue.WLS.block_hours <- wls_total_pv[bs_ind]
        
      }
      #if turnover.time_spec was a regression variable:
      if(TRUE %in%
         (c('I(turnover.time_spec^-1)','turnover.time_spec') %in% regvars_target)){
        
        # missing name returns zero, so pick bigger number (non-zero) as index
        # regvars list for both rvu_target and rvu_total are identical for BS and TO
        to_ind <- max(regvars_target %>% {which(.=='turnover.time_spec')},
                      regvars_target %>% {which(.=='I(turnover.time_spec^-1)')})
        
        output_target$coeff.turnover_time <- coeffs_target[to_ind]
        output_total$coeff.turnover_time <- coeffs_total[to_ind]
        output_target$HCSE.turnover_time <- hce_target_se[to_ind]
        output_total$HCSE.turnover_time <- hce_total_se[to_ind]
        output_target$tvalue.turnover_time <- hce_target_tv[to_ind]
        output_total$tvalue.turnover_time <- hce_total_tv[to_ind]
        output_target$pvalue.turnover_time <- hce_target_pv[to_ind]
        output_total$pvalue.turnover_time <- hce_total_pv[to_ind]
        
        output_target$coeff.WLS.turnover_time <- wls_target_coeff[to_ind]
        output_total$coeff.WLS.turnover_time <- wls_total_coeff[to_ind]
        output_target$SE.WLS.turnover_time <- wls_target_se[to_ind]
        output_total$SE.WLS.turnover_time <- wls_total_se[to_ind]
        output_target$tvalue.WLS.turnover_time <- wls_target_tv[to_ind]
        output_total$tvalue.WLS.turnover_time <- wls_total_tv[to_ind]
        output_target$pvalue.WLS.turnover_time <- wls_target_pv[to_ind]
        output_total$pvalue.WLS.turnover_time <- wls_total_pv[to_ind]
      }
    }
    
    output <- bind_rows(list(output_target, output_total))
    
  }
  else{ #.regress_by_rvu == FALSE
    #regressions which do not regress by rvus, instead blocksize and turnover time only
    # warning(paste("coeffelse",.regress_by_rvu))
    output <- tribble(
      ~surgspec, ~cluster_size, ~casemix, ~turnover_time, ~block_hours,
      ~invert.turnover_time, ~invert.block_hours,
      ~coeff.block_hours, ~coeff.turnover_time, ~intercept,
      #heteroscedasticity consistent standard errors
      ~HCSE.block_hours, ~HCSE.turnover_time, ~HCSE.intercept, 
      ~tvalue.block_hours, ~tvalue.turnover_time, ~tvalue.intercept,
      ~pvalue.block_hours, ~pvalue.turnover_time, ~pvalue.intercept,
      ~rsquared.mult, ~rsquared.adjust, ~stderr.regression,~AIC,
      .specialty, .cluster_size, .casemix, .turnover_time, .block_hours,
      .invert_turnover, .invert_blocksize,
      coeffs["block.size_spec"], coeffs["turnover.time_spec"], coeffs["(Intercept)"],
      NA,NA, hce_se["(Intercept)"],
      NA,NA, hce_tv["(Intercept)"],
      NA,NA, hce_pv["(Intercept)"],
      # NA,NA,summary$coefficients["(Intercept)",4],
      summary$r.squared, summary$adj.r.squared, summary$sigma, AIC_model
    )
    
    #if block.size_spec was a regression variable:
    if(TRUE %in% (c('I(block.size_spec^-1)','block.size_spec') %in% regvars)){
      
      # missing name returns zero, so pick bigger number (non-zero) as index
      bs_ind <- max(regvars %>% {which(.=='block.size_spec')},
                    regvars %>% {which(.=='I(block.size_spec^-1)')})
      
      output$coeff.block_hours <- coeffs[bs_ind]
      output$HCSE.block_hours <- hce_se[bs_ind]
      output$tvalue.block_hours <- hce_tv[bs_ind]
      output$pvalue.block_hours <- hce_pv[bs_ind]
    }
    #if turnover.time_spec was a regression variable:
    if(TRUE %in% (c('I(turnover.time_spec^-1)','turnover.time_spec') %in% regvars)){
      
      # missing name returns zero, so pick bigger number (non-zero) as index
      to_ind <- max(regvars %>% {which(.=='turnover.time_spec')},
                    regvars %>% {which(.=='I(turnover.time_spec^-1)')})
      
      output$coeff.turnover_time <- coeffs[to_ind]
      output$HCSE.turnover_time <- hce_se[to_ind]
      output$tvalue.turnover_time <- hce_tv[to_ind]
      output$pvalue.turnover_time <- hce_pv[to_ind]
    }
  }
  
  #remove block hour and turnover time related columns if not regression variables
  if(is.na(.regress_by_rvu)){
    output %<>% select(-c(coeff.block_hours,coeff.turnover_time,
                          HCSE.block_hours, HCSE.turnover_time,
                          tvalue.block_hours, tvalue.turnover_time,
                          pvalue.block_hours, pvalue.turnover_time,
                          coeff.WLS.block_hours,coeff.WLS.turnover_time,
                          SE.WLS.block_hours, SE.WLS.turnover_time,
                          tvalue.WLS.block_hours, tvalue.WLS.turnover_time,
                          pvalue.WLS.block_hours, pvalue.WLS.turnover_time))
  }
  return(output)
}


####                              ####
####    SUBRS: ROW PREPARATION    ####
####                              ####


####  SUBR: BLOCK REQUIREMENTS   ####

# calculate the linear model for block requirements on data with specified simulation
# characteristics. If .turnover_time and/or .block_hours are NA, then a multivariate
# linear regression is computed on all sims within the variable's range.

# returns a row containing the regression coefficients for the specified simulation(s)

row.prepare.blockreqs <- function(.specialty, .casemix, .cluster_size,
                                  .turnover_time=seq(0,90,10), .block_hours=seq(4,12,2)) {
  
  # accept NA as argument specifying that variable can take values over entire range
  if(length(.turnover_time) == 1 && is.na(.turnover_time)) .turnover_time = seq(0,90,10)
  if(length(.block_hours) == 1 && is.na(.block_hours)) .block_hours = seq(4,12,2)
  
  # Pool data from sims with same specialty, schedule, speed, turnover time, and
  # block size (differ only in RVU target). RVU target is a sim parameter extrinsic
  # to daily block utilization, so blocks from sims that differ only in RVU target
  # are comparable and can be pooled together for stats analysis.
  
  # 1. specify sim data that can be pooled:
  specs <- expand.grid(.specialty = .specialty,
                       .cluster_size = .cluster_size,
                       .casemix = .casemix,
                       .rvu_goal = seq(3000,12000,1000),
                       .turnover_time = .turnover_time,
                       .block_hours = .block_hours,
                       stringsAsFactors = FALSE)
  
  # 2. import the data to be pooled. Cannot use future_pmap since 'row.prepare' is
  #    being called within future_pmap (and hence is already limited to one processor)
  data.raw <- pmap(specs, import.sim) %>%
    rbindlist() %>%
    group_by(trial, specialty_spec, cluster.size_spec,
             casemix_spec, rvu.target_spec, turnover.time_spec, block.size_spec) %>%
    summarize(n.blocks = max(block),
              rvu.total = sum(rvus)) %>%
    ungroup()

  # 3. calculate linear model
  block.tbl <- coeff.table(data.raw, n.blocks,
                           .regress_by_rvu = TRUE, .invert_blocksize = TRUE)
  
  return(block.tbl)
}

####  SUBR: TURNOVERS AND OVERTIME  ####

# Prepares one row summary of linear models for annual overtime (OT) and 
# turnover time (TO). Similar to above, just separated for clarity.

row.prepare.toot <- function(.specialty, .casemix, .cluster_size,
                             .turnover_time=seq(0,90,10), .block_hours=seq(4,12,2)) {
  
  # accept NA as argument specifying that variable can take values over entire range
  if(length(.turnover_time) == 1 && is.na(.turnover_time)) .turnover_time = seq(0,90,10)
  if(length(.block_hours) == 1 && is.na(.block_hours)) .block_hours = seq(4,12,2)
  
  # Pool data from sims with same specialty, schedule, speed, turnover time, and
  # block size (differ only in RVU target). RVU target is a sim parameter extrinsic
  # to daily block utilization, so blocks from sims that differ only in RVU target
  # are comparable and can be pooled together for stats analysis.
  
  # 1. specify sim data that can be pooled:
  specs <- expand.grid(.specialty = .specialty,
                       .casemix = .casemix,
                       .cluster_size = .cluster_size,
                       .rvu_goal = seq(3000,12000,1000),
                       .turnover_time = .turnover_time,
                       .block_hours = .block_hours,
                       stringsAsFactors = FALSE)
  
  # 2. import the data to be pooled. Cannot use future_pmap since 'row.prepare' is
  #    being called within future_pmap (and hence is already limited to one processor)
  data.raw <- pmap(specs, import.sim) %>%
    bind_rows() %>% 
    rename(time.actual = time.or) %>%
    rename(rvus.perblock = rvus) %>% 
    mutate(time.turnover = (n.cases-1)*turnover.time_spec)

  # 3. the below is collection of data by surgeon a.k.a annual schedule, rather
  #    than by block. Linear regression will be perfomed on the results of
  #    surgeon-level summaries.
  data.grouped <- data.raw %>% 
    #amount of time spent over block length
    mutate(ot.actual = pmax(time.actual - block.size_spec,0)) %>%
    #amount of OR time requested (reserved + overtime)
    mutate(OR.actual = pmax(time.actual, block.size_spec)) %>%     
    group_by(trial, specialty_spec, cluster.size_spec,
             casemix_spec, rvu.target_spec, turnover.time_spec, block.size_spec) %>%
    summarize(n.blocks = max(block),
              n.cases = sum(n.cases),
              time.actual.net = sum(time.actual),
              time.turnover.net = sum(time.turnover),
              ot.actual.net = sum(ot.actual),
              rvu.total = sum(rvus.perblock)) %>% 
    ungroup()
  
  # 4. fit five linear models:
  #       annual overtime ~ blocksize + turnover + rvu_target 
  #       annual turnover time ~ ...
  #       annual number of cases  ~ ...
  #       annual observed OR utilization ~ ...
  #       annual RVU total (>= RVU target) ~ ...

  #       annual overtime and turnover time must be pegged to RVU target
  #       rates (percent..., hourly rvus...) can use pooled data from all RVU targets
  
  ot.tbl <- data.grouped %>% coeff.table(ot.actual.net, .regress_by_rvu = TRUE)
  to.tbl <- data.grouped %>% coeff.table(time.turnover.net, .regress_by_rvu = TRUE)
  #n.cases ~ RVU_target only; independent of Block Size and Turnover Time
  ncases.tbl <- data.grouped %>% coeff.table(n.cases, .regress_by_rvu = NA) 
  #n.b. time.actual = time.or != optime
  timeactual.tbl <- data.grouped %>% coeff.table(time.actual.net, .regress_by_rvu = TRUE) 
  rvutotal.tbl <- data.grouped %>% coeff.table(rvu.total, .regress_by_rvu = TRUE)
  
  # 5. format and return output row with identifier columns at left and p values at right
  output <- list(ot.actual = ot.tbl,
                 turnover.net = to.tbl,
                 n.cases = ncases.tbl,
                 time.or = timeactual.tbl,
                 rvu.total = rvutotal.tbl)
  return(output)
}

####  SUBR: RATES (% OVERTIME AND HOURLY RVUS)  ####

# Prepares one row summary of linear models for annual
# percent overtime, hourly RVU production rates

row.prepare.pothrvu <- function(.specialty, .cluster_size, .casemix,
                                .turnover_time=seq(0,90,10), .block_hours=seq(4,12,2)) {
  
  # accept NA as argument specifying that variable can take values over entire range
  if(length(.turnover_time) == 1 && is.na(.turnover_time)) .turnover_time = seq(0,90,10)
  if(length(.block_hours) == 1 && is.na(.block_hours)) .block_hours = seq(4,12,2)
  
  # Pool data from sims with same specialty, schedule, speed, turnover time,
  # and block size (differ only in RVU target). RVU target is a sim parameter
  # extrinsic to daily block utilization, so blocks from sims that differ only
  # in RVU target are comparable and can be pooled together for stats analysis.
  
  # 1. specify sim data that can be pooled:
  specs <- expand.grid(.specialty = .specialty,
                       .cluster_size = .cluster_size,
                       .casemix = .casemix,
                       .rvu_goal = seq(3000,12000,1000),
                       .turnover_time = .turnover_time,
                       .block_hours = .block_hours,
                       stringsAsFactors = FALSE)
  
  # 2. import the data to be pooled. Cannot use future_pmap since 'row.prepare'
  #    is being called within future_pmap (and hence already limited to one processor)
  data.raw <- pmap(specs, import.sim) %>%
    bind_rows() %>% 
    rename(time.actual = time.or) %>%
    rename(rvus.perblock = rvus) %>% 
    mutate(time.turnover = (n.cases-1)*turnover.time_spec)

  # 3. the below is collection of data by surgeon a.k.a annual schedule,
  #    rather than by block. Linear regression will be perfomed on the  
  #    results of surgeon-level summaries
  data.grouped <- data.raw %>% 
    #amount of time spent over block length
    mutate(ot.actual = pmax(time.actual - block.size_spec,0)) %>%
    #amount of OR time requested (reserved + overtime)
    mutate(OR.actual = pmax(time.actual, block.size_spec)) %>%
    #group but do not use rvu_target as a grouping variable
    group_by(trial, specialty_spec, cluster.size_spec,
             casemix_spec, turnover.time_spec, block.size_spec) %>%
    summarize(n.blocks = max(block),
              rvu.total = sum(rvus.perblock),
              n.cases = sum(n.cases),
              #block.size_spec is instantiated during file import (see import.sim)
              p.blocks.ot.actual = sum(time.actual > block.size_spec) / n.blocks,
              ot.actual.percent_util = sum(ot.actual)/sum(time.actual),
              #mean(block.size_spec) is the one value of block.size_spec for the group
              ot.actual.percent_resv = sum(ot.actual)/(n.blocks*mean(block.size_spec)),
              ot.actual.percent_admin = sum(ot.actual)/sum(OR.actual),
              rvus.hourly.util = sum(rvus.perblock)/sum(time.actual)*60,
              rvus.hourly.resv = sum(rvus.perblock)/(n.blocks*mean(block.size_spec))*60,
              rvus.hourly.admin =
                sum(rvus.perblock)/(sum(pmax(time.actual,block.size_spec)))*60) %>% 
    ungroup()
  
  # 4. fit seven linear models:
  #       % of blocks that went overtime ~ blocksize + turnover + rvu_target 
  #       proportion of overtime to utilized OR time ~ ...
  #       proportion of overtime to reserved OR time ~ ...
  #       proportion of overtime to administrative OR time ~ ...
  #       hourly rvus (per utilization) ~ ...
  #       hourly rvus (per reservation) ~ ...
  #       hourly rvus (per admin cost) ~ ...
  
  #       annual overtime and turnover time must be pegged to RVU target
  #       rates (percent..., hourly rvus...) can pool data from all RVU targets
  
  p.ot.actual.tbl <-
    data.grouped %>% coeff.table(p.blocks.ot.actual, .regress_by_rvu = FALSE)
  ot.actual.p.util.tbl <- 
    data.grouped %>% coeff.table(ot.actual.percent_util, .regress_by_rvu = FALSE)
  ot.actual.p.resv.tbl <- 
    data.grouped %>% coeff.table(ot.actual.percent_resv, .regress_by_rvu = FALSE)
  ot.actual.p.admin.tbl <-
    data.grouped %>% coeff.table(ot.actual.percent_admin, .regress_by_rvu = FALSE)
  rvus.hourly.util.tbl <- 
    data.grouped %>% coeff.table(rvus.hourly.util, .regress_by_rvu = FALSE,
                                 .invert_blocksize = TRUE, .invert_turnover = TRUE)
  rvus.hourly.resv.tbl <-
    data.grouped %>% coeff.table(rvus.hourly.resv, .regress_by_rvu = FALSE,
                                 .invert_blocksize = TRUE, .invert_turnover = TRUE)
  rvus.hourly.admin.tbl <-
    data.grouped %>% coeff.table(rvus.hourly.admin, .regress_by_rvu = FALSE,
                                 .invert_blocksize = TRUE, .invert_turnover = TRUE)
  
  
  # 5. format and return output row with identifier columns at left and p values at right
  output <- list(p.ot.actual = p.ot.actual.tbl,
                 ot.actual.p_util = ot.actual.p.util.tbl,
                 ot.actual.p_resv = ot.actual.p.resv.tbl,
                 ot.actual.p_admin = ot.actual.p.admin.tbl,
                 rvus.hourly_util = rvus.hourly.util.tbl,
                 rvus.hourly_resv = rvus.hourly.resv.tbl,
                 rvus.hourly_admin = rvus.hourly.admin.tbl)
  return(output)
}

####  SUBR: OPERATING TIME  ####

# Independent of all factors except for specialty and casemix

row.prepare.optime <- function(.specialty, .casemix){
  
  # 1. specify sim data that can be pooled:
  specs <- expand.grid(.specialty = .specialty,
                       .cluster_size = c(5,10,20,Sched.ANNUAL),
                       .casemix = .casemix,
                       .rvu_goal = seq(3000,12000,1000),
                       .turnover_time = seq(0,90,10),
                       .block_hours = seq(4,12,2),
                       stringsAsFactors = FALSE)
  
  # 2. import the data to be pooled. Cannot use future_pmap since 'row.prepare' is being
  #    called within future_pmap (and hence is already limited to one processor)
  data.raw <- pmap(specs, import.sim) %>% #use same method for ot and rvus 
    bind_rows() %>% 
    rename(time.actual = time.or) %>%
    rename(rvus.perblock = rvus) %>% 
    mutate(time.turnover = (n.cases-1)*turnover.time_spec)
  
  # 3. summarize data
  data.grouped <- data.raw %>% 
    group_by(trial, specialty_spec, cluster.size_spec,
             casemix_spec, turnover.time_spec, block.size_spec, rvu.target_spec) %>%
    summarize(n.cases = sum(n.cases),
              rvu.total = sum(rvus.perblock),
              time.actual.net = sum(time.actual),
              time.turnover.net = sum(time.turnover),
              optime.net = (time.actual.net - time.turnover.net)/60) %>%
    ungroup() %>%
    select(-c(time.actual.net, time.turnover.net))
  
  # 4. compute regression
  # 'NA' is hacky way for specifying regression ONLY by RVU in coeff.table
  opt.tbl <- data.grouped %>% coeff.table(optime.net, .regress_by_rvu = NA)
  
  return(opt.tbl)
}

####  EXECUTE   ####
print(Sys.time())

gc(verbose = TRUE)

options(error=traceback)

workers.total = 14

##  ANNUAL OPERATING TIME  ##

# less workers because regressions require loading nearly all data at same time
plan(multiprocess, workers = 3)

# Array containing combinations of simulation parameters determining total operating time
calc.specs <- expand.grid(.specialty = speclist,
                          .casemix = casemixlist)

print("total optime")

results <-
  future_pmap(calc.specs, row.prepare.optime, .progress = TRUE) %>% bind_rows()

print(results)

print("extracting & writing")

write.csv(results,
          file = paste0(root_results,"linear_models_optime.csv"),
          row.names = FALSE)

print(Sys.time())

##  ANNUAL BLOCK REQUIREMENTS  ##

# Array containing combinations of simulation parameters.
# Each row represents a simulation for which to include in linear model calculations.
calc.specs <- expand.grid(.block_hours = c(4,6,8,10,12,NA),
                          .specialty = speclist,
                          .casemix = casemixlist,
                          .cluster_size = c(5,10,20,Sched.ANNUAL),
                          .turnover_time = c(0,10,20,30,40,50,60,70,80,90,NA),
                          stringsAsFactors = FALSE)

print("blockreqs")

plan(multiprocess, workers = workers.total)

results <-
  future_pmap(calc.specs, row.prepare.blockreqs, .progress = TRUE) %>% bind_rows()

print("extracting & writing")

write.csv(results,
          file = paste0(root_results,"linear_models_blockreqs.csv"),
          row.names = FALSE)

print(Sys.time())

##  ANNUAL COUNTABLE METRICS (NET TURNOVER, NET OVERTIME, NUMBER CASES, ETC.)  ##

print("toot")

plan(multiprocess, workers = workers.total)

results <- future_pmap(calc.specs, row.prepare.toot, .progress = TRUE) %>%
  unlist(recursive = FALSE)

print("extracting")

res.ot <- results %>% .[seq(1,length(results),5)] %>% bind_rows()
res.to <- results %>% .[seq(2,length(results),5)] %>% bind_rows()
res.ncases <- results %>% .[seq(3,length(results),5)] %>% bind_rows()
res.timeor <- results %>% .[seq(4,length(results),5)] %>% bind_rows()
res.rvutotal <- results %>% .[seq(5,length(results),5)] %>% bind_rows()

print("writing")

write.csv(res.ot,
          file = paste0(root_results,"linear_models_overtime.csv"),
          row.names = FALSE)
write.csv(res.to,
          file = paste0(root_results,"linear_models_turnover.csv"),
          row.names = FALSE)
write.csv(res.ncases,
          file = paste0(root_results,"linear_models_ncases.csv"),
          row.names = FALSE)
write.csv(res.timeor,
          file = paste0(root_results,"linear_models_timeor.csv"),
          row.names = FALSE)
write.csv(res.rvutotal,
          file = paste0(root_results,"linear_models_rvutotal.csv"),
          row.names = FALSE)

print(Sys.time())

##  RATE METRICS (% OVERTIME & HOURLY RVUS)  ##
print("pothrvu")

plan(multiprocess, workers = workers.total)

results <- future_pmap(calc.specs, row.prepare.pothrvu, .progress = TRUE) %>%
  unlist(recursive = FALSE)

print("extracting")

res.p.ot.actual <- results %>% .[seq(1,length(results),7)] %>% bind_rows()
res.ot.actual.p_util <- results %>% .[seq(2,length(results),7)] %>% bind_rows()
res.ot.actual.p_resv <- results %>% .[seq(3,length(results),7)] %>% bind_rows()
res.ot.actual.p_admin <- results %>% .[seq(4,length(results),7)] %>% bind_rows()
res.hrvus_util <- results %>% .[seq(5,length(results),7)] %>% bind_rows()
res.hrvus_resv <- results %>% .[seq(6,length(results),7)] %>% bind_rows()
res.hrvus_admin <- results %>% .[seq(7,length(results),7)] %>% bind_rows()

print("writing")

write.csv(res.p.ot.actual,
          file = paste0(root_results,"linear_models_rate_ot_actual.csv"),
          row.names = FALSE)
write.csv(res.ot.actual.p_util,
          file = paste0(root_results,"linear_models_ot_percent_util.csv"),
          row.names = FALSE)
write.csv(res.ot.actual.p_resv,
          file = paste0(root_results,"linear_models_ot_percent_resv.csv"),
          row.names = FALSE)
write.csv(res.ot.actual.p_admin,
          file = paste0(root_results,"linear_models_ot_percent_admin.csv"),
          row.names = FALSE)
write.csv(res.hrvus_util,
          file = paste0(root_results,"linear_models_hrvus_utilized.csv"),
          row.names = FALSE)
write.csv(res.hrvus_resv,
          file = paste0(root_results,"linear_models_hrvus_reserved.csv"),
          row.names = FALSE)
write.csv(res.hrvus_admin,
          file = paste0(root_results,"linear_models_hrvus_admin.csv"),
          row.names = FALSE)

print(Sys.time())
print("done")
#!/usr/bin/env Rscript

##################################################
## Project: RVU Target Feasibility Study
## Script purpose: Model heteroscedasticity in outcome
##                 metrics by fitting standard deviations
##                 as functions of simulation parameters
## Date: 12/25/2020
## Author: Saieesh Rao
##################################################

####  LIBRARIES   ####
library(dplyr)
library(magrittr)
library(purrr)
library(future)
library(furrr)
library(data.table)
library(sandwich)
library(lmtest)

####  CONSTANTS  ####

on.hpc = TRUE  #change to 'TRUE' once running on Gardner HPC

if(on.hpc){
  root = "/scratch/t.sur.rsaieesh/results-binpack/"
} else {
  root = "X:/results-binpack/"
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

casemixlist <- c(Casemix.ALL, Casemix.HIGH, Casemix.MEDIUM, Casemix.LOW)

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

# list of specialties in order, with ir removed
speclist <- c(cardiac, general, gyne, neuro, ortho,
              ent, plastics, thoracic, urology, vascular)

##  READ DATA  ##
alldata <- fread(paste0(root,"summary_table_binpack.csv"))

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
  #print(names(.data))
  var <- substitute(.var)
  .data %<>%
    mutate(block_size = block_size / 60) %>% #convert mins to hours for regression
    rename(rvu_target = rvu_goal)
  
  # algorithmically generate regression formulas from coeff.table arguments.
  # formulas returned as strings and reconstituted as formulas in appropriate environment
  # order of formula creation is RVU, Block Size, Turnover time (order is important later)
  synth.formula <- function(x.regress_by_rvu, x.invert_turnover, x.invert_blocksize) {
    
    #formula base
    fx <- 'eval(var) ~ '
    
    if(is.na(x.regress_by_rvu)) {
      #only regress by RVU; can return formula now
      #invert arguments are irrelevant in this case
      fx_rvu <- 'rvu_target'
      return(paste0(fx,fx_rvu))
      
    } else {
      if(x.regress_by_rvu) {
        #only inlcude rvu var if TRUE, otherwise default blank
        fx_rvu <- 'rvu_target + '
      } else {
        fx_rvu <- ''
      }
      #blocksize and turnover time are always included; variation is
      #whether to invert values or use raw values directly
      if(x.invert_blocksize) {
        fx_bls <- 'I(block_size^-1) + '
      } else {
        fx_bls <- 'block_size + '
      }
      if(x.invert_turnover) {
        fx_tot <- 'I(turnover_time^-1)'
      } else {
        fx_tot <- 'turnover_time'
      }
      
      #combine formula terms as string
      fx <- paste0(fx, paste0(fx_rvu,fx_bls,fx_tot))
      
      return(fx)
    }
  }
  
  formula <- synth.formula(x.regress_by_rvu = .regress_by_rvu,
                            x.invert_turnover = .invert_turnover,
                            x.invert_blocksize = .invert_blocksize)
  
  # remove turnover time == 0 if inverting turnover_time (can't divide by zero)
  # similar logic applies for inverting block_size, but no sims have block_size == 0
  if(.invert_turnover){
    .data %<>% filter(turnover_time > 0)
    
    #exit if no data because turnover_time == 0 for all data
    if(.data %>% dim() %>% .[1] == 0){
      return(NULL) #will result in dropped row during row binding
    }
  }
  
  # calculate ordinary least-squares (OLS) linear models with given formula
  model_target <- lm(as.formula(formula), .data)
  
  # obtain model parameters and related values, depending on form of formula
  if(is.na(.regress_by_rvu) || .regress_by_rvu){
    
    #additional weighted least-squares (WLS) regression, weighting by RVU metric
    wls_target <- lm(as.formula(formula), data = .data, weights = .data$rvu_target^-1)

    # OLS model using RVU target
    coeffs_target <- model_target$coefficients
    summary_target <- summary(model_target)
    regvars_target <- summary_target$coefficients %>% row.names()
    #heteroscedasticity consistent estimators
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
    
  }
  else{ #.regress_by_rvu == FALSE
    #formula_target and formula_total are identical when .regress_by_rvu == FALSE
    #(there's no RVU term in the formula) so just pick one for model calculations.
    #additionally no need for weighting by RVU metric if not regressing by RVUs.
    
    coeffs <- model_target$coefficients
    summary <- summary(model_target)
    regvars <- summary$coefficients %>% row.names()
    #heteroscedasticity consistent estimators
    hce <- model_target %>% coeftest(vcov=vcovHC(.,type='HC3')) 
    #for each of below, [1='(Intercept)', 2=rvu_target,
    #3=block_size OR I(block_size^-1), 4=turnover_time OR I(turnover_time^-1)]
    hce_tv <- hce[,'t value']
    hce_pv <- hce[,'Pr(>|t|)']
    hce_se <- hce[,'Std. Error']
    AIC_model <- AIC(model_target)
    
  }
  
  .specialty <- .data$specialty %>% unique()         #should be one value
  .cluster_size <- .data$cluster_size %>% unique()   #should be one value
  .casemix <- .data$casemix %>% unique()             #should be one value
  .turnover_time <- .data$turnover_time %>% unique()
  .block_hours <- .data$block_size %>% unique()
  .rvu_target <- .data$rvu_target %>% unique()

  #if list OR null, need to convert to single value in output table
  if(length(.specialty) != 1) .specialty <- NA         
  if(length(.cluster_size) != 1) .cluster_size <- NA  
  if(length(.casemix) != 1) .casemix <- NA            
  if(length(.block_hours) != 1) .block_hours <- NA     
  if(length(.turnover_time) != 1) .turnover_time <- NA
  if(length(.rvu_target) != 1) .rvu_target <- NA
  
  # if regression contains RVU as independent var
  # n.b. RVU total regression results were omitted in final analysis since
  #      unclear how to calculate SD for RVU totals with few obs. at each unique total
  if(is.na(.regress_by_rvu) || .regress_by_rvu){ 
    
    output_target <- tribble(
      ~surgspec, ~cluster_size, ~casemix, ~turnover_time, ~block_hours,
      ~rvu_target, ~invert.turnover_time, ~invert.block_hours,
      ~coeff.rvu_target, ~coeff.block_hours, ~coeff.turnover_time, ~intercept,
      #heteroscedasticity consistent standard errors
      ~HCSE.rvu_target, ~HCSE.block_hours, ~HCSE.turnover_time, ~HCSE.intercept,
      ~tvalue.rvu_target, ~tvalue.block_hours, ~tvalue.turnover_time, ~tvalue.intercept,
      ~pvalue.rvu_target, ~pvalue.block_hours, ~pvalue.turnover_time, ~pvalue.intercept,
      ~rsquared.mult, ~rsquared.adjust,
      ~stderr.regression,~AIC_target,
      ~coeff.WLS.rvu_target, ~coeff.WLS.block_hours,
      ~coeff.WLS.turnover_time, ~intercept.WLS,
      ~SE.WLS.rvu_target, ~SE.WLS.block_hours,
      ~SE.WLS.turnover_time, ~SE.WLS.intercept,
      ~tvalue.WLS.rvu_target, ~tvalue.WLS.block_hours,
      ~tvalue.WLS.turnover_time, ~tvalue.WLS.intercept,
      ~pvalue.WLS.rvu_target, ~pvalue.WLS.block_hours,
      ~pvalue.WLS.turnover_time, ~pvalue.WLS.intercept,
      ~rsquared.mult.WLS, ~rsquared.adjust.WLS,
      ~stderr.regression.WLS,~AIC_target.WLS,
      
      .specialty, .cluster_size, .casemix, .turnover_time, .block_hours,
      .rvu_target, .invert_turnover, .invert_blocksize,
      coeffs_target["rvu_target"], NA,
      NA, coeffs_target["(Intercept)"],
      NA,NA,
      NA, hce_target_se["(Intercept)"],
      NA,NA,
      NA, hce_target_tv["(Intercept)"],
      NA,NA,
      NA, hce_target_pv["(Intercept)"],
      # NA,NA,NA,NA,summary_target$coefficients["(Intercept)",4],
      summary_target$r.squared, summary_target$adj.r.squared,
      summary_target$sigma, AIC_target,
      wls_target_coeff["rvu_target"], NA, NA, wls_target_coeff["(Intercept)"],
      NA,NA,NA, wls_target_se["(Intercept)"],
      NA,NA,NA, wls_target_tv["(Intercept)"],
      NA,NA,NA, wls_target_pv["(Intercept)"],
      summary_wls_target$r.squared, summary_wls_target$adj.r.squared,
      summary_wls_target$sigma, AIC_wls_target
    )
    
    #add rvu associated coeffs, since these regressions have rvus as independent variables
    output_target$HCSE.rvu_target <- hce_target_se["rvu_target"]
    # output_total$HCSE.rvu_total <- hce_total_se["rvu.total"]
    output_target$tvalue.rvu_target <- hce_target_tv["rvu_target"]
    # output_total$tvalue.rvu_total <- hce_total_tv["rvu.total"]
    output_target$pvalue.rvu_target <- hce_target_pv["rvu_target"]
    # output_total$pvalue.rvu_total <- hce_total_pv["rvu.total"]
    
    output_target$SE.WLS.rvu_target <- wls_target_se["rvu_target"]
    # output_total$SE.WLS.rvu_total <- wls_total_se["rvu.total"]
    output_target$tvalue.WLS.rvu_target <- wls_target_tv["rvu_target"]
    # output_total$tvalue.WLS.rvu_total <- wls_total_tv["rvu.total"]
    output_target$pvalue.WLS.rvu_target <- wls_target_pv["rvu_target"]
    # output_total$pvalue.WLS.rvu_total <- wls_total_pv["rvu.total"]
    
    #if regressed by RVU + things besides RVU
    if(!is.na(.regress_by_rvu)){
      
      #if block.size_spec was a regression variable:
      if(TRUE %in%
         (c('I(block_size^-1)','block_size') %in% regvars_target)){
        
        # missing name returns zero, so pick bigger number (non-zero) as index
        bs_ind <- max(regvars_target %>% {which(.=='block_size')},
                      regvars_target %>% {which(.=='I(block_size^-1)')})
        
        output_target$coeff.block_hours <- coeffs_target[bs_ind]
        # output_total$coeff.block_hours <- coeffs_total[bs_ind]
        output_target$HCSE.block_hours <- hce_target_se[bs_ind]
        # output_total$HCSE.block_hours <- hce_total_se[bs_ind]
        output_target$tvalue.block_hours <- hce_target_tv[bs_ind]
        # output_total$tvalue.block_hours <- hce_total_tv[bs_ind]
        output_target$pvalue.block_hours <- hce_target_pv[bs_ind]
        # output_total$pvalue.block_hours <- hce_total_pv[bs_ind]
        
        output_target$coeff.WLS.block_hours <- wls_target_coeff[bs_ind]
        # output_total$coeff.WLS.block_hours <- wls_total_coeff[bs_ind]
        output_target$SE.WLS.block_hours <- wls_target_se[bs_ind]
        # output_total$SE.WLS.block_hours <- wls_total_se[bs_ind]
        output_target$tvalue.WLS.block_hours <- wls_target_tv[bs_ind]
        # output_total$tvalue.WLS.block_hours <- wls_total_tv[bs_ind]
        output_target$pvalue.WLS.block_hours <- wls_target_pv[bs_ind]
        # output_total$pvalue.WLS.block_hours <- wls_total_pv[bs_ind]
        
      }
      #if turnover.time_spec was a regression variable:
      if(TRUE %in%
         (c('I(turnover_time^-1)','turnover_time') %in% regvars_target)){
        
        # missing name returns zero, so pick bigger number (non-zero) as index
        # regvars list for both rvu_target and rvu_total are identical for BS and TO
        to_ind <- max(regvars_target %>% {which(.=='turnover_time')},
                      regvars_target %>% {which(.=='I(turnover_time^-1)')})
        
        output_target$coeff.turnover_time <- coeffs_target[to_ind]
        # output_total$coeff.turnover_time <- coeffs_total[to_ind]
        output_target$HCSE.turnover_time <- hce_target_se[to_ind]
        # output_total$HCSE.turnover_time <- hce_total_se[to_ind]
        output_target$tvalue.turnover_time <- hce_target_tv[to_ind]
        # output_total$tvalue.turnover_time <- hce_total_tv[to_ind]
        output_target$pvalue.turnover_time <- hce_target_pv[to_ind]
        # output_total$pvalue.turnover_time <- hce_total_pv[to_ind]
        
        output_target$coeff.WLS.turnover_time <- wls_target_coeff[to_ind]
        # output_total$coeff.WLS.turnover_time <- wls_total_coeff[to_ind]
        output_target$SE.WLS.turnover_time <- wls_target_se[to_ind]
        # output_total$SE.WLS.turnover_time <- wls_total_se[to_ind]
        output_target$tvalue.WLS.turnover_time <- wls_target_tv[to_ind]
        # output_total$tvalue.WLS.turnover_time <- wls_total_tv[to_ind]
        output_target$pvalue.WLS.turnover_time <- wls_target_pv[to_ind]
        # output_total$pvalue.WLS.turnover_time <- wls_total_pv[to_ind]
      }
    }
    
    output <- output_target
    
  }
  else{ #.regress_by_rvu == NA
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
      coeffs["block_size"], coeffs["turnover_time"], coeffs["(Intercept)"],
      NA,NA, hce_se["(Intercept)"],
      NA,NA, hce_tv["(Intercept)"],
      NA,NA, hce_pv["(Intercept)"],
      # NA,NA,summary$coefficients["(Intercept)",4],
      summary$r.squared, summary$adj.r.squared, summary$sigma, AIC_model
    )
    
    #if block.size_spec was a regression variable:
    if(TRUE %in% (c('I(block_size^-1)','block_size') %in% regvars)){
      
      # missing name returns zero, so pick bigger number (non-zero) as index
      bs_ind <- max(regvars %>% {which(.=='block_size')},
                    regvars %>% {which(.=='I(block_size^-1)')})
      
      output$coeff.block_hours <- coeffs[bs_ind]
      output$HCSE.block_hours <- hce_se[bs_ind]
      output$tvalue.block_hours <- hce_tv[bs_ind]
      output$pvalue.block_hours <- hce_pv[bs_ind]
    }
    #if turnover.time_spec was a regression variable:
    if(TRUE %in% (c('I(turnover_time^-1)','turnover_time') %in% regvars)){
      
      # missing name returns zero, so pick bigger number (non-zero) as index
      to_ind <- max(regvars %>% {which(.=='turnover_time')},
                    regvars %>% {which(.=='I(turnover_time^-1)')})
      
      output$coeff.turnover_time <- coeffs[to_ind]
      output$HCSE.turnover_time <- hce_se[to_ind]
      output$tvalue.turnover_time <- hce_tv[to_ind]
      output$pvalue.turnover_time <- hce_pv[to_ind]
    }
  }
  
  #remove block hour and turnover time related columns if not regression variables
  if(is.na(.regress_by_rvu)){
    output %<>% select(-c(
      coeff.block_hours,coeff.turnover_time,
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


####  SUBR: MULTIVARIATE REGRESSION OF STANDARD DEVIATIONS  ####

# regress standard deviations of several variables against RVUs

row.prepare.sd <- function(.specialty, .cluster_size, .casemix,
                                  .turnover_time=seq(0,90,10), .block_hours=seq(4,12,2)){
  
  if(length(.turnover_time) == 1 && is.na(.turnover_time)) {
    .turnover_time = seq(0,90,10)
  }
  if(length(.block_hours) == 1 && is.na(.block_hours)) {
    .block_hours = seq(4,12,2)
  }
  
  .subset <- alldata %>%
    filter(specialty %in% .specialty) %>%
    filter(cluster_size %in% .cluster_size) %>%
    filter(casemix %in% .casemix) %>%
    filter(turnover_time %in% .turnover_time) %>%
    filter(block_size %in% (.block_hours*60)) %>%
    filter(!is.na(rvu_goal))

  n.blocks.sd.tbl <- .subset %>%
    coeff.table(n.blocks.sd, .regress_by_rvu = TRUE, .invert_blocksize = TRUE)
  n.cases.sd.tbl <- .subset %>%
    coeff.table(n.cases.sd, .regress_by_rvu = NA)
  turnover.sd.tbl <- .subset %>%
    coeff.table(time.turnover.net.sd, .regress_by_rvu = TRUE)
  optime.sd.tbl <- .subset %>%
    coeff.table(optime.actual.net.sd, .regress_by_rvu = NA)
  ot.actual.sd <- .subset %>%
    coeff.table(ot.actual.net.sd, .regress_by_rvu = TRUE)
  timeor.sd.tbl <- .subset %>%
    coeff.table(time.actual.sd, .regress_by_rvu = TRUE)
  rvusperblock.sd.tbl <- .subset %>%
    coeff.table(rvus.perblock.sd, .regress_by_rvu = FALSE)
  rateotactual.sd.tbl <- .subset %>%
    coeff.table(p.blocks.ot.actual.sd, .regress_by_rvu = FALSE)
  otpresv.sd.tbl <- .subset %>%
    coeff.table(ot.actual.percent_resv.sd, .regress_by_rvu = FALSE)
  otputil.sd.tbl <- .subset %>%
    coeff.table(ot.actual.percent_util.sd, .regress_by_rvu = FALSE)
  otpadmin.sd.tbl <- .subset %>%
    coeff.table(ot.actual.percent_admin.sd, .regress_by_rvu = FALSE)
  timeresv.p.util.sd.tbl <- .subset %>%
    coeff.table(time.resv.percent_util.sd, .regress_by_rvu = FALSE)
  timeor.p.resv.sd.tbl <- .subset %>%
    coeff.table(time.or.percent_resv.sd, .regress_by_rvu = FALSE)
  timeor.p.admin.sd.tbl <- .subset %>%
    coeff.table(time.or.percent_admin.sd, .regress_by_rvu = FALSE)
  hrvus.util.sd.tbl <- .subset %>%
    coeff.table(rvus.hourly.util.sd, .regress_by_rvu = FALSE,
                .invert_blocksize = TRUE, .invert_turnover = TRUE)
  hrvus.resv.sd.tbl <- .subset %>%
    coeff.table(rvus.hourly.resv.sd, .regress_by_rvu = FALSE,
                .invert_blocksize = TRUE, .invert_turnover = TRUE)
  hrvus.admin.sd.tbl <- .subset %>%
    coeff.table(rvus.hourly.admin.sd, .regress_by_rvu = FALSE,
                .invert_blocksize = TRUE, .invert_turnover = TRUE)
  
  result <- list(n.blocks.sd.tbl,
                 n.cases.sd.tbl,
                 turnover.sd.tbl,
                 optime.sd.tbl,
                 ot.actual.sd,
                 timeor.sd.tbl,
                 rvusperblock.sd.tbl,
                 rateotactual.sd.tbl,
                 otpresv.sd.tbl,
                 otputil.sd.tbl,
                 otpadmin.sd.tbl,
                 timeresv.p.util.sd.tbl,
                 timeor.p.resv.sd.tbl,
                 timeor.p.admin.sd.tbl,
                 hrvus.util.sd.tbl,
                 hrvus.resv.sd.tbl,
                 hrvus.admin.sd.tbl)
  
  return(result)
}

####  EXECUTE  ####

plan(multiprocess, workers=8)

calc.specs <- expand.grid(.block_hours = c(4,6,8,10,12,NA),
                          .specialty = speclist,
                          .cluster_size = c(5,10,20,Sched.ANNUAL),
                          .casemix = casemixlist,
                          .turnover_time = c(0,10,20,30,40,50,60,70,80,90,NA),
                          stringsAsFactors = FALSE)

results <- future_pmap(calc.specs, row.prepare.sd, .progress = TRUE) %>%
  unlist(recursive = FALSE)

# extract collated result into data frames specific to a single metric
wrap <- 17 #number of different types of results to be extracted

n.blocks.sd.tbl <- results %>% .[seq(1,length(results),wrap)] %>% bind_rows()
n.cases.sd.tbl <- results %>% .[seq(2,length(results),wrap)] %>% bind_rows()
turnover.sd.tbl <- results %>% .[seq(3,length(results),wrap)] %>% bind_rows()
optime.sd.tbl <- results %>% .[seq(4,length(results),wrap)] %>% bind_rows()
ot.actual.sd.tbl <- results %>% .[seq(5,length(results),wrap)] %>% bind_rows()
time.or.sd.tbl <- results %>% .[seq(6,length(results),wrap)] %>% bind_rows()
rvus.perblock.sd.tbl <- results %>% .[seq(7,length(results),wrap)] %>% bind_rows()
p.ot.actual.sd.tbl <- results %>% .[seq(8,length(results),wrap)] %>% bind_rows()
ot.actual.percent_resv <- results %>% .[seq(9,length(results),wrap)] %>% bind_rows()
ot.actual.percent_util <- results %>% .[seq(10,length(results),wrap)] %>% bind_rows()
ot.actual.percent_admin <- results %>% .[seq(11,length(results),wrap)] %>% bind_rows()
time.resv.percent_util <- results %>% .[seq(12,length(results),wrap)] %>% bind_rows()
time.or.percent_resv <- results %>% .[seq(13,length(results),wrap)] %>% bind_rows()
time.or.percent_admin <- results %>% .[seq(14,length(results),wrap)] %>% bind_rows()
hrvus.util <- results %>% .[seq(15,length(results),wrap)] %>% bind_rows()
hrvus.resv <- results %>% .[seq(16,length(results),wrap)] %>% bind_rows()
hrvus.admin <- results %>% .[seq(17,length(results),wrap)] %>% bind_rows()

# save results in individual spreadsheets
write.csv(n.blocks.sd.tbl,
          file = paste0(root,"linear_models_blockreqs_sd.csv"),
          row.names = FALSE)
write.csv(n.cases.sd.tbl,
          file = paste0(root,"linear_models_ncases_sd.csv"),
          row.names = FALSE)
write.csv(turnover.sd.tbl,
          file = paste0(root,"linear_models_turnover_sd.csv"),
          row.names = FALSE)
write.csv(optime.sd.tbl,
          file = paste0(root,"linear_models_optime_sd.csv"),
          row.names = FALSE)
write.csv(ot.actual.sd.tbl,
          file = paste0(root,"linear_models_overtime_sd.csv"),
          row.names = FALSE)
write.csv(time.or.sd.tbl,
          file = paste0(root,"linear_models_timeor_sd.csv"),
          row.names = FALSE)
write.csv(rvus.perblock.sd.tbl,
          file = paste0(root,"linear_models_rvusperblock_sd.csv"),
          row.names = FALSE)
write.csv(p.ot.actual.sd.tbl,
          file = paste0(root,"linear_models_rate_ot_actual_sd.csv"),
          row.names = FALSE)
write.csv(ot.actual.percent_resv,
          file = paste0(root,"linear_models_ot_percent_reserved_sd.csv"),
          row.names = FALSE)
write.csv(ot.actual.percent_util,
          file = paste0(root,"linear_models_ot_percent_utilized_sd.csv"),
          row.names = FALSE)
write.csv(ot.actual.percent_admin,
          file = paste0(root,"linear_models_ot_percent_admin_sd.csv"),
          row.names = FALSE)
write.csv(time.resv.percent_util,
          file = paste0(root,"linear_models_timeor_percent_utilized_sd.csv"),
          row.names = FALSE)
write.csv(time.or.percent_resv,
          file = paste0(root,"linear_models_timeor_percent_reserved_sd.csv"),
          row.names = FALSE)
write.csv(time.or.percent_admin,
          file = paste0(root,"linear_models_timeor_percent_admin_sd.csv"),
          row.names = FALSE)
write.csv(hrvus.util,
          file = paste0(root,"linear_models_hrvus_utilized_sd.csv"),
          row.names = FALSE)
write.csv(hrvus.resv,
          file = paste0(root,"linear_models_hrvus_reserved_sd.csv"),
          row.names = FALSE)
write.csv(hrvus.admin,
          file = paste0(root,"linear_models_hrvus_admin_sd.csv"),
          row.names = FALSE)

print('done')
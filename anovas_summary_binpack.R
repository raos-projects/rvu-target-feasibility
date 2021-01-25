#!/usr/bin/env Rscript

# ANOVA of select outcome measures (optime, turnover time, OR time,
# and number of blocks) vs each of the six simulation parameters.
# Results reported with Tukey HSD correction

####  LIBRARIES   ####
library(dplyr)
library(magrittr)
library(tibble)
library(tidyr)  # drop_na
library(purrr)
library(future)
library(furrr)
library(data.table)

on.hpc <- TRUE   #change to true on HPC

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

root_result <- '/scratch/t.sur.rsaieesh/results-binpack/'

####  SUBR: SUMMARIZE ANOVA RESULTS  ####

# Report the percentage of p-values that are NOT statistically significant,
# suggesting that there is no significant difference between simulations.
# Although unusual, it is helpful to report the percentage because there
# are too many p-values to report individually and there is not a clear
# pattern aside from higher RVU targets being associated with more significant
# differences.

summarize.anova.results <- function(.metric, .indepvar){
  
  #generate anova filename algorithmically 
  filename <- paste0(root_result,'anova_',.metric,'_',.indepvar,'.csv')
  
  #read anova data (p-values)
  anova.data <- fread(filename)
  
  #get non-param column names; first 6 columns are SPECs and thus dropped
  cols <- names(anova.data) %>% .[7:length(.)]
  
  params <- c("specialty")
  
  if(.indepvar == 'byspecialty') {
    params <- c('specialty','rvu_target')
  }

  # count p-values > 0.05 (Tukey HSD significance level)
  anova.result <-
    anova.data %>%
    filter(cluster_size != Sched.ANNUAL) %>%
    # indepvar column is all 'NA' in anova data,
    # hence indepvar will be ignored in grouping
    group_by_at(params) %>%
    mutate_all(function(x) ifelse(x>0.05,1,0)) %>%
    summarise_at(cols,mean)
  
  write.csv(anova.result, paste0(root_result,'summary_anova_',
                                 .metric,'_',.indepvar,'.csv'), row.names = FALSE)
}

####  EXECUTE   ####

# file parameters
metrics <- c('blockreqs','optime','time-actual','turnover-net')
indepvar <- c('byblocksize','bycasemix',
              'byrvutarget','byschedule','byspecialty','byturnover')
file.specs <- expand.grid(.metric = metrics, .indepvar = indepvar)

print('Starting')
print(Sys.time())

file.specs %>% future_pmap(summarize.anova.results, .progress = TRUE)

print('Done')
print(Sys.time())
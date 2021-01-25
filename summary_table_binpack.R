#!/usr/bin/env Rscript

##################################################
## Project: RVU Target Feasibility Study
## Script purpose: Create Summary Table of All Sims
## Date: 12/12/2020
## Author: Saieesh Rao
##################################################

####  LIBRARIES   ####
library(dplyr)
library(magrittr)
library(purrr)
library(future)
library(furrr)
library(gtools)
library(data.table)

####  CONSTANTS  ####

on.hpc <- TRUE

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

# list of specialties; removed ir
speclist <- c(cardiac, general, gyne, neuro, ortho, ent,
              plastics, thoracic, urology, vascular)

##  SCHEDULE  ##

#no longer using scheduling patterns; instead, cluster of
#cases to be grouped at one time (5,10,20,30,99999) by the
#bin-packing algorithm. Sched.ANNUAL is an arbitrarily large
#constant meant to encompass all cases in an annual caselog.

Sched.ANNUAL = 99999

# Array containing combinations of simulation parameters.
# Each row represents a simulation for which to summarize .
calc.specs <- expand.grid(.block_hours = seq(from=4,to=12,by=2),
                          .specialty = speclist,
                          .cluster_size = c(5,10,20,Sched.ANNUAL),
                          .casemix = casemixlist,
                          .turnover_time = seq(from=0,to=90,by=10),
                          .rvu_goal = c(seq(from=3000,to=12000,by=1000),NA),
                          stringsAsFactors = FALSE)


####   SUBR: IMPORT FULL SIM RESULTS   ####

# defines a new 'import' function, which is used to read simulation data stored in
# the scratch folder. the filename is identified algorithmically from the parameters.

import <- function(.specialty, .rvu_goal,
                   .block_hours, .cluster_size, .casemix, .turnover_time) {
  
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
  data <- fread(path)  #read in trial data
  data %<>% mutate(rvu.target = .rvu_goal)  #mark rvu target assoc. with each simulation
  return(data)    #return data.frame containing imported trial data
}


####   SUBR: GENERIC SUMMARY STATS   ###

# calculate summary stats for a single variable in sim results.
# summary stats are calculated for the specified variable (.var).
# if variable is a percentage value, then set .sum=FALSE

summary.stats <- function(.data, .var, .sum=TRUE) {
  var <- substitute(.var)
  output <- .data %>% summarize(.min = min(eval(var, .data)),
                                '.10th' = quantile(eval(var, .data), probs=0.10, type=2),
                                '.25th' = quantile(eval(var, .data), probs=0.25, type=2),
                                '.50th' = median(eval(var, .data)),
                                '.75th' = quantile(eval(var, .data), probs=0.75, type=2),
                                '.90th' = quantile(eval(var, .data), probs=0.90, type=2),
                                .max = max(eval(var, .data)),
                                .mean = mean(eval(var, .data)),
                                .sd = sd(eval(var, .data)),
                                .sum = sum(eval(var, .data)))
  if(!.sum) output %<>% select(-.sum)
  output %<>% setNames(paste0(var, names(.))) #rename columns as var name + stat name
  return(output)
}

####  SUBR: PREPARE ROW OF SUMMARY STATS FOR ALL VARS  ####

# calculate summary stats for multiple variables associated with a single trial.
# result is returned as a table with a single row; these rows will be combined
# to generate the final summary table containing all sims.

row.prepare <- function(.specialty, .cluster_size,
                        .casemix, .turnover_time, .block_hours, .rvu_goal) {
  
  # If .rvu_goal = NA, then pool data from sims with same specialty, schedule, speed,
  # turnover time, and block size (differ only in RVU target). RVU target is a sim
  # parameter extrinsic to daily block utilization, so blocks from sims that differ
  # only in RVU target are comparable and can be pooled together for stats analysis.
  
  if(is.na(.rvu_goal)) .rvu_goal = seq(3000,12000,1000)
  
  ## 1. specify sim data that can be pooled:
  specs <- expand.grid(.specialty = .specialty,
                       .cluster_size = .cluster_size,
                       .casemix = .casemix,
                       .rvu_goal = .rvu_goal, #either one value or all values
                       .turnover_time = .turnover_time,
                       .block_hours = .block_hours,
                       stringsAsFactors = FALSE)
  
  ## 2. import the data to be pooled. Cannot use future_pmap since 'row.prepare' is being
  ##    called within future_pmap (and hence is already limited to one processor)
  data <- pmap(specs, import) %>%
    bind_rows()
  
  # calculate block_size in minutes for convenience
  block_size = .block_hours*60
  
  # more changes for convenience
  data %<>% 
    rename(time.actual = time.or) %>%
    rename(rvus.perblock = rvus) %>% 
    mutate(time.turnover = (n.cases-1)*UQ(.turnover_time))
    
  # create filtered datasets for only cases that were estimated overtime or undertime
  data.ot.actual <-
    data %>% filter(time.actual > block_size) %>%
    mutate(ot.actual = time.actual - UQ(block_size)) %>% dplyr::select(ot.actual)
  
  ## 3. start analysis. these first results are for an entire cohort of surgeons
  ##    whose schedules were simulated using the same sim parameters.
  ##    variables prefixed by 'n.' are summable values, whereas those prefixed by 'p.'
  ##    are proportions between two or more variables and therefore not summable.
  results <- data %>% summarise(n.blocks = n(),
                                n.blocks.ot.actual = sum(time.actual > block_size),
                                p.blocks.ot.actual = n.blocks.ot.actual / n.blocks) %>%
    bind_cols(list(summary.stats(data,time.actual),
                   summary.stats(data,rvus.perblock))
              )
  
  # calculate net overtime only in blocks which ran overtime
  results.ot.actual <- data.ot.actual %>% summary.stats(ot.actual)
  
  # now analyze data by surgeon a.k.a. annual schedule.
  # first, we summarize each of the 100 or 1000 surgeon schedules for set of sim
  # parameters, which is stored as a data.table (data.grouped).
  # then, we summarize the summaries of the surgeon schedules to condense surgeon-level 
  # metrics into a single summary row for each set of sim parameters (results.grouped).

  data.grouped <- data %>% 
    #amount of time spent over block length for each block; zero if no overtime
    mutate(ot.actual = pmax(time.actual - UQ(block_size),0)) %>%
    #amount of OR time requested (reserved + overtime)
    mutate(OR.actual = pmax(time.actual, UQ(block_size))) %>%
    group_by(trial,rvu.target) %>%
    summarize(#number of blocks used by a single surgeon
              n.blocks = max(block),
              #number of cases performed by a single surgeon
              n.cases = sum(n.cases),
              #percentage of blocks that went overtime for a single surgeon
              p.blocks.ot.actual = sum(time.actual > block_size) / n.blocks,
              #net time spent in turnovers by a single surgeon
              time.turnover.net = sum(time.turnover),
              #net time spent operating by a single surgeon
              optime.actual.net = (sum(time.actual)-sum(time.turnover)),
              #net overtime block utilization for a single surgeon
              ot.actual.net = sum(ot.actual),
              #proportion of overtime to overall OR utilization for a single surgeon
              ot.actual.percent_util = sum(ot.actual)/sum(time.actual),
              #proportion of overtime to reserved block time for a single surgeon
              ot.actual.percent_resv = sum(ot.actual)/(n.blocks*block_size),
              #proportion of overtime to administrative block time for a single surgeon
              #n.b. administrative time is at minimum block reservation, plus overtime
              ot.actual.percent_admin = sum(ot.actual)/sum(OR.actual),
              #total OR time, inclusive of turnovers, for a single surgeon
              time.or.net = sum(time.actual),
              #proportion of total OR time to reserved block time, for a single surgeon
              time.or.percent_resv = sum(time.actual)/(n.blocks*block_size),
              #proportion of total OR time to administrative time, for a single surgeon
              time.or.percent_admin = sum(time.actual)/sum(OR.actual),
              #proportion of reserved block time that was utilized, for a single surgeon
              time.resv.percent_util =
                (sum(time.actual)-sum(ot.actual))/(n.blocks*block_size),
              #rate of RVUs generated per hour of utilized block time for a single surgeon
              rvus.hourly.util = sum(rvus.perblock)/sum(time.actual)*60,
              #rate of RVUs generated per hour of reserved block time for a single surgeon
              rvus.hourly.resv = sum(rvus.perblock)/(n.blocks*block_size)*60,
              #rate of RVUs generated per hour of admin block time for a single surgeon
              rvus.hourly.admin = 
                sum(rvus.perblock)/(sum(pmax(time.actual,block_size)))*60) %>% 
    ungroup()

  results.grouped <-
    bind_cols(list(summary.stats(data.grouped, n.blocks),
                   summary.stats(data.grouped, n.cases),
                   summary.stats(data.grouped, p.blocks.ot.actual, .sum=FALSE),
                   summary.stats(data.grouped, time.turnover.net),
                   summary.stats(data.grouped, optime.actual.net),
                   summary.stats(data.grouped, ot.actual.net),
                   summary.stats(data.grouped, ot.actual.percent_util, .sum=FALSE),
                   summary.stats(data.grouped, ot.actual.percent_resv, .sum=FALSE),
                   summary.stats(data.grouped, ot.actual.percent_admin, .sum=FALSE),
                   summary.stats(data.grouped, time.or.net),
                   summary.stats(data.grouped, time.or.percent_resv, .sum=FALSE),
                   summary.stats(data.grouped, time.or.percent_admin, .sum=FALSE),
                   summary.stats(data.grouped, time.resv.percent_util, .sum=FALSE),
                   summary.stats(data.grouped, rvus.hourly.util, .sum=FALSE),
                   summary.stats(data.grouped, rvus.hourly.resv, .sum=FALSE),
                   summary.stats(data.grouped, rvus.hourly.admin, .sum=FALSE))
              )

  # 4. format and return output row with identifier columns at left and p values at right
  if(length(.rvu_goal) > 1) .rvu_goal <- NA   #convert [3000,12000] to NA for .csv
  output <- tribble(
    ~specialty, ~cluster_size, ~casemix, ~turnover_time, ~block_size, ~rvu_goal,
    .specialty, .cluster_size, .casemix, .turnover_time, block_size, .rvu_goal
  )
  
  output %<>% bind_cols(list(results,
                             results.ot.actual,
                             results.grouped))
  return(output)
}

####  EXECUTE   ####

print(Sys.time())

plan(multiprocess, workers = 8)

# calculate summary row for every sim, in parallel
output <- future_pmap(calc.specs, row.prepare, .progress = TRUE)

# stack rows to form table from rows
output %<>% bind_rows()

# save summary table and finish
write.csv(output,
          file = "/scratch/t.sur.rsaieesh/results-binpack/summary_table_binpack.csv",
          row.names = FALSE)

print(Sys.time())
print("done")
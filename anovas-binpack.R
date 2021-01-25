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

import <- function(.specialty, .rvu_target,
                   .block_hours, .cluster_size, .casemix, .turnover_time) {
  #systematically generate filename and path based on trial settings
  
  if(.cluster_size == Sched.ANNUAL){
    middle_path = 'binpack-global/bp_'
  }
  else{
    middle_path = 'binpack/bp_'
  }
  
  if(on.hpc){
    path <- paste0("/scratch/t.sur.rsaieesh/",middle_path,.specialty,
                   "_RVUs-",.rvu_target,
                   "_blocksize-",.block_hours,
                   "_casemix-",.casemix,
                   "_turnover-",.turnover_time,
                   "-cluster-",.cluster_size,
                   ".csv")
  } else {
    path <- paste0("X:/",middle_path,.specialty,
                   "_RVUs-",.rvu_target,
                   "_blocksize-",.block_hours,
                   "_casemix-",.casemix,
                   "_turnover-",.turnover_time,
                   "-cluster-",.cluster_size,
                   ".csv")
  }
  print(paste("importing",path))
  data <- read.csv(path)  #read in trial data
  data %<>% mutate(rvu_target = .rvu_target) %>%  #mark rvu target for each sim
    mutate(specialty = as.factor(.specialty)) %>%
    mutate(casemix = as.factor(.casemix)) %>%
    mutate(cluster_size = as.factor(.cluster_size)) %>%
    mutate(turnover_time = .turnover_time) %>%
    mutate(block_hours = .block_hours)
  return(data) #return data.frame containing requested trial data as output from function
}

####    SUBR: GENERIC ANOVA FUNCTION    ####

# generalizable one-way ANOVA function.
# independent var is passed as 'NA' in argument list
row.prepare.anovas <- function(.specialty, .cluster_size,
                               .casemix, .block_hours, .turnover_time, .rvu_target) {
  #1. Identify independent variable for ANOVA from parameter list
  if(length(.specialty) == 1 && is.na(.specialty)){
    .specialty = speclist
    var <- substitute(specialty)
  }
  if(length(.casemix) == 1 && is.na(.casemix)){
    .casemix = casemixlist
    var <- substitute(casemix)
  } 
  if(length(.cluster_size) == 1 && is.na(.cluster_size)){
    .cluster_size = c(5,10,20,Sched.ANNUAL)
    var <- substitute(cluster_size)
  }
  if(length(.turnover_time) == 1 && is.na(.turnover_time)){
    .turnover_time = seq(0,90,10)
    var <- substitute(turnover_time)
  }
  if(length(.block_hours) == 1 && is.na(.block_hours)){
    .block_hours = seq(4,12,2)
    var <- substitute(block_hours)
  }
  if(length(.rvu_target) == 1 && is.na(.rvu_target)){
    .rvu_target = seq(3000,12000,1000)
    var <- substitute(rvu_target)
  }
  
  # 2. import the data to be pooled. Cannot use future_pmap since 'row.prepare'
  #    is being called within future_pmap (and hence is already limited to one processor)
  specs <- expand.grid(.specialty = .specialty,
                       .cluster_size = .cluster_size,
                       .casemix = .casemix,
                       .rvu_target = .rvu_target,
                       .turnover_time = .turnover_time,
                       .block_hours = .block_hours,
                       stringsAsFactors = FALSE)
  data <- pmap(specs, import) %>%
    bind_rows() %>%
    group_by(specialty, cluster_size, casemix,
             turnover_time, block_hours, rvu_target, trial) %>%
    mutate(turnover.time.block = (n.cases-1)*turnover_time) %>%
    #define metrics to be compared using ANOVA
    summarise(n.blocks = max(block),
              time.turnover.net = sum(turnover.time.block),
              time.actual.net = sum(time.or),
              optime.net = time.actual.net - time.turnover.net)
  
  environment(row.prepare.anovas.construct_anon_row) <- environment()
  
  output <- list(
    n.blocks = data %>% row.prepare.anovas.construct_anon_row(
      .metric = substitute(n.blocks),.var = var),
    time.turnover = data %>% row.prepare.anovas.construct_anon_row(
      .metric = substitute(time.turnover.net),.var = var),
    time.actual = data %>% row.prepare.anovas.construct_anon_row(
      .metric = substitute(time.actual.net),.var = var),
    optime = data %>% row.prepare.anovas.construct_anon_row(
      .metric = substitute(optime.net),.var = var))
  return(output)
}

####  SUBR: ANOVA TABLE HELPER (ROW CREATOR)  ####

# helper function for row.prepare.anovas to do heavy lifting for row creation

row.prepare.anovas.construct_anon_row <- function(.data, .metric, .var){
  #.metric = n.block, time.turnover.net, time.actual.net, optime.net
  #.var = specialty, casemix, rvu_target, etc.
  
  result.aov <- aov(eval(.metric) ~ as.factor(eval(.var)), data = .data)
  result.tukeyHSD <- TukeyHSD(result.aov)
  p.val.aov <- summary(result.aov)[[1]][1,5]
  
  specialty <- .data$specialty %>% unique() 
  cluster_size <- .data$cluster_size %>% unique()   
  casemix <- .data$casemix %>% unique()
  turnover_time <- .data$turnover_time %>% unique()
  block_hours <- .data$block_hours %>% unique()
  rvu_target <- .data$rvu_target %>% unique()
  
  #if list OR null, need to convert to single value in output table
  if(length(specialty) != 1) .specialty <- NA   
  if(length(cluster_size) != 1) .cluster_size <- NA
  if(length(casemix) != 1) .casemix <- NA  
  if(length(block_hours) != 1) .block_hours <- NA
  if(length(turnover_time) != 1) .turnover_time <- NA
  if(length(rvu_target) != 1) .rvu_target <- NA
  
  output <- tribble(
    ~specialty, ~casemix, ~cluster_size, ~block_hours,
    ~turnover_time, ~rvu_target, ~p_value.anova, 
    .specialty, .casemix, .cluster_size, .block_hours,
    .turnover_time, .rvu_target, p.val.aov,
  )
  
  #add Tukey HSDs programmatically (because so many combinations of specialties (55))
  combos <- combn(eval(.var),2)
  for(i in 1:dim(combos)[2]){
    #create new column name spelling out combo with specialty names
    colname <- paste(combos[2,i],combos[1,i],sep="-")
    #get p value from Tukey HSD result
    p.val <- result.tukeyHSD[[1]][colname,"p adj"]
    #create new column with p.value of Tukey HSD
    output %<>% add_column(!!(paste0("p_value.",colname)) := p.val) 
  }
  
  return(output)
}

####    SUBR: SAVE ANOVA RESULTS    ####

# write anova results to .csv files on HPC

save_anovas <- function(output_raw, byvar){
  
  result <- output_raw %>% unlist(recursive = FALSE)
  
  n.blocks <- result %>% .[seq(1,length(.),4)] %>% bind_rows()
  turnover <- result %>% .[seq(2,length(.),4)] %>% bind_rows()
  time.actual <- result %>% .[seq(3,length(.),4)] %>% bind_rows()
  optime <- result %>% .[seq(4,length(.),4)] %>% bind_rows()
  
  write.csv(n.blocks,
            file =paste0(root_result,"anova_blockreqs_by",byvar,".csv"),
            row.names = FALSE)
  write.csv(turnover, 
            file = paste0(root_result,"anova_turnover-net_by",byvar,".csv"),
            row.names = FALSE)
  write.csv(time.actual,
            file = paste0(root_result,"anova_time-actual_by",byvar,".csv"),
            row.names = FALSE)
  write.csv(optime,
            file = paste0(root_result,"anova_optime_by",byvar,".csv"),
            row.names = FALSE)

  return(TRUE)
}

####  EXECUTE   ####

workers.total <- 8

print(Sys.time())
print("starting - anovas by casemix")

plan(multiprocess, workers = workers.total)

calc.specs <- expand.grid(.block_hours = seq(4,12,2),
                          .specialty = speclist,
                          .cluster_size = c(5,10,20,Sched.ANNUAL),
                          .casemix = casemixlist,
                          .turnover_time = seq(from=0,to=90,by=10),
                          .rvu_target = seq(3000,12000,1000),
                          stringsAsFactors = FALSE)
save_anovas(
  future_pmap(calc.specs %>% mutate(.casemix = NA) %>% distinct(),
              row.prepare.anovas, .progress = TRUE),
  byvar = "casemix"
)

print("done - anovas by casemix")
print(Sys.time())
print("starting - anovas by specialty")

plan(multiprocess, workers = workers.total)

save_anovas(
  future_pmap(calc.specs %>% mutate(.specialty = NA) %>% distinct(),
              row.prepare.anovas, .progress = TRUE),
  byvar = "specialty"
)

print("done - anovas by specialty")
print(Sys.time())
print("starting - anovas by turnover time")

plan(multiprocess, workers = workers.total)

save_anovas(
  future_pmap(calc.specs %>% mutate(.turnover_time = NA) %>% distinct(),
              row.prepare.anovas, .progress = TRUE),
  byvar = "turnover"
)

print("done - anovas by turnover time")
print(Sys.time())
print("starting - anova by schedule")

plan(multiprocess, workers = workers.total)

save_anovas(
  future_pmap(calc.specs %>% mutate(.cluster_size = NA) %>% distinct(),
              row.prepare.anovas, .progress = TRUE),
  byvar = "schedule"
)

print("done - anovas by schedule")
print(Sys.time())
print("starting - anovas by blocksize")

plan(multiprocess, workers = workers.total)

save_anovas(
  future_pmap(calc.specs %>% mutate(.block_hours = NA) %>% distinct(),
              row.prepare.anovas, .progress = TRUE),
  byvar = "blocksize"
)

print("done - anovas by blocksize")
print(Sys.time())
print("starting - anovas by RVU Target")

plan(multiprocess, workers = workers.total)

save_anovas(
  future_pmap(calc.specs %>% mutate(.rvu_target = NA) %>% distinct(),
              row.prepare.anovas, .progress = TRUE),
  byvar = "rvutarget"
)

print("done - anovas by blocksize")
print("done")
print(Sys.time())

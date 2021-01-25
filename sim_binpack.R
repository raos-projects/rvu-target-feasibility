#!/usr/bin/env Rscript

##################################################
## Project: RVU Target Feasibility Study
## Script purpose: Generate annual block schedules
## Date: 12/12/2020
## Author: Saieesh Rao
##################################################

args = commandArgs(trailingOnly = TRUE)

specialty.arg = strtoi(args[1])
BLOCK.MINS = strtoi(args[2])
TURNOVER_TIME = strtoi(args[3])

on.hpc = TRUE

if(specialty.arg == 4) stop("escape IR simulations")

#### LIBRARIES ####
library(gbp)
library(data.table)
library(haven)
library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)			#required for drop_na
library(future)
library(furrr)

####    CONSTANTS    ####

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

Casemix.LOW = 1
Casemix.MEDIUM = 2
Casemix.HIGH = 3
Casemix.ALL = 0

##  SPECIALTY   ##

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

speclist <- c(cardiac, general, gyne, ir, neuro, ortho,
              ent, plastics, thoracic, urology, vascular)

##  SCHEDULE  ##

#no longer using scheduling patterns; instead, cluster of
#cases to be grouped at one time (5,10,20,30,99999) by the
#bin-packing algorithm. Sched.ANNUAL is an arbitrarily large
#constant meant to encompass all cases in an annual caselog.

Sched.ANNUAL = 99999

##  SIM OPTIONS ##

GEN_TRIALS <- TRUE          #generate new data (as opposed to reading saved data)
SAVE_TRIALS <- TRUE         #write data from each trial (will overwrite existing data)
SEARCH_PREV_TRIALS <- TRUE  #if prior data, then skip sim (even if GEN_TRIALS = TRUE)

####    SUBR: REMOVE COMPLETED SIMS FROM SIM SPEC GRID    ####

#saves time because removed sims are not scheduled on processors,
#which would cancel the sim after finding completed sim data anyways.

filter.simspecs <- function(.sim.specs, .directory, .prefix) {
  
  # algorithmically generate filepath based on sim parameters, for use in next step
  filepath <- function(x.surgspec, x.rvu_goal, x.casemix,
                       x.turnover_time,x.cluster_size, x.directory, x.prefix) {
    
    path = paste0("/scratch/t.sur.rsaieesh",x.directory,
                  x.prefix,"_",x.surgspec,
                  "_RVUs-",x.rvu_goal,
                  "_blocksize-",BLOCK.MINS/60,
                  "_casemix-",x.casemix,
                  "_turnover-",x.turnover_time,
                  "-cluster-",x.cluster_size,    #typo hypen instead of underscore -
                  ".csv")                        #kept as is
    
    return(path)
  }
  
  # keep only those sim parameters for which simulation data exists inthe scratch folder
  .sim.specs %<>% filter(!file.exists(filepath(x.surgspec = .surgspec,
                                               x.rvu_goal = .rvu_goal,
                                               x.casemix = .casemix,
                                               x.turnover_time = .turnover_time,
                                               x.cluster_size = .cluster_size,
                                               x.directory = .directory,
                                               x.prefix = .prefix)))
  return(.sim.specs)
}

##  SIM SPECS ##

#create sim.specs grid and filter already completed sims.
#.block_mins and .turnover_time specified as command line arguments.
#expand.grid(...) creates a table of all combinations of arguments.

sim.specs <-
  expand.grid(.rvu_goal = seq(from=3000,to=12000,by=1000),
              .iterations = 1000,
              .surgspec = speclist[specialty.arg],
              .casemix = c(Casemix.ALL,Casemix.LOW,Casemix.MEDIUM,Casemix.HIGH),
              .block_mins = BLOCK.MINS,
              .turnover_time = TURNOVER_TIME,
              .cluster_size = c(5,10,20),
              stringsAsFactors = FALSE) %>%
  filter.simspecs("/binpack/","bp")

sim.specs.global_optimum <-
  expand.grid(.rvu_goal = seq(from=3000,to=12000,by=1000),
              .iterations = 100,
              .surgspec = speclist[specialty.arg],
              .casemix = c(Casemix.ALL,Casemix.LOW,Casemix.MEDIUM,Casemix.HIGH),
              .block_mins = BLOCK.MINS,
              .turnover_time = TURNOVER_TIME,
              .cluster_size = Sched.ANNUAL, #arbitrarily large number 99999
              stringsAsFactors = FALSE) %>%
  filter.simspecs("/binpack-global/","bp")


####  SHORT-CIRCUIT  ####

#exit program if all sims have been completed

if(dim(sim.specs)[1] == 0 && dim(sim.specs.global_optimum)[1] == 0) {
  print("Found data for all sims in scratch folder, exiting now")
  print(Sys.time())
  stop("Found data for all sims in scratch folder, exiting now")
} else {
  print(paste(dim(sim.specs)[1],"jobs remaining"))
  print(paste(dim(sim.specs.global_optimum)[1],"jobs remaining (global binpack)"))
}


####    READ DATA    ####

##  ON CLUSTER (GARDNER HPC)
print("Reading rvu brief.dta")
data <- read_dta("rvu brief.dta")

####    CREATE DIRECTORIES IF NONEXISTENT   ####

if(!dir.exists(paste(getwd(),"/binpack",sep=""))) {
  dir.create(paste(getwd(),"/binpack",sep=""))
}
if(!dir.exists(paste(getwd(),"/binpack-global",sep=""))) {
  dir.create(paste(getwd(),"/binpack-global",sep=""))
}


####    SUBR: ASSIGN CASEMIX COMPLEXITY TERTILES    ####

assign.tertiles <- function(.cases) {
  .cases %<>%
    #group by both specialty and prncptx at same time, so surgeons in same specialty
    #are compared against each other for procedures performed by multiple specialties
    group_by_at(vars(surgspec,prncptx)) %>% 
    mutate(tertile = ntile(optime,3)) %>%
    ungroup()
  return(.cases)
}


####    SUBR: SINGLE SIMULATION - PER ACTUAL OPTIME - BIN PACK    ####

sim.act.bp <- function(.rvu_goal, .iterations, .surgspec,
                       .casemix, .block_mins, .turnover_time, .cluster_size) {
  
  ## set seed for random number generation; primes to reduce seed overlap
  SEED_i <- SEED + .rvu_goal*7 + .iterations*11 + .casemix*13 +
    nchar(.surgspec)*17 + .block_mins*19 + .turnover_time*23 + .cluster_size * 29
  
  set.seed(SEED_i)
  
  #generate filename for sim and check if saved data already exists
  
  if(.cluster_size == Sched.ANNUAL){
    middle_path = 'binpack-global/bp_'
  }
  else{
    middle_path = 'binpack/bp_'
  }
  
  if(on.hpc){
    filename <- paste0(
      "/scratch/t.sur.rsaieesh/",middle_path,.surgspec,
      "_RVUs-",.rvu_goal,
      # "_sched-",.schedule,
      "_blocksize-",.block_mins/60,
      "_casemix-",.casemix,
      "_turnover-",.turnover_time,
      "-cluster-",.cluster_size,
      ".csv")
    print(paste(filename, file.exists(filename)))
  } else {
    filename <- paste0(
      "X:/",middle_path,.surgspec,
      "_RVUs-",.rvu_goal,
      # "_sched-",.schedule,
      "_blocksize-",.block_mins/60,
      "_casemix-",.casemix,
      "_turnover-",.turnover_time,
      "-cluster-",.cluster_size,
      ".csv")
    print(paste(filename, file.exists(filename)))
  }

  
  #only generate new trial data if global env variable GEN_TRIALS = TRUE
  #and there is no previous data if SEARCH_PREV_TRIALS = TRUE
  if(GEN_TRIALS && (!SEARCH_PREV_TRIALS || !file.exists(filename))){
    
    #filter only cases performed by specific subspecialty
    .cases <- data %>% filter(surgspec == .surgspec)
    
    #filter by surgeon casemix
    if(.casemix > 0){
      .cases %<>% filter(tertile == .casemix)
    }
    
    #initialize some vars for sim
    trial <- 1                      #nth trial out of 'iterations' trials
    ctr <- 1                        #row index for trials.data (never resets)    
    
    #initialize table that will record every block utilization in every trial.
    #1500*.iterations is an excessive guesstimate of table size; trimmed afterwards
    trials.data <- setDT(as.data.frame(matrix(NA, nrow = 1500*.iterations, ncol = 6)))
    trials.data <- trials.data[,lapply(.SD, as.numeric)]
    
    #set appropriate column names for binpacking analogy
    names(trials.data) <- c("oid","sku","l","d","h","w")
    #arbitrarily small number to preclude packing in depth and height dimensions.
    #only need 1D packing along length; logically depth and height would be zero
    trials.data[,'d' := 1e-5]
    trials.data[,'h' := 1e-5]
    
    
    #begin simulation - randomly sample set of cases with total RVUs meeting RVU target
    while(trial <= .iterations){
      
      #reset values for each trial
      rvus <- 0
      case_num <- 0
      
      while(rvus < .rvu_goal){
        
        case.index <- runif(1,1,dim(.cases)[1])
        
        #combine trial (ten-thousands place) and cluster (single/tens/hundreds digit).
        #necessary for packing algorithm to identify clusters to be packed (uses oid)
        trials.data[ctr, "oid" := trial*10000 + (case_num %/% .cluster_size)]
        #use sku to store cluster number (redundant for convenience)
        trials.data[ctr, "sku" := (case_num %/% .cluster_size)]
        #case length plus turnover time, packed together as a unit
        trials.data[ctr, "l" := .cases$optime[case.index] + .turnover_time]
        #weight corresponds to associated RVUs
        trials.data[ctr, "w" := .cases$totrvu[case.index]]
        
        #store case
        ctr <- ctr + 1
        rvus <- rvus + .cases$totrvu[case.index]
        case_num <- case_num + 1
      }
      
      trial <- trial+1
    }
    
    #remove empty rows created during initialization of dataframe;
    #oid picked arbitrarily since always specified
    trials.data %<>% na.omit(cols="oid")
    
    #specify dimension of the bin; length is logically block time,
    #and width/height are arbitrarily small (logically zero but
    #bin-pack algorithm requires a non-zero value for dimensions)
    bin <- data.table(
      id = c(1),
      #add turnover time to accommodate the extra turnover after last case;
      #extra turnover will be subtracted later
      l = c(.block_mins)+.turnover_time,
      d = c(1e-5),
      h = c(1e-5),
      w = c(99999) #arbitrarily high number --> no rvu max per block
    )
    
    #extract block assignment of each case over course of year;
    #convert bin-packing analogy to block packing by renaming columns
    pack_result <- (bpp_solver(trials.data, bin))$it %>%
      rename(trial = oid) %>%
      mutate(cluster = as.integer(sku)) %>%
      rename(block = tid) %>%
      rename(time.or = l) %>%
      rename(rvus = w) %>%
      select(c(trial,cluster,block,time.or,rvus))
    
    #Determine how many blocks were packed. Using the bpp_solver, only cases that fit
    #entirely within the block were packed. Cases that extended over block were not
    #packed and instead assigned a block # of zero. These cases need to be manually
    #assigned to their own unique block (not all zero)
    
    ## 1. Need to first determine the # of packed blocks (e.g., 49) and then replace
    ##    every block == 0 with a consecutive number (e.g., 50, 51, 52, etc.)
    blocks_packed <- pack_result %>%
      group_by(trial,cluster) %>%
      summarise(n = n(),
                n.pack = max(block),
                n.zero = sum(block == 0)) %>%
      ungroup()
    
    ## 2. highest block id in each trial, repeated for each block == 0 in each trial.
    ##    will use in next step
    base_ids <- rep(blocks_packed$n.pack, blocks_packed$n.zero)
    
    ## 3. generate new block ids using 'ave' with 'seq_along' function;
    
    # label zeros in consecutive sequence (1,2,3,...n)
    replacement_block_ids <- ave(pack_result[pack_result$block == 0,]$block,
                                 c(pack_result[pack_result$block == 0,]$trial,
                                   pack_result[pack_result$block == 0,]$cluster),
                                 FUN = seq_along) 
    
    # add base_ids count, to continue count without overlapping with
    # previous blocks ids (1,2,3.. --> 50, 51, 52...)
    replacement_block_ids <- replacement_block_ids + base_ids 
    
    ## 4. Next replace all the zeros with values in the sequence generated
    ##    in step 3, beginning with "n.pack+1"
    
    # flag all unscheduled cases/blocks (same thing since 1 case per block)
    # with NA for convenience; again, these are cases that were longer than block
    pack_result_cleaned <- pack_result %>%
      mutate(block = na_if(block,0))    

    # replace NAs with replacement block ids
    pack_result_cleaned[is.na(pack_result_cleaned$block),]$block <- replacement_block_ids 

    #generate summary table condensing individual cases in the same block.
    #each table row becomes one block
    pack_result_summary <- pack_result_cleaned %>%
      ungroup() %>%
      mutate(trial = floor(trial/10000)) %>%
      group_by(trial,cluster,block) %>%
      summarise(n.cases = n(),
                time.or = sum(time.or)-.turnover_time,
                rvus = sum(rvus))
    
    #renumber blocks according to trial alone, so numbering is continuous
    #across clusters (doesn't restart within each cluster).
    #the actual number of the block is not important
    pack_result_summary$block <- ave(pack_result_summary$block,
                                     pack_result_summary$trial, FUN = seq_along)
    
    #save results. finish.
    if(SAVE_TRIALS){
      write.csv(pack_result_summary,file = filename,
                row.names = FALSE)
      print(paste0("Saved ",filename))
    }
  }
}


####    SUBR: AGGREGATE SIMULATION - SCHED PER ACTUAL OPTIME - BIN PACK   ####

sim.agg.act.bp <- function(.gen_trials, .save_trials, .search_prev_trials, .sim.specs) {

  GEN_TRIALS <- .gen_trials
  SAVE_TRIALS <- .save_trials
  SEARCH_PREV_TRIALS <- .search_prev_trials

  invisible(future_pmap(.sim.specs, sim.act.bp, .progress = TRUE))
}

####    PRELIM DATA CLEANING   ####

data %<>%
  filter(surgspec == speclist[specialty.arg]) %>%
  filter(optime > 0) %>%     #remove negative / zero optimes
  filter(totrvu > 0) %>%     #remove negative / zero RVUs
  filter(workrvu > 0) %>%    #wRVU for prncptx > 0
  filter(!(otherproc1 != "NULL" & otherwrvu1 == 0)) %>%  #remove cases with procedures that
  filter(!(otherproc2 != "NULL" & otherwrvu2 == 0)) %>%  #happened but have zero RVUs.
  filter(!(otherproc3 != "NULL" & otherwrvu3 == 0)) %>% 
  filter(!(otherproc4 != "NULL" & otherwrvu4 == 0)) %>%
  filter(!(otherproc5 != "NULL" & otherwrvu5 == 0)) %>%
  filter(!(otherproc6 != "NULL" & otherwrvu6 == 0)) %>%
  filter(!(otherproc7 != "NULL" & otherwrvu7 == 0)) %>%
  filter(!(otherproc8 != "NULL" & otherwrvu8 == 0)) %>%
  filter(!(otherproc9 != "NULL" & otherwrvu9 == 0)) %>%
  filter(!(otherproc10 != "NULL" & otherwrvu10 == 0)) %>%
  select(c(prncptx, cpt, surgspec, optime, totrvu)) %>%
  assign.tertiles()


####    EXECUTE AND SAVE RESULTS   #####


GEN_TRIALS <- TRUE
SAVE_TRIALS <- TRUE
SEARCH_PREV_TRIALS <- TRUE

plan(multiprocess, workers = 4)             #4 processors per job; each job has 120 tasks
options(future.globals.maxSize = 60*1024^3) #arbitrary 60GB
future_options(seed = TRUE)                 #RNG safe for parallel computing

print(paste("sim.agg.act.bp",Sys.time()))

#cluster packing in parallel
sim.agg.act.bp(GEN_TRIALS, SAVE_TRIALS, SEARCH_PREV_TRIALS, sim.specs)

#global packing benchmark (annual caseload single-step packing) in parallel
sim.agg.act.bp(GEN_TRIALS, SAVE_TRIALS, SEARCH_PREV_TRIALS, sim.specs.global_optimum)

print(paste("done",Sys.time())) 
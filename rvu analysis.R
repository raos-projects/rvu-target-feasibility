####    LIBRARIES   ####

#must execute these lines first
library(haven)
library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)

####    READ DATA    ####
data <- read_dta(choose.files())   #select "rvu brief.dta" file containing NSQIP data (this will take a while to import)
setwd(choose.dir())                #select working directory in which "results" folder will be created
if(!dir.exists(paste(getwd(),"/results",sep=""))) {
  dir.create(paste(getwd(),"/results",sep=""))
}
if(!dir.exists(paste(getwd(),"/trials",sep=""))) {
  dir.create(paste(getwd(),"/trials",sep=""))
}

####    CONSTANTS    ####

##    RVU GOALS   ##

RVU_Goal <- 5000

##  SPEED  ##

Speed.FAST = 1 #"Fast"
Speed.MEDIUM = 2 #"Medium"
Speed.SLOW = 3 #"Slow"
Speed.ALL = 0 #"All"

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

##  RVU GOAL - SPECIALTY TABLE  ##

rvu_goal.byspecialty <-
  tribble(
    ~specialty, ~rvu_goal,
    cardiac, 5000,
    general, 5000,
    gyne, 5000,
    ir, 5000,
    neuro, 5000,
    ortho, 5000,
    ent, 5000,
    plastics, 5000,
    thoracic, 5000,
    urology, 5000,
    vascular, 5000
  )

##  SCHEDULE  ##

Sched.NEVER_OVER = "Never overtime"
Sched.SOFTCAP = "No new cases after 6hrs"
Sched.ALWAYS_OVER = "Always overtime"

BLOCK.MINS <- 480 #8hrs

SOFTCAP.MINS = 360 #<- function(){return(BLOCK.MINS - 120)};
TURNOVER_TIME = #46, the median difference between anesthesia time and surgery time when both exist
  median(
    data[!is.na(data$anetime) & !is.na(data$optime) & data$anetime > 0,]$anetime-
      data[!is.na(data$anetime) & !is.na(data$optime) & data$anetime > 0,]$optime
  );

##  SIM OPTIONS ##

GEN_TRIALS <- TRUE    #generate new data (as opposed to reading saved data)
SAVE_TRIALS <- TRUE   #write data from each trial (will overwrite existing trial data files)


####    SUBR: ASSIGN TERTILES    ####
assign.tertiles <- function(cases) {
  cases %<>%
    group_by(prncptx) %>%
    mutate(tertile = ntile(optime,3))
  return(cases)
}


####    SUBR: SINGLE SIMIULATION    ####
sim <- function(cases, rvu_goal, iterations, surgspec, schedule, speed, has_turnover_time) {
  
  ## setup for various simulation modes
  
  #add turnover time or not
  if(has_turnover_time) {
    turnover_time = TURNOVER_TIME
  }
  else{
    turnover_time = 0;
  }
  
  #only generate new trial data if global env variable GEN_TRIALS = TRUE
  if(GEN_TRIALS){

    trial <- 1                #nth trial out of 'iterations' trials
    time.max <- BLOCK.MINS    #minutes
    ctr <- 1                  #index for time.or.allblocks (never resets)    
    
    trials.data <- as.data.frame(matrix(NA, nrow = 150*iterations, ncol = 3))   #initialize table that will record every block utilization in every trial. 150*iterations is a guesstimate of table size
    names(trials.data) <- c("trial","block","time.or")                          #set appropriate column names
    
    #filter by surgeon speed
    if(speed > 0){
      cases %<>% filter(tertile == speed)
    }
    
    #condition for ending each OR block varies by schedule
    if(schedule == Sched.NEVER_OVER){
      meets.condition <- function(){
        return(time.or + cases$optime[case.index] > time.max)
      }
    }
    else if(schedule ==  Sched.ALWAYS_OVER) {
      meets.condition <- function(){
        return(time.or > time.max)
      }
    }
    else if(schedule == Sched.SOFTCAP) {
      meets.condition <- function(){
        return(time.or > SOFTCAP.MINS)
      }
    }
    else {
      meets.condition <- function(){
        stop("Scheduling condition is not properly specified")
      }
    }
    
    while(trial <= iterations){
      
      blocks <- 1
      time.or <- 0
      rvus <- 0
      while(rvus < rvu_goal) {
        
        case.index <- runif(1,1,dim(cases)[1])
        if(meets.condition()) {
          trials.data[ctr,] <- c(trial,blocks,time.or)
          blocks <- blocks+1
          ctr <- ctr + 1
          time.or <- 0
        }
        time.or <- time.or + cases$optime[case.index] + turnover_time
        rvus <- rvus + cases$totrvu[case.index]
      }
      
      trial <- trial+1
      
    }
    
    trials.data %<>% drop_na()  #remove empty rows created during initialization of dataframe
    
  }
  else {  #GEN_TRIALS == FALSE, so read from file. Will throw error if file doesn't exist
    
    trials.data <- read.csv(file = paste(getwd(),
                                         "/trials/trials_",surgspec,
                                         "_RVUs-",rvu_goal,
                                         "_sched-",schedule,
                                         "_blocksize-",BLOCK.MINS/60,
                                         "_speed-",speed,
                                         "_turnover-",turnover_time,
                                         ".csv",sep = ""),
                            header = TRUE)
    
  }
  
  if(SAVE_TRIALS && GEN_TRIALS){  #only need to save data if it was just generated ;)
    write.csv(trials.data,file = paste(getwd(),
                                       "/trials/trials_",surgspec,
                                       "_RVUs-",rvu_goal,
                                       "_sched-",schedule,
                                       "_blocksize-",BLOCK.MINS/60,
                                       "_speed-",speed,
                                       "_turnover-",turnover_time,
                                       ".csv",sep = ""),
              row.names = FALSE)
  }
  
  ##  ANALSYIS OF TRIALS.DATA
  trials.data %<>% group_by(trial)
  sumstats <- trials.data %>%
    summarise(blocks = max(block), time.or = mean(time.or))

  return(
    tribble(
      ~specialty, ~scheduling, ~speed, ~block.size, ~rvu.goal, ~blocks.mean, ~blocks.sd, ~time.or.mean, ~time.or.sd, ~turnover.time,
      surgspec, schedule, speed, BLOCK.MINS, rvu_goal, mean(sumstats$blocks), sd(sumstats$blocks), mean(trials.data$time.or), sd(trials.data$time.or), turnover_time
    )
  )
  
}

####    SUBR: SIMULATION BY SPECIALTY   ####
sim.byspecialty <- function(cases.x, iterations.x, surgspec.x) {
  time.start <- Sys.time()
  cases.x %<>%
    filter(.$surgspec == surgspec.x) %>%
    assign.tertiles()
  rvu_goal.x <- rvu_goal.byspecialty %>%
    filter(.$specialty == surgspec.x) %>%
    select(rvu_goal)
  results <- cases.x %>%
  {
    rbind(
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.ALL, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.ALL, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.ALL, has_turnover_time = FALSE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.FAST, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.FAST, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.FAST, has_turnover_time = FALSE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.MEDIUM, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.MEDIUM, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.MEDIUM, has_turnover_time = FALSE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.SLOW, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.SLOW, has_turnover_time = FALSE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.SLOW, has_turnover_time = FALSE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.ALL, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.ALL, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.ALL, has_turnover_time = TRUE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.FAST, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.FAST, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.FAST, has_turnover_time = TRUE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.MEDIUM, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.MEDIUM, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.MEDIUM, has_turnover_time = TRUE),
      sim(cases.x,rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.NEVER_OVER, speed = Speed.SLOW, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.SOFTCAP, speed = Speed.SLOW, has_turnover_time = TRUE),
      sim(cases.x, rvu_goal.x, iterations.x, surgspec.x, schedule = Sched.ALWAYS_OVER, speed = Speed.SLOW, has_turnover_time = TRUE)
    )
  }
  print(paste("Time Elapsed: ",(Sys.time()-time.start)))
  return(results)
}

####    SUBR: AGGREGATE SIMULATION    ####
sim.all <- function(cases,iterations) {
  time.start <- Sys.time()
  iter = iterations;
  results <- rbind(
    sim.byspecialty(cases,iter,cardiac),
    sim.byspecialty(cases,iter,general),
    sim.byspecialty(cases,iter,gyne),
    sim.byspecialty(cases,iter,ir),
    sim.byspecialty(cases,iter,neuro),
    sim.byspecialty(cases,iter,ortho),
    sim.byspecialty(cases,iter,ent),
    sim.byspecialty(cases,iter,plastics),
    sim.byspecialty(cases,iter,thoracic),
    sim.byspecialty(cases,iter,urology),
    sim.byspecialty(cases,iter,vascular)
  )
  print(paste("Time Elapsed for Total Simulation:",Sys.time()-time.start))
  return(results)
}

####    SUBR: META-AGGREGATE SIMULATION   ####

sim.run <- function() {
  ## vary block size, RVU targets
  
  rvu_goal.x <- 5000
  block_sizes <- c(3,4,6,8)#hours
  n.surgeons <- 10000
  
  while(rvu_goal.x <= 12000){
    for(size.block in block_sizes){
      BLOCK.MINS <- size.block*60 #convert hours to mins
      results <- sim.all(data,n.surgeons)
      write.csv(results, file = paste(getwd(),"/results/results_blocksize",size.block,"hrs_RVUtarget",rvu_goal.x,"_",".csv",sep=""),row.names = F)
    }
    rvu_goal.x <- rvu_goal.x + 500
  }
}

####    PRELIM DATA CLEANING   ####
data <- data %>%
  filter(.$optime > 0) %>%    #remove negative / zero optimes
  #filter(.$optime < 480) %>%  #remove extremely long cases
  filter(.$totrvu > 0) %>%    #remove negative / zero RVUs
  filter(!grepl('UNL',.$prncptx)) %>%     #remove "unlisted" procedures
  filter(!grepl('UNL',.$otherproc1)) %>%
  filter(!grepl('UNL',.$otherproc2)) %>%
  filter(!grepl('UNL',.$otherproc3)) %>%
  filter(!grepl('UNL',.$otherproc4)) %>%
  filter(!grepl('UNL',.$otherproc5)) %>%
  filter(!grepl('UNL',.$otherproc6)) %>%
  filter(!grepl('UNL',.$otherproc7)) %>%
  filter(!grepl('UNL',.$otherproc8)) %>%
  filter(!grepl('UNL',.$otherproc9)) %>%
  filter(!grepl('UNL',.$otherproc10)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc1)) %>%   #remove specific procedure
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc2)) %>%   #that has zero RVUs (mistake?)
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc3)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc4)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc5)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc6)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc7)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc8)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc9)) %>%
  filter(!grepl('LAPS GSTRC RSTRICTIV PX LONGITUDINAL GASTRECTOMY',.$otherproc10))

####    EXECUTE AND SAVE RESULTS   #####

results <- sim.all(data,1000)
write.csv(results, file = paste(getwd(),"/results/results_aggregate_",Sys.Date(),".csv",sep=""),row.names = F)

####    PLOT RESULTS    ####
#barplot.default(results[results]$blocks.mean,names.arg = paste(results$specialty,results$speed))

####    SANDBOX     ####
#use this space only for testing/debugging, will be deleted in final version

GEN_TRIALS <- TRUE
SAVE_TRIALS <- TRUE
sim(data,10000,100,cardiac,Sched.NEVER_OVER,Speed.ALL,FALSE)

library(haven)
library(dplyr)
library(magrittr)
library(purrr)

####    READ DATA    ####
data <- read_dta(choose.files())
setwd(choose.dir())
if(!dir.exists(paste(getwd(),"/results",sep=""))) {
  dir.create(paste(getwd(),"/results",sep=""))
}

####    CONSTANTS    ####

##    RVU GOALS   ##

#source: https://www.beckershospitalreview.com/compensation-issues/2015-physician-compensation-work-rvu-by-specialty.html

# RVU_goal.cardiac <- 10072
# RVU_goal.general <- 6736
# RVU_goal.obgyn <- 6853
# RVU_goal.urology <- 7649
# RVU_goal.otolaryng <- 6903
# 
# RVU_goal.thoracic <- 10000
# RVU_goal.neuro <- 10000
# RVU_goal.ortho <- 10000
# RVU_goal.plastics <- 10000
# RVU_goal.vascular <- 10000
# 
# RVU_goal.ir <- 10000

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

# rvu_goal.byspecialty <-
#   tribble(
#     ~specialty, ~rvu_goal,
#     cardiac, RVU_goal.cardiac,
#     general, RVU_goal.general,
#     gyne, RVU_goal.obgyn,
#     ir, RVU_goal.ir,
#     neuro, RVU_goal.neuro,
#     ortho, RVU_goal.ortho,
#     ent, RVU_goal.otolaryng,
#     plastics, RVU_goal.plastics,
#     thoracic, RVU_goal.thoracic,
#     urology, RVU_goal.urology,
#     vascular, RVU_goal.vascular
#   )

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

####    SUBR: ASSIGN TERTILES    ####
assign.tertiles <- function(cases) {
  cases %<>%
    group_by(prncptx) %>%
    mutate(tertile = ntile(optime,3))
  return(cases)
}


####    SUBR: SINGLE SIMIULATION    ####
sim <- function(cases, rvu_goal, iterations, surgspec, schedule, speed, has_turnover_time) {
  
  days.results <- c()             #total work days to reach RVU target
  time.or.q24.results <- c()      #mean hours spent in OR each day (per trial)
  time.or.q24 <- c()              #total hrs spent in OR each day
  
  time.or.alldays <- rep(NA,150*iterations)        #time spent in OR everyday, for calculating SD of daily OR time
  
  trial <- 1        #nth trial out of 'iterations' trials
  time.max <- BLOCK.MINS   #minutes
  ctr <- 1          #index for time.or.alldays (never resets)
  
  ## setup for various simulation modes
  #filter by surgeon speed
  
  if(speed > 0){
    cases %<>% filter(tertile == speed)
  }
  
  #add turnover time or not
  if(has_turnover_time) {
    turnover_time = TURNOVER_TIME
  }
  else{
    turnover_time = 0;
  }
  
  #condition for ending each OR day varies by schedule
  if(schedule == Sched.NEVER_OVER){
    meets.condition <- function(){
      return(time.or + cases$optime[case.index] > time.max)
    }
  } else if(schedule ==  Sched.ALWAYS_OVER) {
    meets.condition <- function(){
      return(time.or > time.max)
    }
  } else if(schedule == Sched.SOFTCAP) {
    meets.condition <- function(){
      return(time.or > SOFTCAP.MINS)
    }
  } else {
    meets.condition <- function(){
      stop("Scheduling condition is not properly specified")
    }
  }
  
  while(trial <= iterations){
    
    days <- 1
    time.or <- 0
    rvus <- 0
    while(rvus < rvu_goal) {
      
      case.index <- runif(1,1,dim(cases)[1])
      if(meets.condition()) {
        days <- days+1
        time.or.alldays[ctr] <- time.or/60
        ctr <- ctr + 1
        time.or <- 0
      }
      time.or <- time.or + cases$optime[case.index] + turnover_time
      rvus <- rvus + cases$totrvu[case.index]
    }
    
    days.results[[trial]] <- days
    trial <- trial+1
    
  }
  
  time.or.alldays %<>% na.omit()
  
  return(
    tribble(
      ~specialty, ~scheduling, ~speed, ~days.mean, ~days.sd, ~hrs.mean, ~hrs.sd, ~turnover.time,
      surgspec, schedule, speed, mean(days.results), sd(days.results), mean(time.or.alldays), sd(time.or.alldays), turnover_time
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

## vary block size, RVU targets

rvu_goal.x <- 5000
block_sizes <- c(3,4,6,8)#hours
n.surgeons <- 10000

while(rvu_goal.x <= 12000){
  for(size.block in block_sizes){
    BLOCK.MINS <- size.block*60 #convert hours to mins
    results <- sim.all(data,n.surgeons)
    write.csv(results, file = paste(getwd(),"/results/results_blocksize",size.block,"hrs_RVUtarget",rvu_goal.x,".csv",sep=""),row.names = F)
  }
  rvu_goal.x <- rvu_goal.x + 500
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
results %>% barplot.default(results[results]$days.mean,names.arg = paste(results$specialty,results$speed))
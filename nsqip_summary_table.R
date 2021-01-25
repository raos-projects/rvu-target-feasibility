####  LIBRARIES   ####
library(dplyr)
library(magrittr)
library(tidyr)  #for drop_na
library(haven)
library(ggplot2)
library(gtools) #for na.replace
library(rJava)
library(rChoiceDialogs)

####    READ DATA    ####

# path for files on local computer (Windows)
root <- "C:/Users/Saieesh Rao/Documents/_Pritzker School of Medicine/Miscellaneous/Misc Research/Turaga RVU Project/nsqip-analysis/"

##  ON CLUSTER
#print("Reading rvu brief.dta")
#data <- read_dta("rvu brief.dta")

##  ON WINDOWS (LAPTOP)
setwd(rchoose.dir())                #select working directory in which "results" and "trials" folders will be created
data <- read_dta(choose.files())   #select "rvu brief.dta" file containing NSQIP data (this will take a while to import)

##  ON MAC  (CRERAR)
#setwd(rchoose.dir())                #select working directory in which "results" folder will be created
#data <- read_dta(file.choose())   #select "rvu brief.dta" file containing NSQIP data (this will take a while to import)
#setwd("/users/Saieesh/Desktop/NSQIP Analysis")


##  ON LINUX.CS.UCHICAGO.EDU
#print("Reading rvu brief.dta")
#data <- read_dta("/home/saieesh/Desktop/rvu brief.dta") #on linux.cs.uchicago.edu
#Must already be in /home/Saieesh/Desktop to work, cannot set home directory programmatically in cmd line

####    CONSTANTS    ####

##  SPECIALTY  ##

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

####    SUBR: TERTILES BY SURGSPEC & PRNCPTX    ####
assign.tertiles <- function(.cases) {
  .cases %<>%
    #group by both specialty and prncptx at same time, so surgeons in same specialty
    #are compared against each other for procedures performed by multiple specialties
    group_by_at(vars(surgspec,prncptx)) %>%
    mutate(tertile = ntile(optime,3)) %>%
    ungroup()
  return(.cases)
}

####    SUBR: MEDIAN OPTIMES BY SURGSPEC & PRNCPTX    ####

standard.optimes <- function(.cases, .casemix) {      #assumes global var 'data' will be assigned median optimes

  if(.casemix != Casemix.ALL) {
    .cases %<>% filter(tertile == .casemix)
  }
  
  optime.medians <- .cases %>%
    group_by(surgspec,prncptx) %>%
    summarise(n.cases = n(),
              std.optime = median(optime),
              rvus.percase.min = min(totrvu),
              rvus.percase.10th = quantile(totrvu, probs=0.10, type=2),
              rvus.percase.25th = quantile(totrvu, probs=0.25, type=2),
              rvus.percase.50th = median(totrvu),
              rvus.percase.75th = quantile(totrvu, probs=0.75, type=2),
              rvus.percase.90th = quantile(totrvu, probs=0.90, type=2),
              rvus.percase.max = max(totrvu),
              rvus.percase.mean = mean(totrvu),
              rvus.percase.sd = sd(totrvu)) %>%  #std.optime = 'standard optime'
    mutate(casemix = switch(.casemix+1,"All","Low","Medium","High"))
  return(optime.medians)
}

assign.standard.optimes <- function(.cases) {   #assume tertiles already assigned, casemix already filtered; must pass subset as '.cases'
  .cases %<>%
    #group_by(prncptx) %>% #amended from group_by(prncptx) on 1/13/2020 (just to be safe)
    group_by_at(vars(surgspec,prncptx)) %>%
    mutate(std.optime = median(optime)) %>%   #std.optime = 'standard optime'
    ungroup()
  return(.cases)
}

####  SUBR: RECORD EXCLUSION COUNTS FOR STROBE DIAGRAM   ####

# table for tracking N's after each step of cleaning
forksum.trail <- data %>% group_by(surgspec, .drop=FALSE) %>% summarise()

# function for saving N's to table (forksum.trail) after each step of cleaning
forksum <- function(.data, .colname) {

  forksum.trail <<-
    .data %>% 
    summarise(!!.colname := n()) %>%
    dplyr::select(!!.colname) %>%
    {bind_cols(forksum.trail,.)} %>%
    print()
  return(.data)
}

####  SUMMARY STATS   ####

##  1. STROBE DIAGRAM (DATA CLEANING)
# count cases by specialty during data cleaning
data2 <- data 
data2 %<>% mutate(surgspec = factor(surgspec, levels = c(cardiac, general, gyne, ir, neuro,
                                                         ortho, "Other", ent, plastics,
                                                         thoracic, "Unknown", urology, vascular)))
data2 %<>%
  group_by(surgspec, .drop=FALSE) %>%
  forksum("n.cases") %>%
  filter(optime > 0) %>%
  forksum("optime.gtz") %>%     #gtz = greater than zero
  filter(workrvu > 0) %>%
  forksum("workrvu.gtz") %>%
  filter(totrvu > 0) %>%
  forksum("totrvu.gtz") %>%
  filter(!(otherproc1 != "NULL" & otherwrvu1 == 0)) %>%  #remove cases with procedures that
  filter(!(otherproc2 != "NULL" & otherwrvu2 == 0)) %>%  #happened but have zero RVUs (unsure
  filter(!(otherproc3 != "NULL" & otherwrvu3 == 0)) %>%  #why that happens)
  filter(!(otherproc4 != "NULL" & otherwrvu4 == 0)) %>%
  filter(!(otherproc5 != "NULL" & otherwrvu5 == 0)) %>%
  filter(!(otherproc6 != "NULL" & otherwrvu6 == 0)) %>%
  filter(!(otherproc7 != "NULL" & otherwrvu7 == 0)) %>%
  filter(!(otherproc8 != "NULL" & otherwrvu8 == 0)) %>%
  filter(!(otherproc9 != "NULL" & otherwrvu9 == 0)) %>%
  filter(!(otherproc10 != "NULL" & otherwrvu10 == 0)) %>%
  forksum("rm.zerorvu") %>%
  filter(surgspec != "Unknown" && surgspec != "Other" && surgspec != ir) %>%
  forksum("specfilter")
write.csv(forksum.trail,paste0(root,"results/nsqip_data_clean.csv"), row.names = FALSE)

# code below can be used to demonstrate efficacy of original code in sim_generic.R (missing ! to select for mismatch)
# cannot use below code on 'grouped' dataframe (e.g. group_by(surgspec)), else throws error
# sim_generic.R unintentionally worked around this by filtering only cases for one surgspec (rather than grouping)
# data %>% filter((!.$otherproc1 == "NULL" & .$otherwrvu1 == 0)) %>% select(otherproc1,otherwrvu1) %>% View()

# data %>%
#   group_by(surgspec) %>%
#   filter(.$optime > 0) %>%    #remove negative / zero optimes
#   #filter(.$optime < 480) %>%  #remove extremely long cases
#   filter(.$totrvu > 0) %>%    #remove negative / zero RVUs
#   filter(.$workrvu > 0) %>%
#   filter(!(!.$otherproc1 == "NULL" & .$otherwrvu1 == 0)) %>%  #remove cases with procedures that
#   filter(!(!.$otherproc2 == "NULL" & .$otherwrvu2 == 0)) %>%  #happened but have zero RVUs (unsure
#   filter(!(!.$otherproc3 == "NULL" & .$otherwrvu3 == 0)) %>%  #why that happens)
#   filter(!(!.$otherproc4 == "NULL" & .$otherwrvu4 == 0)) %>%
#   filter(!(!.$otherproc5 == "NULL" & .$otherwrvu5 == 0)) %>%
#   filter(!(!.$otherproc6 == "NULL" & .$otherwrvu6 == 0)) %>%
#   filter(!(!.$otherproc7 == "NULL" & .$otherwrvu7 == 0)) %>%
#   filter(!(!.$otherproc8 == "NULL" & .$otherwrvu8 == 0)) %>%
#   filter(!(!.$otherproc9 == "NULL" & .$otherwrvu9 == 0)) %>%
#   filter(!(!.$otherproc10 == "NULL" & .$otherwrvu10 == 0)) %>%
#   select(c(prncptx, cpt, surgspec, optime, totrvu)) %>%
#   assign.tertiles()

# 2. ESTIMATED OPERATING TIMES BY PRNCPTX WITH RVU DISTRIBUTION
# estimated times = median times within each tertile. RVUs vary because prncptx does not account for other procedures in case
std.optimes <- 
  data2 %>% 
  assign.tertiles() %>%
  {bind_rows(list(
    standard.optimes(.,Casemix.ALL),
    standard.optimes(.,Casemix.LOW),
    standard.optimes(.,Casemix.MEDIUM),
    standard.optimes(.,Casemix.HIGH)
  ))} %>%
  pivot_wider(names_from="casemix", values_from=c("n.cases",
                                                "std.optime",
                                                "rvus.percase.min",
                                                "rvus.percase.10th",
                                                "rvus.percase.25th",
                                                "rvus.percase.50th",
                                                "rvus.percase.75th",
                                                "rvus.percase.90th",
                                                "rvus.percase.max",
                                                "rvus.percase.mean",
                                                "rvus.percase.sd")) %>%
  replace_na(list(n.cases_Medium = 0, n.cases_High = 0))
write.csv(std.optimes, file=paste0(root,"results/nsqip_standard_optimes.csv"), row.names=FALSE)

# 3. DISTRIBUTION OF OPERATING TIMES, TOTAL RVUS, AND HOURLY RVUS ACROSS SPECIALTIES, TERTILES

#View(case_stats_speed %>% filter(is.na(rvu_hourly)))

case_stats_casemix <-
  data2 %>%
  assign.tertiles() %>%
  group_by(surgspec,tertile) %>%
  mutate(rvu_hourly = totrvu/optime*60) %>%
  summarise(n.cases = n(),
            n.prncptx = n_distinct(prncptx),
            optime_min = min(optime),
            optime_10th = quantile(optime, probs = 0.10, type = 2),
            optime_25th = quantile(optime, probs = 0.25, type = 2),
            optime_50th = median(optime),
            optime_75th = quantile(optime, probs = 0.75, type = 2),
            optime_90th = quantile(optime, probs = 0.90, type = 2),
            optime_max = max(optime),
            optime_mean = mean(optime),
            optime_sd = sd(optime),
            rvu_min = min(totrvu),
            rvu_10th = quantile(totrvu, probs = 0.10, type = 2),
            rvu_25th = quantile(totrvu, probs = 0.25, type = 2),
            rvu_50th = median(totrvu),
            rvu_75th = quantile(totrvu, probs = 0.75, type = 2),
            rvu_90th = quantile(totrvu, probs = 0.90, type = 2),
            rvu_max = max(totrvu),
            rvu_mean = mean(totrvu),
            rvu_sd = sd(totrvu),
            rvu_hourly_min = min(rvu_hourly),
            rvu_hourly_10th = quantile(rvu_hourly, probs = 0.10, type = 2),
            rvu_hourly_25th = quantile(rvu_hourly, probs = 0.25, type = 2),
            rvu_hourly_50th = median(rvu_hourly),
            rvu_hourly_75th = quantile(rvu_hourly, probs = 0.75, type = 2),
            rvu_hourly_90th = quantile(rvu_hourly, probs = 0.90, type = 2),
            rvu_hourly_max = max(rvu_hourly),
            rvu_hourly_mean = mean(rvu_hourly),
            rvu_hourly_sd = sd(rvu_hourly))
case_stats_nocasemix <- 
  data2 %>%
  group_by(surgspec) %>%
  mutate(rvu_hourly = totrvu/optime*60) %>%
  summarise(n.cases = n(),
            n.prncptx = n_distinct(prncptx),
            optime_min = min(optime),
            optime_10th = quantile(optime, probs = 0.10, type = 2),
            optime_25th = quantile(optime, probs = 0.25, type = 2),
            optime_50th = median(optime),
            optime_75th = quantile(optime, probs = 0.75, type = 2),
            optime_90th = quantile(optime, probs = 0.90, type = 2),
            optime_max = max(optime),
            optime_mean = mean(optime),
            optime_sd = sd(optime),
            rvu_min = min(totrvu),
            rvu_10th = quantile(totrvu, probs = 0.10, type = 2),
            rvu_25th = quantile(totrvu, probs = 0.25, type = 2),
            rvu_50th = median(totrvu),
            rvu_75th = quantile(totrvu, probs = 0.75, type = 2),
            rvu_90th = quantile(totrvu, probs = 0.90, type = 2),
            rvu_max = max(totrvu),
            rvu_mean = mean(totrvu),
            rvu_sd = sd(totrvu),
            rvu_hourly_min = min(rvu_hourly),
            rvu_hourly_10th = quantile(rvu_hourly, probs = 0.10, type = 2),
            rvu_hourly_25th = quantile(rvu_hourly, probs = 0.25, type = 2),
            rvu_hourly_50th = median(rvu_hourly),
            rvu_hourly_75th = quantile(rvu_hourly, probs = 0.75, type = 2),
            rvu_hourly_90th = quantile(rvu_hourly, probs = 0.90, type = 2),
            rvu_hourly_max = max(rvu_hourly),
            rvu_hourly_mean = mean(rvu_hourly),
            rvu_hourly_sd = sd(rvu_hourly))
case_stats <- 
  bind_rows(case_stats_casemix, case_stats_nocasemix) %>%
  na.replace(0) %>%
  rename(casemix = tertile) %>%
  arrange(surgspec, casemix)

#t.test(data2[data2$tertile == 1,]$totrvu, data2[data2$tertile == 2,]$totrvu)
write.csv(case_stats, file=paste0(root,"results/nsqip_dataset_stats.csv"), row.names=FALSE)

# 4. PLOT SIMILARITY OF RVU MIX (AND DIFFERENCE IN HOURLY RVU MIX) ACROSS TERTILES
# histogram and CDF plots of RVU values, hourly RVU rates across tertiles

pdf(file = paste0(root,"results/NSQIP Optime-RVU Distributions.pdf"), onefile = TRUE)
for(spec in speclist) {
  #calc approx xlims 
  optime_xlim = case_stats %>% subset(surgspec == spec) %>% subset(casemix == 0) %>% pull(optime_90th) * 3
  rvu_xlim = case_stats %>% subset(surgspec == spec) %>% subset(casemix == 0) %>% pull(rvu_90th) * 3
  hourly_xlim = case_stats %>% subset(surgspec == spec) %>% subset(casemix == 0) %>% pull(rvu_hourly_90th) * 3
  
  # hist: totrvu
  print(
    data2 %>% assign.tertiles() %>% 
      filter(surgspec == general) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + geom_histogram(binwidth=2,position="identity",alpha=0.2,aes(x=totrvu,y=..density..,color = Casemix, fill = Casemix)) + xlim(0,rvu_xlim) +
      labs(title = spec, subtitle = "Accesible RVU distribution across casemixes") + xlab("RVUs per case") + ylab("Density")
  )
  #ecdf: totrvu
  print(
    data2 %>% assign.tertiles() %>%
      filter(surgspec == general) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + stat_ecdf(aes(totrvu, colour = Casemix)) + xlim(0, rvu_xlim) +
      labs(title = spec, subtitle = "Accesible RVU distribution across casemixes") + xlab("RVUs per case") + ylab("Percentile")
  )
  #hist: optime
  print(
    data2 %>% assign.tertiles() %>% 
      filter(surgspec == general) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + geom_histogram(binwidth=10,position="identity",alpha=0.2,aes(x=optime,y=..density..,color = Casemix, fill = Casemix)) + xlim(0,optime_xlim) +
      labs(title = spec, subtitle = "Case lengths across casemixes") + xlab("Operating time per case") + ylab("Density")
  )
  #ecdf: optime
  print(
    data2 %>% assign.tertiles() %>%
      filter(surgspec == general) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + stat_ecdf(aes(optime, colour = Casemix)) + xlim(0, optime_xlim) +
      labs(title = spec, subtitle = "Case lengths across casemixes") +xlab("Operating time per case") + ylab("Percentile")
  )
  #hist: rvu_hourly
  print(
    data2 %>% assign.tertiles() %>%
      filter(surgspec == general) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      mutate(rvu_hourly = totrvu/optime*60) %>%
      ggplot() + geom_histogram(binwidth=2,position="identity",alpha=0.2,aes(x=rvu_hourly,y=..density..,color = Casemix, fill = Casemix)) + xlim(0,hourly_xlim) +
      labs(title = spec, subtitle = "Hourly RVU production in cases accessible to each casemix") + xlab("Hourly RVU production during case") + ylab("Density")
  )
  #ecdf: rvu_hourly
  print(
    data2 %>% assign.tertiles() %>%
      filter(surgspec == general) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      mutate(rvu_hourly = totrvu/optime*60) %>%
      ggplot() + stat_ecdf(aes(x = rvu_hourly, colour = Casemix)) + xlim(0,hourly_xlim) +
      labs(title = spec, subtitle = "Hourly RVU production in cases accessible to each casemix") + xlab("Hourly RVU production during case") + ylab("Percentile")
  )
}
dev.off()
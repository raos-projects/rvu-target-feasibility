####  LIBRARIES   ####
library(dplyr)
library(magrittr)
library(tidyr)  #for drop_na
library(haven)
library(ggplot2)
library(gtools) #for na.replace

####    READ DATA    ####

##  ON CLUSTER
print("Reading rvu brief.dta")
data <- read_dta("rvu brief.dta")
root <- '/scratch/t.sur.rsaieesh/results-binpack/'

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

####  SUBR: RECORD EXCLUSION COUNTS FOR STROBE DIAGRAM   ####

# count cases by specialty during data cleaning
data2 <- data #duplicate data to leave original unmodified
data2 %<>% mutate(surgspec = factor(surgspec,
                                    levels = c(cardiac, general, gyne, ir, neuro,
                                               ortho, "Other", ent, plastics,
                                               thoracic, "Unknown", urology, vascular)))

# table for tracking N's after each step of cleaning
forksum.trail <- data2 %>% group_by(surgspec, .drop=FALSE) %>% summarise()

# function for saving N's to table (forksum.trail) after each step of cleaning
forksum <- function(.data, .colname) {

  forksum.trail <<-
    .data %>% 
    summarise(!!.colname := n()) %>%
    dplyr::select(!!.colname) %>%
    {bind_cols(forksum.trail,.)} #append to right side
  
  print(forksum.trail)
  
  .data %>%
    ungroup() %>%
    group_by(surgspec, .drop=FALSE)
  
  return(.data)
}

####  SUMMARY STATS   ####

##  1. STROBE DIAGRAM (DATA CLEANING)

## NOTE: The following code requires at least
## dplyr 0.8.0 for use of the .drop=FALSE argument
## in group_by(...)

# clean data while recording N's after each step
data2 %<>%
  group_by(surgspec, .drop=FALSE) %>%
  forksum("n.cases") %>%
  filter(optime > 0) %>%
  forksum("optime.gtz") %>%     #gtz = greater than zero
  filter(workrvu > 0) %>%
  forksum("workrvu.gtz") %>%
  filter(totrvu > 0) %>%
  forksum("totrvu.gtz") %>%
  filter(!(otherproc1 != "NULL" & otherwrvu1 == 0)) %>%  #remove cases with procedures
  filter(!(otherproc2 != "NULL" & otherwrvu2 == 0)) %>%  #that have zero RVUs
  filter(!(otherproc3 != "NULL" & otherwrvu3 == 0)) %>%
  filter(!(otherproc4 != "NULL" & otherwrvu4 == 0)) %>%
  filter(!(otherproc5 != "NULL" & otherwrvu5 == 0)) %>%
  filter(!(otherproc6 != "NULL" & otherwrvu6 == 0)) %>%
  filter(!(otherproc7 != "NULL" & otherwrvu7 == 0)) %>%
  filter(!(otherproc8 != "NULL" & otherwrvu8 == 0)) %>%
  filter(!(otherproc9 != "NULL" & otherwrvu9 == 0)) %>%
  filter(!(otherproc10 != "NULL" & otherwrvu10 == 0)) %>%
  forksum("rm.zerorvu") %>%
  filter(surgspec != "Unknown" && surgspec != "Other" && surgspec != ir) %>%
  group_by(surgspec, .drop=FALSE) %>%
  forksum("specfilter")

write.csv(forksum.trail,
          paste0(root,"nsqip_data_clean.csv"),
          row.names = FALSE)

# 2. DISTRIBUTION OF OPERATING TIMES, TOTAL RVUS, AND
#    HOURLY RVUS ACROSS SPECIALTIES, TERTILES, PRNCPTX

# 2A. Generate table describing stats for each prncptx with
#     consideration for surgspec and casemix classification
prncptx_stats_casemix <-
  data2 %>%
  assign.tertiles() %>%
  group_by(surgspec,tertile,prncptx) %>%
  mutate(rvu_hourly = totrvu/optime*60) %>%
  summarise(n.cases = n(),
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

# 2B. Generate table describing stats for each prncptx with
#     consideration for surgspec but include all casemixes
prncptx_stats_nocasemix <- 
  data2 %>%
  group_by(surgspec,prncptx) %>%
  mutate(rvu_hourly = totrvu/optime*60) %>%
  summarise(n.cases = n(),
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

# 2C. Merge 2A and 2B for combined summary table, and save
prncptx_stats <- 
  bind_rows(prncptx_stats_casemix, prncptx_stats_nocasemix) %>%
  na.replace(0) %>%
  rename(casemix = tertile) %>%
  arrange(surgspec, casemix, prncptx)

write.csv(prncptx_stats,
          file=paste0(root,"nsqip_prncptx_stats.csv"),
          row.names=FALSE)

# 3. DISTRIBUTION OF OPERATING TIMES, TOTAL RVUS, AND
#    HOURLY RVUS ACROSS SPECIALTIES, TERTILES

# 3A. Generate table describing stats for each casemix with
#     consideration for surgspec

surgeon_stats_casemix <-
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

# 3B. Generate table describing stats for each surgspec
#     independent of casemix considerations

surgeon_stats_nocasemix <- 
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

# 3C. Merge 3A and 3B for combined summary table, and save
surgeon_stats <- 
  bind_rows(surgeon_stats_casemix, surgeon_stats_nocasemix) %>%
  na.replace(0) %>%
  rename(casemix = tertile) %>%
  arrange(surgspec, casemix)

write.csv(surgeon_stats,
          file=paste0(root,"nsqip_casemix_stats.csv"),
          row.names=FALSE)

# 4. PLOT OPTIMES, RVU MIX, AND HOURLY RVU MIX ACROSS TERTILES

data2 %<>%
  assign.tertiles() %>%
  mutate(rvu_hourly = totrvu/optime*60)

# declare pdf file to which figures will be saved
pdf(file = paste0(root,"NSQIP Optime-RVU Distributions.pdf"), onefile = TRUE)

for(spec in speclist) {
  
  #calculate approx xlims 
  optime_xlim = surgeon_stats %>%
    subset(surgspec == spec) %>% subset(casemix == 0) %>% pull(optime_90th) * 3
  rvu_xlim = surgeon_stats %>%
    subset(surgspec == spec) %>% subset(casemix == 0) %>% pull(rvu_90th) * 3
  hourly_xlim = surgeon_stats %>%
    subset(surgspec == spec) %>% subset(casemix == 0) %>% pull(rvu_hourly_90th) * 3
  
  # hisogram of total wRVUs per case (totrvu)
  print(
    data2 %>%
      filter(surgspec == spec) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() +
      geom_histogram(binwidth=2,position="identity",alpha=0.2,
                     aes(x=totrvu,y=..density..,color = Casemix, fill = Casemix)) +
      xlim(0,rvu_xlim) +
      labs(title = spec, subtitle = "Accesible RVU distribution across casemixes") +
      xlab("RVUs per case") + ylab("Density")
  )
  
  #empirical CDF of total wRVUs per case (totrvu)
  print(
    data2 %>%
      filter(surgspec == spec) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + stat_ecdf(aes(totrvu, colour = Casemix)) + xlim(0, rvu_xlim) +
      labs(title = spec, subtitle = "Accesible RVU distribution across casemixes") +
      xlab("RVUs per case") + ylab("Percentile")
  )
  
  #histogram of optime per case
  print(
    data2 %>%
      filter(surgspec == spec) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + geom_histogram(binwidth=10,position="identity",alpha=0.2,
                                aes(x=optime,y=..density..,
                                    color = Casemix, fill = Casemix)) +
      xlim(0,optime_xlim) +
      labs(title = spec, subtitle = "Case lengths across casemixes") +
      xlab("Operating time per case") + ylab("Density")
  )
  
  #empirical CDF of optime per case
  print(
    data2 %>%
      filter(surgspec == spec) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + stat_ecdf(aes(optime, colour = Casemix)) +
      xlim(0, optime_xlim) +
      labs(title = spec, subtitle = "Case lengths across casemixes") +
      xlab("Operating time per case") + ylab("Percentile")
  )
  
  #histogram of hourly RVUs per case
  print(
    data2 %>%
      filter(surgspec == spec) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + geom_histogram(binwidth=2,position="identity",alpha=0.2,
                                aes(x=rvu_hourly,y=..density..,
                                    color = Casemix, fill = Casemix)) +
      xlim(0,hourly_xlim) +
      labs(title = spec,
           subtitle = "Hourly RVU production in cases accessible to each casemix") +
      xlab("Hourly RVU production during case") + ylab("Density")
  )
  
  #empirical CDF of hourly RVUs per case
  print(
    data2 %>%
      filter(surgspec == spec) %>%
      mutate(Casemix = case_when(tertile == 1 ~ "Low",
                               tertile == 2 ~ "Medium",
                               tertile == 3 ~ "High")) %>%
      ggplot() + stat_ecdf(aes(x = rvu_hourly, colour = Casemix)) +
      xlim(0,hourly_xlim) +
      labs(title = spec,
           subtitle = "Hourly RVU production in cases accessible to each casemix") +
      xlab("Hourly RVU production during case") + ylab("Percentile")
  )
}
dev.off()
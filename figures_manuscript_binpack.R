#!/usr/bin/env Rscript

##################################################
## Project: RVU Target Feasibility Study
## Script purpose: Generate Figures
## Date: 1/5/2021
## Author: Saieesh Rao
##################################################

####  LIBRARIES  ####
library(dplyr)
library(magrittr)
library(stringr)
library(future)
library(furrr)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(data.table)

####  CONSTANTS  ####

on.hpc = TRUE  #change to 'TRUE' once running on Gardner HPC

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

####  LOAD SUMMARY TABLE  ####

if(on.hpc){
  root_results = '/scratch/t.sur.rsaieesh/results-binpack/'
  root_figures = '/scratch/t.sur.rsaieesh/figures-binpack/'
} else {
  root_results = 'X:/results-binpack/'
  root_figures = 'X:/figures-binpack/'
}

alldata <- fread(paste0(root_results,"summary_table_binpack.csv"))
hrvus.util <- fread(paste0(root_results,"linear_models_hrvus_utilized.csv"))
hrvus.admin <- fread(paste0(root_results,"linear_models_hrvus_admin.csv"))

####  UTILITIES  ####

import <- function(.specialty, .rvu_goal,
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

# import each simulation's result summarized as a single row for each simulation
import.summary <- function(.specialty, .rvu_goal,
                           .block_hours, .cluster_size, .casemix, .turnover_time) {
  data <- import(.specialty, .rvu_goal,
                 .block_hours, .cluster_size, .casemix, .turnover_time)
  summ <- data %>% 
    group_by(trial, specialty_spec, cluster.size_spec,
             casemix_spec, rvu.target_spec, block.size_spec, turnover.time_spec) %>%
    summarize(n.blocks = max(block),
              optime.net = sum(time.or - UQ(.turnover_time)*(n.cases-1)),
              rvus.net = sum(rvus)) %>%
    ungroup()
  return(summ)
}

####  GLOBAL COLOR SCALES  ####

# Modify color scale definitions here to globally
# change color scheme across figures. Functions
# return color schemes from the 'ggsci' package.

# For arguments to passed to ggsci functions,
# color schemes must be specified in figure
# code manually.

scale_fill <- function() {
  return(scale_fill_jama())
}
scale_color <- function() {
  return(scale_color_jama())
}

####  PAPER 1 - FIGURE 1 - Net Optime versus RVU Target by Specialty Cluster  ####

# import relevant data for 3 panel figure using parallel processing.
# gNdata nomenclature signifies grob1, grob2, grob3 corresponding to each of the panels
plan(multiprocess, workers = 3)

g1data <- future_pmap(expand.grid(.specialty = c(cardiac, vascular, neuro, ortho),
                                  .cluster_size = 10,
                                  .casemix = casemixlist,
                                  .block_hours = 8,
                                  .turnover_time = 60,
                                  .rvu_goal = seq(3000,12000,1000)),
                      import.summary, .progress = TRUE) %>% bind_rows()
g2data <- future_pmap(expand.grid(.specialty = c(general, gyne, urology),
                                  .cluster_size = 10,
                                  .casemix = casemixlist,
                                  .block_hours = 8,
                                  .turnover_time = 60,
                                  .rvu_goal = seq(3000,12000,1000)),
                      import.summary, .progress = TRUE) %>% bind_rows()
g3data <- future_pmap(expand.grid(.specialty = c(ent, plastics, thoracic),
                                  .cluster_size = 10,
                                  .casemix = casemixlist,
                                  .block_hours = 8,
                                  .turnover_time = 60,
                                  .rvu_goal = seq(3000,12000,1000)),
                      import.summary, .progress = TRUE) %>% bind_rows()

# format data for ggplot usage
gdata <- list(g1data,g2data,g3data)
gplot <- list()
names(gdata[[1]])

# create same figure for all casemixes
for(j in casemixlist) {
  # create each of the 3 panels separately and store in gplot
  for(i in 1:3) {
    print(paste("CYCLE:",i,j))
    params <- gdata[[i]] %>%
      filter(casemix_spec == j) %>%
      group_by(specialty_spec, rvu.target_spec, casemix_spec) %>%
      summarise(sd = sd(optime.net/60),
                mean = mean(optime.net/60),
                ymin = mean - 1.96*sd,
                ymax = mean + 1.96*sd) %>%
      ungroup() %>%
      mutate(specialty_spec = as.character(specialty_spec)) %>%
      # rename "Otoloaryngology (ENT)" as just "ENT"
      mutate(specialty_spec =
               str_replace_all(as.character(specialty_spec),
                               "Otolaryngology \\(ENT\\)", c("ENT")))

    # store panel in gplot (list)
    gplot[[i]] <- ggplot() + xlim(3000,12000) + ylim(50,1800) +
      geom_line(data = params,
                aes(x = rvu.target_spec, y = mean, color = specialty_spec),
                size = 1) +
      geom_ribbon(data = params,
                  aes(x = rvu.target_spec, ymin = ymin,
                      ymax = ymax, fill = specialty_spec),
                  alpha = 0.2) +
      geom_vline(xintercept = 12000) +
      geom_vline(xintercept = 3000) +
      labs(x = "RVU Target",
           y = "Net Operating Time",
           color = "Specialty",
           fill = "Specialty") +
      scale_fill() +
      scale_color() +
      annotate("text",x=3500,y=1500,label=c("A","B","C")[[i]]) +
      theme_bw() +
      theme(legend.position = c(0.435, c(0.82,0.845,0.845)[[i]]),
            legend.text=element_text(size=rel(0.7)),
            legend.title = element_text(size = rel(0.7)))
  }
  
  # save PNG file of all three panels
  png(paste0(root_figures,'paper1/p1figure1_optime_',
             c("All","Low","Medium","High")[[j+1]],'.png'),
      width = 6.8, height=5.8, units = "in", res = 300)
  # arrange all three panels in a 1x3 grid
  grid.arrange(grobs = gplot, ncol = 3)
  graphics.off()
}

####  PAPER 1 - FIGURE 2 - Block Reqs vs RVU target by specialty cluster  ####

# import relevant data for 3 panel figure using parallel processing.
# gNdata nomenclature signifies grob1, grob2, grob3 corresponding to each of the panels
plan(multiprocess, workers = 3)

g1data <- future_pmap(expand.grid(.specialty = c(cardiac, vascular, neuro, ortho),
                                  .cluster_size = 10,
                                  .casemix = casemixlist,
                                  .block_hours = 8,
                                  .turnover_time = 60,
                                  .rvu_goal = seq(3000,12000,1000)),
                      import.summary, .progress = TRUE) %>% bind_rows()
g2data <- future_pmap(expand.grid(.specialty = c(general, gyne, urology),
                                  .cluster_size = 10,
                                  .casemix = casemixlist,
                                  .block_hours = 8,
                                  .turnover_time = 60,
                                  .rvu_goal = seq(3000,12000,1000)),
                      import.summary, .progress = TRUE) %>% bind_rows()
g3data <- future_pmap(expand.grid(.specialty = c(ent, plastics, thoracic),
                                  .cluster_size = 10,
                                  .casemix = casemixlist,
                                  .block_hours = 8,
                                  .turnover_time = 60,
                                  .rvu_goal = seq(3000,12000,1000)),
                      import.summary, .progress = TRUE) %>% bind_rows()

# format data for ggplot usage
gdata <- list(g1data,g2data,g3data)
gplot <- list()
names(gdata[[1]])

# create same figure for all casemixes
for(j in casemixlist) {
  # create each of the 3 panels separately and store in gplot
  for(i in 1:3) {
    print(paste("CYCLE:",i,j))
    params <- gdata[[i]] %>%
      filter(casemix_spec == j) %>%
      group_by(specialty_spec, rvu.target_spec, casemix_spec) %>%
      summarise(sd = sd(n.blocks),
                mean = mean(n.blocks),
                ymin = mean - 1.96*sd,
                ymax = mean + 1.96*sd) %>%
      ungroup() %>%
      mutate(specialty_spec = as.character(specialty_spec)) %>%
      # rename "Otoloaryngology (ENT)" as just "ENT"
      mutate(specialty_spec =
               str_replace_all(as.character(specialty_spec),
                               "Otolaryngology \\(ENT\\)", c("ENT")))

    # store panel in gplot (list)
    gplot[[i]] <- ggplot() + xlim(3000,12000) + ylim(25,290) +
      geom_line(data = params,
                aes(x = rvu.target_spec, y = mean, color = specialty_spec),
                size = 1) +
      geom_ribbon(data = params,
                  aes(x = rvu.target_spec, ymin = ymin,
                      ymax = ymax, fill = specialty_spec),
                  alpha = 0.2) +
      geom_vline(xintercept = 12000) +
      geom_vline(xintercept = 3000) +
      labs(x = "RVU Target",
           y = "Blocks Required",
           color = "Specialty",
           fill = "Specialty") +
      scale_fill() +
      scale_color() +
      annotate("text",x=3500,y=290,label=c("A","B","C")[[i]]) +
      theme_bw() +
      theme(legend.position = c(0.42, c(0.82,0.845,0.845)[[i]]),
            legend.text=element_text(size=rel(0.7)),
            legend.title = element_text(size = rel(0.7)))
  }

  # save PNG file of all three panels
  png(paste0(root_figures,'paper1/p1figure2_blocks_',
             c("All","Low","Medium","High")[[j+1]],'.png'),
      width = 6.8, height=5.8, units = "in", res = 300)
  # arrange all three panels in a 1x3 grid
  grid.arrange(grobs = gplot, ncol = 3)
  graphics.off()
}


####  PAPER 1 - FIGURE 3 - Block Reqs to reach RVU Target by specialty, ####
####  Variation with Turnover Time                                      ####

# import relevant data for 3 panel figure using parallel processing.
# gNdata nomenclature signifies grob1, grob2, grob3 corresponding to each of the panels
plan(multiprocess, workers = 3)

f13data <- future_pmap(expand.grid(.specialty = speclist,
                                   .cluster_size = 10,
                                   .casemix = casemixlist,
                                   .block_hours = 8,
                                   .turnover_time = c(50,60),
                                   .rvu_goal = seq(3000,12000,1000)),
                       import.summary, .progress = TRUE) %>% bind_rows()

for(j in casemixlist) {
  f13data %>%
    filter(casemix_spec == j) %>%
    filter(rvu.target_spec == 6000) %>%
    mutate(turnover_time = as.factor(turnover.time_spec)) %>%
    ggplot() + geom_boxplot(aes(group=interaction(specialty_spec,
                                                  as.factor(turnover.time_spec)),
                                x = specialty_spec,
                                y = n.blocks,
                                color = as.factor(turnover.time_spec))) +
    labs(y = "Blocks", x = "",
         colour = "Turnover Time")+
    geom_hline(yintercept = 48, linetype="dashed")+
    geom_hline(yintercept = 96, linetype="dashed")+
    geom_hline(yintercept = 144, linetype="dashed")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(plot.margin = unit(c(20,20,10,10),"pt")) +
    scale_x_discrete(breaks=c(cardiac, general, gyne, neuro, ortho,
                              ent, plastics, thoracic, urology, vascular),
                     labels=c("CARD", "GEN", "GYNE","NEUR","ORTH",
                              "ENT","PRS","THOR","URO","VASC")) +
    scale_color()+
    theme_bw() +
    theme(legend.position = c(0.1,0.85))
  ggsave(paste0(root_figures,'paper1/p1figure3_turnover_',
                c("All","Low","Medium","High")[[j+1]],'.png'),
         width = 6.8, height=5.8, units = "in")
}


####  PAPER 1 - FIGURE 4 - Block Reqs to Reach RVU target by specialty, casemix  ####

png(paste0(root_figures,'paper1/p1figure4_casemix.png'),
    width = 6.8, height=5.8, units = "in", res = 300)

f13data %>%
  filter(turnover.time_spec == 60) %>%
  filter(rvu.target_spec == 6000) %>%
  filter(casemix_spec != Casemix.ALL) %>%
  mutate(turnover_time = as.factor(turnover.time_spec)) %>%
  ggplot() + geom_boxplot(aes(group=interaction(specialty_spec,
                                                as.factor(casemix_spec)),
                              x = specialty_spec,
                              y = n.blocks,
                              color = as.factor(casemix_spec))) +
  labs(y = "Blocks", x = "",
       colour = "Complexity")+ 
  geom_hline(yintercept = 48, linetype="dashed")+
  geom_hline(yintercept = 96, linetype="dashed")+
  geom_hline(yintercept = 144, linetype="dashed")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(plot.margin = unit(c(20,20,10,10),"pt")) +
  scale_x_discrete(breaks=c(cardiac, general, gyne, neuro, ortho,
                            ent, plastics, thoracic, urology, vascular),
                   labels=c("CARD", "GEN", "GYNE","NEUR","ORTH",
                            "ENT","PRS","THOR","URO","VASC")) +
  scale_color_jama(labels = c("Low", "Medium","High")) +
  theme_bw() +
  theme(legend.position = c(0.1,0.86))

graphics.off()

##############################################
####  SUPPLEMENTAL FIGURES                ####
##############################################

####  SUPPLEMENTAL FIGURE SET 1A - HEATMAP SURGEON TIME HRVUS  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    alldata %>%
      filter(is.na(rvu_goal)) %>%
      filter(casemix == j) %>%
      filter(cluster_size == i) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
      geom_raster(aes(fill = rvus.hourly.util.mean))+
      geom_text(aes(label=round(rvus.hourly.util.mean,1)), size = 2) +
      scale_fill_gradient(low = "#0091ff", high = "#f0650e",
                          limits=c(0, 25), name = 'Hourly RVUs') +
      labs(x = "Block Length (Hours)",
           y = "Turnover Time (Minutes)")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5)
    ggsave(paste0(root_figures,'paper2/p2figure1A_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 1B - HEATMAP ADMIN TIME HRVUS  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    alldata %>%
      filter(is.na(rvu_goal)) %>%
      filter(casemix == j) %>%
      filter(cluster_size == i) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
      geom_raster(aes(fill = rvus.hourly.admin.mean))+
      geom_text(aes(label=round(rvus.hourly.admin.mean,1)), size = 2) +
      scale_fill_gradient(low = "#0091ff", high = "#f0650e",
                          limits=c(0, 25), name = 'Hourly RVUs') +
      labs(x = "Block Length (Hours)",
           y = "Turnover Time (Minutes)")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5)
    ggsave(paste0(root_figures,'paper2/p2figure1B_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 1C - HEATMAP RESERVED TIME HRVUS  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    alldata %>%
      filter(is.na(rvu_goal)) %>%
      filter(casemix == j) %>%
      filter(cluster_size == i) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
      geom_raster(aes(fill = rvus.hourly.resv.mean))+
      geom_text(aes(label=round(rvus.hourly.resv.mean,1)), size = 2) +
      scale_fill_gradient(low = "#0091ff", high = "#f0650e",
                          limits=c(0, 25), name = 'Hourly RVUs') +
      labs(x = "Block Length (Hours)",
           y = "Turnover Time (Minutes)")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5)
    ggsave(paste0(root_figures,'paper2/p2figure1C_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 2 - ADMIN AND SURGEON HRVUS VS BLOCK LENGTH  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    
    #get appropriate sample size for error bar calculations
    n_obs <- ifelse(i == Sched.ANNUAL, 100,1000)
    
    hrvus.admin2 <- hrvus.admin %>%
      filter(cluster_size == i) %>%
      filter(casemix == j) %>%
      filter(is.na(turnover_time)) %>%
      filter(is.na(block_hours))
    
    alldata %>%
      filter(cluster_size == i) %>%
      filter(casemix == j) %>%
      filter(rvu_goal == 12000) %>%
      filter(turnover_time == 60) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot() +
      scale_fill() +
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = rvus.hourly.util.mean - 1.96*rvus.hourly.util.sd/sqrt(n_obs),
                        ymax = rvus.hourly.util.mean + 1.96*rvus.hourly.util.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = rvus.hourly.util.mean),
                 size = 5, color = "#df8f44ff") + #orange
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = rvus.hourly.admin.mean - 1.96*rvus.hourly.admin.sd/sqrt(n_obs),
                        ymax = rvus.hourly.admin.mean + 1.96*rvus.hourly.admin.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = rvus.hourly.admin.mean),
                 size = 5, color = "#00a1d5ff") + #blue
      # geom_errorbar(aes(x = as.factor(block_size/60),
      #                   ymin = rvus.hourly.resv.mean - 1.96*rvus.hourly.resv.sd/sqrt(n_obs),
      #                   ymax = rvus.hourly.resv.mean + 1.96*rvus.hourly.resv.sd/sqrt(n_obs)),
      #               width = 0.5) +
      # geom_point(aes(x = as.factor(block_size/60), y = rvus.hourly.resv.mean),
      #            size = 5, color = "#6a6599ff") + #purple
      
      labs(x = "Block Length (Hours)",
           y = "Hourly RVUs",
           color = "Turnover Time")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5) +
      theme(legend.position = c(0.15, 0.2),
            legend.text=element_text(size=rel(1.0)),
            legend.title = element_text(size = rel(1.0)))
    ggsave(paste0(root_figures,'paper2/p2figure2_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 3A - HEATMAP PERCENT EFFICIENCY OF BLOCK UTILIZATION  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    alldata %>%
      filter(is.na(rvu_goal)) %>%
      filter(casemix == j) %>%
      filter(cluster_size == i) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
      geom_raster(aes(fill = time.or.percent_admin.mean))+
      geom_text(aes(label=round(time.or.percent_admin.mean,2)), size = 2) +
      scale_fill_gradient(low = "#0091ff", high = "#f0650e",
                          limits=c(0.5, 1.5), name = '% Utilization') +
      labs(x = "Block Length (Hours)",
           y = "Turnover Time (Minutes)")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5)
    ggsave(paste0(root_figures,'paper2/p2figure3A_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 3B - HEATMAP PERCENT EFFICIENCY OF BLOCK UTILIZATION  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    alldata %>%
      filter(is.na(rvu_goal)) %>%
      filter(casemix == j) %>%
      filter(cluster_size == i) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
      geom_raster(aes(fill = time.or.percent_resv.mean))+
      geom_text(aes(label=round(time.or.percent_resv.mean,2)), size = 2) +
      scale_fill_gradient(low = "#0091ff", high = "#f0650e",
                          limits=c(0.5, 1.5), name = '% Utilization') +
      labs(x = "Block Length (Hours)",
           y = "Turnover Time (Minutes)")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5)
    ggsave(paste0(root_figures,'paper2/p2figure3B_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 3C - HEATMAP PERCENT EFFICIENCY OF BLOCK UTILIZATION  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    alldata %>%
      filter(is.na(rvu_goal)) %>%
      filter(casemix == j) %>%
      filter(cluster_size == i) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
      geom_raster(aes(fill = time.resv.percent_util.mean))+
      geom_text(aes(label=round(time.resv.percent_util.mean,2)), size = 2) +
      scale_fill_gradient(low = "#0091ff", high = "#f0650e",
                          limits=c(0.5, 1.5), name = '% Utilization') +
      labs(x = "Block Length (Hours)",
           y = "Turnover Time (Minutes)")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5)
    ggsave(paste0(root_figures,'paper2/p2figure3C_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE SET 4 - PERCENT EFFICIENCY OF BLOCK UTILIZATION  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    
    #get appropriate sample size for error bar calculations
    n_obs <- ifelse(i == Sched.ANNUAL, 100,1000)
    
    alldata %>%
      filter(cluster_size == i) %>%
      filter(casemix == j) %>%
      filter(is.na(rvu_goal)) %>%
      filter(turnover_time == 60) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot() +
      scale_fill() +
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = time.or.percent_admin.mean - 1.96*time.or.percent_admin.sd/sqrt(n_obs),
                        ymax = time.or.percent_admin.mean + 1.96*time.or.percent_admin.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = time.or.percent_admin.mean),
                 size = 5, color = "#00a1d5ff") + #blue
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = time.resv.percent_util.mean - 1.96*time.resv.percent_util.sd/sqrt(n_obs),
                        ymax = time.resv.percent_util.mean + 1.96*time.resv.percent_util.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = time.resv.percent_util.mean),
                 size = 5, color = "#df8f44ff") + #orange
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = time.or.percent_resv.mean - 1.96*time.or.percent_resv.sd/sqrt(n_obs),
                        ymax = time.or.percent_resv.mean + 1.96*time.or.percent_resv.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = time.or.percent_resv.mean),
                 size = 5, color = "#b24745ff") + #red
      labs(x = "Block Length (Hours)",
           y = "Rate",
           color = "Turnover Time")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5) +
      theme(legend.position = c(0.15, 0.2),
            legend.text=element_text(size=rel(1.0)),
            legend.title = element_text(size = rel(1.0)))
    ggsave(paste0(root_figures,'paper2/p2figure4_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}


####  SUPPLEMENTAL FIGURE SET 5 - PERCENT OVERTIME  ####

for(j in casemixlist) {
  for(i in c(5,10,20,Sched.ANNUAL)) {
    
    #get appropriate sample size for error bar calculations
    n_obs <- ifelse(i == Sched.ANNUAL, 100,1000)
    
    alldata %>%
      filter(cluster_size == i) %>%
      filter(casemix == j) %>%
      filter(is.na(rvu_goal)) %>%
      filter(turnover_time == 60) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot() +
      scale_fill() +
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = p.blocks.ot.actual.mean - 1.96*p.blocks.ot.actual.sd/sqrt(n_obs),
                        ymax = p.blocks.ot.actual.mean + 1.96*p.blocks.ot.actual.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = p.blocks.ot.actual.mean),
                 size = 5, color = "#b24745ff") + #red
      geom_errorbar(aes(x = as.factor(block_size/60),
                        ymin = ot.actual.percent_util.mean - 1.96*ot.actual.percent_util.sd/sqrt(n_obs),
                        ymax = ot.actual.percent_util.mean + 1.96*ot.actual.percent_util.sd/sqrt(n_obs)),
                    width = 0.5) +
      geom_point(aes(x = as.factor(block_size/60), y = ot.actual.percent_util.mean),
                 size = 5, color = "#79af97ff") + #green
      labs(x = "Block Length (Hours)",
           y = "Rate",
           color = "Turnover Time")+
      theme_bw() +
      facet_wrap(.~specialty, nrow = 2, ncol = 5) +
      theme(legend.position = c(0.15, 0.2),
            legend.text=element_text(size=rel(1.0)),
            legend.title = element_text(size = rel(1.0)))
    ggsave(paste0(root_figures,'paper2/p2figure5_cluster-',i,'_',
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.8, height=5.8, units = "in")
  }
}

# ####  PAPER 2 - FIGURE 4ALT - HEATMAP PERCENT OVERTIME  ####
## percent overtime is invariant with turnover time in the binpack sim,
## so no point using heatmap when dot plot works (x=block size)
# 
# for(j in casemixlist) {
#   print(alldata %>%
#           # filter(specialty == general) %>%
#           # filter(schedule == Sched.SOFTCAP) %>%
#           # filter(speed == Speed.MEDIUM) %>%
#           filter(rvu_goal == 12000) %>%
#           filter(casemix == j) %>%
#           filter(cluster_size == 10) %>%
#           # filter(turnover_time < 30 | mod(turnover_time,30)==0) %>%
#           ggplot(aes(as.factor(block_size/60), as.factor(turnover_time))) +
#           ### geom_raster(aes(fill = rvus.hourly.util.mean - rvus.hourly.admin.mean))+
#           geom_raster(aes(fill = p.blocks.ot.actual.mean))+
#           geom_text(aes(label=round(p.blocks.ot.actual.mean,1),  size = 0.8)) +
#           # scale_fill_jama() +
#           scale_fill_gradient(low = "#0091ff", high = "#f0650e", limits=c(0, 1), name = '% Blocks Overtime') +
#           # scale_y_reverse(limits=c(90,0), expand=c(0,0)) +
#           labs(#title = 'RVUs per Hour of Block Utilization',
#             #subtitle = paste(c("All","Low","Medium","High")[[j+1]],'Casemix Complexity'),
#             x = "Block Length (Hours)",
#             y = "Turnover Time (Minutes)",
#             color = "Turnover Time")+
#           theme_bw() +
#           facet_wrap(.~specialty, nrow = 2, ncol = 5))
#   # theme(legend.position = c(0.15, 0.2),
#   #       legend.text=element_text(size=rel(1.0)),
#   #       legend.title = element_text(size = rel(1.0)),))
#   # ggsave(paste0(root_figures,'paper2/p2figure4-1_',c("All","Low","Medium","High")[[j+1]],'.png'), width = 6.8, height=5.8, units = "in")
# }

####  SUPPLEMENTAL FIGURE SET 6 - EFFECT OF CLUSTER SIZE ON BLOCK REQS  ####

for(j in casemixlist) {
  for(k in speclist) {
    alldata %>%
      filter(casemix == j) %>%
      filter(specialty == k) %>%
      filter(!is.na(rvu_goal)) %>%
      mutate(cluster_size = as.factor(cluster_size)) %>%
      filter(turnover_time %% 30 == 0) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot() + xlim(3000,12000) + ylim(0,600) +
      geom_line(aes(x = rvu_goal,
                    y = n.blocks.mean,
                    color = cluster_size,
                    group = cluster_size),
                size = 1) +
      geom_ribbon(aes(x = rvu_goal,
                      ymin = n.blocks.mean - 1.96*n.blocks.sd,
                      ymax = n.blocks.mean + 1.96*n.blocks.sd,
                      group = cluster_size,
                      fill = cluster_size),
                  alpha = 0.2) +
      geom_vline(xintercept = 12000) +
      geom_vline(xintercept = 3000) +
      
      labs(title = 'Block Requirement Vs. Cluster Size',
           subtitle = paste0(k,', ',c("All","Low","Medium","High")[[j+1]],
                             ' Casemix Complexity'),
           x = "RVU Target",
           y = "Blocks",
           color = "Cluster Size",
           fill = "Cluster Size") +
      scale_fill() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      facet_grid(turnover_time ~ block_size) +
      theme(legend.text=element_text(size=rel(1.0)),
            legend.title = element_text(size = rel(1.0)))
    ggsave(paste0(root_figures,'supplement/sfigure1_',k,"_",
                  c("All","Low","Medium","High")[[j+1]],'.png'),
           width = 6.5, height=8, units = "in")
  }
}

####  SUPPLEMENTAL FIGURE 7 - EFFECT OF CASEMIX ON BLOCK REQS  ####

for(j in c(5,10,20,Sched.ANNUAL)) {
  # for(i in casemixlist) {
  for(k in speclist) {
    alldata %>%
      filter(cluster_size == j) %>%
      filter(specialty == k) %>%
      filter(turnover_time %% 30 == 0) %>%
      filter(!is.na(rvu_goal)) %>%
      mutate(cluster_size = as.factor(cluster_size)) %>%
      mutate(casemix = as.factor(recode(casemix, '0'='All', '1'='Low',
                                        '2'='Medium', '3'='High'))) %>%
      
      # filter(turnover_time == 60) %>%
      # filter(block_size == i*60) %>%
      mutate(specialty = as.character(specialty)) %>%
      mutate(specialty =
               str_replace_all(as.character(specialty),
                               "Otolaryngology \\(ENT\\)", c("Otolaryngology"))) %>%
      ggplot() + xlim(3000,12000) + ylim(0,600) +
      geom_line(aes(x = rvu_goal,
                    y = n.blocks.mean,
                    color = casemix,
                    group = casemix),
                size = 1) +
      geom_ribbon(aes(x = rvu_goal,
                      ymin = n.blocks.mean - 1.96*n.blocks.sd,
                      ymax = n.blocks.mean + 1.96*n.blocks.sd,
                      group = casemix,
                      fill = casemix),
                  alpha = 0.2) +
      geom_vline(xintercept = 12000) +
      geom_vline(xintercept = 3000) +
      
      labs(title = 'Block Requirement Vs. Case-Mix Complexity',
           subtitle = paste0(k,', Cluster Size ',j),
           x = "RVU Target",
           y = "Blocks",
           color = "Casemix Complexity",
           fill = "Casemix Complexity") +
      scale_fill() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      facet_grid(turnover_time ~ block_size) +
      theme(legend.text=element_text(size=rel(1.0)),
            legend.title = element_text(size = rel(1.0)))
    ggsave(paste0(root_figures,'supplement/sfigure2_',k,'_cluster-',j,'.png'),
           width = 6.5, height=8, units = "in")
  }
}

# ####  SUPPLEMENTAL FIGURES 3 - DOT PLOT HRVUS ####
# 
# for(j in casemixlist) {
#   for(i in speclist){
#     print(alldata %>%
#             filter(specialty == i) %>%
#             # filter(schedule == "No new cases after soft cap") %>% #Sched.SOFTCAP
#             # filter(speed == 2) %>% #Speed.MEDIUM
#             filter(cluster_size == 10) %>%
#             filter(casemix == j) %>%
#             filter(rvu_goal == 12000) %>%
#             filter(turnover_time < 30 | mod(turnover_time,30)==0) %>%
#             ggplot() + 
#             geom_errorbar(aes(x = as.factor(block_size/60), ymin = rvus.hourly.util.mean - 1.96*rvus.hourly.util.sd, ymax = rvus.hourly.util.mean + 1.96*rvus.hourly.util.sd), width = 0.2) +
#             geom_point(aes(x = as.factor(block_size/60), y = rvus.hourly.util.mean,color = as.factor(turnover_time)), size = 10) +
#             scale_color() +
#             labs(x = "Block Length (Hours)",
#                  y = "Hourly RVUs (Surgeon Time)",
#                  color = "Turnover Time")+
#             theme_bw() +
#             theme(legend.position = c(0.1, 0.1), legend.text=element_text(size=rel(0.6)), legend.title = element_text(size = rel(0.6))))
#     # ggsave(paste0(root_figures,'paper2/p2figure1_',i,'_',c("All","Low","Medium","High")[[j+1]],'.png'), width = 6.8, height=5.8, units = "in")
#   }
# }



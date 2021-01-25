# nsqip-analysis
Clean, filter, and simulate OR usage for surgeons in varying specialties and practice environments based on cases in the ACS NSQIP dataset. Shared in advance of expected publication of Rao et al. (2021).

Project workflow:
1. sim_binpack.R
  No project dependencies.
  Run the simulation on an HPC cluster to generate operating schedules for surgeons in various practice environments.
2. summary_table_binpack.R
  Dependent on output from sim_binpack.R
  Calculate summary statistics for data from each simulation.
  Calculate summary statistics on subgroups sharing one or more simulation parameters in common.
3. linear_models_binpack.R
  Dependent on output from sim_binpack.R
  Calculate multivariate linear models for means of several metrics of interest.
4. linear_models_sd_binpack.R
  Dependent on output from summary_table_binpack.R
  Calculate multivariate linear models for standard deviations of several metrics of interest.
5. anovas_binpack.R
  Dependent on output from sim_binpack.R
  Calculate ANOVA across simulation data differeing in one or more simulation parameters.
6. anovas_summary_binpack.R
  Dependent on output from anovas_binpack.R
  Summarize frequency of significant differences across ANOVA tests on simulation data sharing one or more parameters.
7. figures_manuscript_binpack.R
  Dependent on output from sim_binpack.R AND summary_table_binpack.R
  Construct figures for manuscript and supplement.
8. nsqip_summary_table_binpack.R
  No project dependencies.
  Record numbers for STROBE diagram and reporting guidelines.

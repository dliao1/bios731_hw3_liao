The workflow of the my simulation study is as follows:

- All the main simulation code is in simulations/run_sim.R. Running that R script should generate all simulation data across all given parameter combinations and save intermediate results for the current parameter combination to a .RDA file in results/sim_[method]/scenario_[i].RDA . For example, the 3rd parameter combination calculating bootstrap percentile confidence intervals would be in results/sim_boot_percentile/scenario_3.RDA.
- When you run run_sim.R it will call 5 other helper functions to actually generate data and extract estimates.
- Helper functions are stored in the source folder. There are 6 R scripts in there: gen_data.R generates data given an n, beta, and error distribution, fit_lm.R fits a linear regression model to a given dataset, run_bootstrap.R generates a specified number of bootstrap samples, get_estimates.R extracts Wald/Bootstrap percentile estimates given a dataset, process_results.R aggregates results from each parameter combination into one dataframe
  - Tables and figures summarizing bias, coverage, and computation time from my simulation study can be found in analysis/HW3_simulations.pdf. I have included both the knitted pdf and the RMD file in the analysis folder.




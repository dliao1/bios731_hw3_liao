#!/bin/bash
#SBATCH --array=1-4%2
#SBATCH --job-name=run_sim_job_dliao6
#SBATCH --partition=encore
#SBATCH --output=output_run_sim
#SBATCH --error=error_run_sim

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript run_sim.R $JOBID
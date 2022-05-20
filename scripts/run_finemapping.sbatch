#!/bin/bash

#SBATCH --job-name=finemapping
#SBATCH --output=finemapping_%A_%a.out
#SBATCH --error=finemapping_%A_%a.err
#SBATCH --array=1-122
#SBATCH --partition=broadwl
#SBATCH --account=pi-xinhe
#SBATCH --mem=10G

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load R/4.0.4

cd /home/kaixuan/projects/analysis_pipelines/scripts

Rscript run_finemapping.R $SLURM_ARRAY_TASK_ID

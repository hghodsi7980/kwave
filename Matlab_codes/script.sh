#!/bin/sh
#
#SBATCH --job-name="clot_generator"
#SBATCH --account=Research-ME-BME
#SBATCH --partition=compute
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=96G
#SBATCH --array=1-210

module load matlab

# Calculate the iteration range for this job
startIter=$(( ($SLURM_ARRAY_TASK_ID - 1) * 10 + 1 ))
endIter=$(( $SLURM_ARRAY_TASK_ID * 10 ))

srun matlab -batch "PAsim(${startIter},${endIter});exit;"



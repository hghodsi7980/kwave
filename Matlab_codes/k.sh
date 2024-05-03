#!/bin/sh
#
#SBATCH --job-name="k"
#SBATCH --account=Research-ME-BME
#SBATCH --partition=compute-p2
#SBATCH --time=60:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=3G


module load matlab

srun matlab -batch "run('/home/hghodsi/Matlab_codes/kwavetest_1.m');exit;"
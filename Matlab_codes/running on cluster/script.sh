#!/bin/sh
#
#SBATCH --job-name="clot"
#SBATCH --account=Research-ME-BME
#SBATCH --partition=compute-p2
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=196G
#SBATCH --array=1-685  # Adjust this based on half the number of missing files

module load matlab

# Read missing files into an array
missing_files=($(cat missing_files.txt))

# Get the index of the current task
index=$((SLURM_ARRAY_TASK_ID - 1))

# Calculate the start and end indices for this task
start_index=$((index * 2))
end_index=$((start_index + 1))

# Get the file indices to process
if [ $start_index -lt ${#missing_files[@]} ]; then
    file_index1=${missing_files[$start_index]}
    if [ $end_index -lt ${#missing_files[@]} ]; then
        file_index2=${missing_files[$end_index]}
        srun matlab -batch "PAsim(${file_index1},${file_index1}); PAsim(${file_index2},${file_index2});exit;"
    else
        srun matlab -batch "PAsim(${file_index1},${file_index1});exit;"
    fi
fi

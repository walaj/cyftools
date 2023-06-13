#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8G                          # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

input_file=$1
output_file=$2

module load R

if ! command -v Rscript &> /dev/null
then
    echo "Rscript could not be found"
    exit
fi

## parameters
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "...running: Rscript /home/jaw34/git/cysift/R/frame_process.R ${input_file} ${output_file}" 
    Rscript /home/jaw34/git/cysift/R/frame_process.R ${input_file} ${output_file}
fi

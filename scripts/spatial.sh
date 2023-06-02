#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=32G                          # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

input_file=$1
output_file=$2

## parameters

# max number of KNN nodes
K=10000
# max distance in pixels
D=1000
# max number of cores
T=1

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "...running: cysift spatial -k $K -d $D -t $T -v ${input_file} - | gzip > ${output_file}"
    cysift spatial -k $K -d $D -t $T -v ${input_file} - | gzip > ${output_file}
fi

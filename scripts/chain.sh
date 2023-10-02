#!/bin/bash

#SBATCH -c 4                               # Request four cores
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=12G                          # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

input_file=$1
output_file=$2
pheno_file=$3

function check_file_exists {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "File not found: $file"
        exit 1
    fi
}

T=4
V="-v"

base="${input_file%%.*}"
echo "BASE $base"

# radial file
RAD=/home/jaw34/projects/orion/radial.csv

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "...running: cysift chain on ${base}"
    echo "cysift pheno ${V} $input_file -t $pheno_file - | cysift tumor ${V} - - -o 131072 -f 0.50 -k 25 -t ${T} | cysift radialdens ${V} - - -t ${T} -f ${RAD} | cysift delaunay - ${output_file} -l 20" 
    cysift pheno ${V} $input_file -t $pheno_file - |\
	cysift tumor ${V} - - -o 4 -f 0.50 -k 25 -t ${T} |\
	cysift radialdens ${V} - - -t ${T} -f ${RAD} |\
	cysift delaunay -l 20 - ${output_file}
fi

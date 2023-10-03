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

source ~/git/cysift/scripts/config.sh

T=4
V="-v"

base="${input_file%%.*}"

# radial file

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
elif contains_string "$orion41_73" "$input_file"; then #[[ $orion41_73 == *"$base"* ]]; then
    echo "chain.sh: detected Orion 41-73"
    RAD=/home/jaw34/projects/orion/radial.csv    
elif contains_string "$orion1_40" "$input_file"; then #[[ $orion1_40 == *"$base"* ]]; then
    echo "chain.sh: detected Orion 1-40"
    RAD=/home/jaw34/projects/orion/radial.csv
elif [[ "$input_file" == *"immune"* ]]; then
    echo "chain.sh: detected CyCIF Immune. Need radial file"
    exit 1
elif [[ "$input_file" == *"tumor"* ]]; then
    echo "chain.sh: detected CyCIF Tumor. Need radial file"
    exit 1
elif contains_string "$prostate" "$input_file"; then
    echo "chain.sh: detected Prostate"
    RAD=/home/jaw34/projects/prostate/radial.csv    
else
    echo "header.sh: Warning: $input_file doesn't fit into cycif, prostate, orion, etc"
    exit 1
fi


if [[ ! -f "$input_file" ]]; then
    echo "Error in chain.sh: File '$input_file' does not exist."
    exit 1
else
    echo "...running: cysift chain on ${base}"
    cmd="cysift pheno ${V} $input_file -t $pheno_file - |\
    		      cysift filter - - -a 4 ${V} |\
		      cysift tumor - - -f 0.50 -k 25 -t ${T} ${V} |\
		      cysift radialdens ${V} - - -t ${T} -f ${RAD} |\
		      cysift delaunay -l 20 - ${output_file}"
    echo "$cmd"
    eval "$cmd"
fi

#!/bin/bash

#SBATCH -c 4                               # Request four cores
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                           # Memory total in MiB (for all cores)
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
    cysift pheno ${V} $input_file -t $pheno_file - |\
	cysift tumor ${V} - - -o 131072 -f 0.50 -k 25 -t ${T} |\
	cysift radialdens ${V} - - -t ${T} -f ${RAD} |\
        cysift delaunay - ${base}.ptrd.cys -l 20

    # check_file_exists ${base}.ptrd.cys
    
    # echo "...running: cysift chain tumor on ${base}"    
    # cysift log10 ${base}.ptrd.cys - |\
    # 	cysift select ${V} - - -O 1 |\
    # 	cysift umap - ${V} ${base}.tumor.cys -t ${T}
    # echo "...running: cysift chain 2 on $base"
    # cysift select ${V} ${base}.ptrd.cys - -o 2048 -O 1 |\
    # 	cysift log10 - - |\
    # 	cysift umap - ${V} ${base}.PDL1.cys -t ${T}
    # echo "...running: cysift chain 3 on $base"    
    # cysift select ${V} ${base}.ptrd.cys - -a 133120 -O 1 |\
    # 	cysift log10 - - |\
    # 	cysift umap - ${V} ${base}.PDL1_PanCK.cys -t ${T}
    # echo "...running: cysift chain 4 on $base"    
    # cysift select ${V} ${base}.ptrd.cys - -O 1 -a 2064 |\
    # 	cysift log10 - - |\
    # 	cysift umap - ${V} ${base}.PDL1_CD68.cys -t ${T}
fi

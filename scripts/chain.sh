#!/bin/bash

#SBATCH -c 4                               # Request four cores
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=12G                          # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

input_file=$1
output_file=$2
pheno_file=$3
roi_file=$4

source ~/git/cyftools/scripts/config.sh

TLS_FLAG=32
T=6
#V="-v"

base="${input_file%%.*}"

# radial file

if [ -d /Users ]; then
    PROJ_HOME=/Users/jeremiahwala/Sorger/projects
    PROJ_DATA=/Users/jeremiahwala/Sorger/projects    
else
    PROJ_HOME=/home/jaw34/projects/
    PROJ_DATA=/n/scratch3/users/j/jaw34/
fi    

if [[ ! -f "$input_file" ]]; then
    parallel_echo "Error: File '$input_file' does not exist."
    exit 1
elif contains_string "$orion41_73" "$input_file"; then
    parallel_echo "...chain.sh: detected Orion 41-73"
    RAD=${PROJ_HOME}/orion/radial.csv
    TUMOR_MARKER=131072
    TCELL_MARKER=4096
elif contains_string "$orion1_40" "$input_file"; then
    parallel_echo "...chain.sh: detected Orion 1-40"
    RAD=${PROJ_HOME}/orion/radial.csv    
    TUMOR_MARKER=131072
    TCELL_MARKER=4096    
elif [[ "$input_file" == *"immune"* ]]; then
    parallel_echo "...chain.sh: detected CyCIF Immune. Need radial file"
    exit 1
elif [[ "$input_file" == *"tumor"* ]]; then
    parallel_echo "...chain.sh: detected CyCIF Tumor. Need radial file"
    exit 1
elif contains_string "$prostate" "$input_file"; then
    parallel_echo "...chain.sh: detected Prostate"
    RAD=${PROJ_HOME}/prostate/radial.csv
    TUMOR_MARKER=4
    TCELL_MARKER=2048
    BCELL_MARKER=64
    IMMUNE_MARKER=268435408

else
    parallel_echo "...header.sh: Warning: $input_file doesn't fit into cycif, prostate, orion, etc"
    exit 1
fi

## set the roi file
if [[ -f "$roi_file" ]]; then
    roicmd="cyftools roi - - -b -m 0.325 -r $roi_file |"
else
    roicmd=""
fi

if [[ ! -f "$input_file" ]]; then
    parallel_echo "Error in chain.sh: File '$input_file' does not exist."
    exit 1
elif [[ "$input_file" == *"LSP10388"* || "$input_file" == *"LSP10353"* || "$input_file" == *"LSP10375"* || "$input_file" == *"LSP10364"* ]]; then
    # Your slightly different command here
    parallel_echo "Running the rescale version for LSP10388, LSP10353, LSP10375, or LSP10364"
    cmd="cyftools magnify $input_file -f 1.015625 - | cyftools pheno ${V} - -t $pheno_file - |
    		      cyftools filter - - -a $TUMOR_MARKER ${V} -M |
		      cyftools annotate - - -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 | $roicmd
		      cyftools filter - - -a $TCELL_MARKER ${V} -M |
		      cyftools annotate - - -f 0.50 -k 25 -d 100000 -t ${T} ${V} -F 16 | 
		      cyftools island - - -n 5000 -T | cyftools island - - -S -n 5000 |
		      cyftools margin -d 100 - - |
		      cyftools radialdens ${V} - - -t ${T} -f ${RAD} |
		      cyftools delaunay -l 20 - ${output_file}"
    echo "$cmd" | tr '\n' ' '
    echo ""    
    eval "$cmd"
else
    parallel_echo "...running: cyftools chain on ${base}"
    cmd="cyftools pheno ${V} $input_file -t $pheno_file - |
    		      cyftools filter - - -a $TUMOR_MARKER ${V} -M | cyftools annotate - - -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 | $roicmd
		      cyftools tls - - -b $BCELL_MARKER -i $IMMUNE_MARKER -m 300 -d 35 ${V} |
 		      cyftools island - - -n 5000 -T | cyftools island - - -S -n 5000 | cyftools margin -d 100 - - |
           	      cyftools radialdens ${V} - - -t ${T} -f ${RAD} |  cyftools delaunay -l 20 - ${output_file}"
    echo "$cmd" | tr '\n' ' '
    echo ""
    eval "$cmd"
fi

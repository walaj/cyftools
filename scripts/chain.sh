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
V="-v"

base="${input_file%%.*}"

# radial file

if [ -d /Users ]; then
    PROJ_HOME=/Users/jeremiahwala/Dropbox/Sorger/projects
    PROJ_DATA=/Users/jeremiahwala/Dropbox/Sorger/projects    
else
    PROJ_HOME=/home/jaw34/projects/
    PROJ_DATA=/n/scratch3/users/j/jaw34/
fi    

xoffset=0;
yoffset=0;
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
    yoffset=0
    distann="cyftools filter - - -A 8 -M | cyftools dist -i tumor - - |
 cyftools filter - - -a 2048 -M | cyftools dist -i CD3 - - |
 cyftools filter - - -a 4096 -M | cyftools dist -i CD8 - - |
 cyftools filter - - -a 32768 -M | cyftools dist -i PD1 - - |
 cyftools filter - - -a 16384 -M | cyftools dist -i FOXP3 - - |
 cyftools filter - - -a 4 -M | cyftools dist -i AMCAR - - |
 cyftools filter - - -a 8 -M | cyftools dist -i HMWCK - - |
 cyftools filter - - -a 65536 -M | cyftools dist -i CD57 - - |
 cyftools filter - - -A 4096 -M | cyftools dist -i GG1 - - |
 cyftools filter - - -A 8192 -M | cyftools dist -i GG2 - - |
 cyftools filter - - -A 16384 -M | cyftools dist -i GG3 - - |
 cyftools filter - - -A 32768 -M | cyftools dist -i GG4 - - |
 cyftools filter - - -A 65536 -M | cyftools dist -i GG5 - - |
 cyftools filter - - -A 131072 -M | cyftools dist -i PNI - - |
 cyftools filter - - -A 262144 -M | cyftools dist -i SV - - |
 cyftools filter - - -A 32 -M | cyftools dist -i TLS - - |"
else
    parallel_echo "...header.sh: Warning: $input_file doesn't fit into cycif, prostate, orion, etc"
    exit 1
fi

# ## set the offsets
# declare -A xoffsetlut=(
#     [LSP12601]=0
#     [LSP12603]=0
#     [LSP12605]=0
#     [LSP12607]=0
#     [LSP12609]=0
#     [LSP12611]=0
#     [LSP12613]=0
#     [LSP12615]=0    
# )
# declare -A yoffsetlut=(
#     [LSP12601]=-275
#     [LSP12603]=-275
#     [LSP12605]=-275
#     [LSP12607]=-275
#     [LSP12609]=-275
#     [LSP12611]=-275
#     [LSP12613]=-275    
# )

# # Extract the offset from the string
# if [[ $input_file =~ (LSP[0-9]+) ]]; then
#     PATTERN=${BASH_REMATCH[1]}
#     echo "Extracted pattern: $PATTERN"

#     # Lookup the pattern in the LUT
#     if [[ -n "${lut[$PATTERN]}" ]]; then
#         xoffset=${xoffsetlut[$PATTERN]}
#         yoffset=${yoffsetlut[$PATTERN]}	
#         echo "Found in LUT: $X"
#     else
#         # Pattern not found in LUT, default X to zero
#         xoffset=0
# 	yoffset=0
#         echo "Pattern not found in LUT. Defaulting X to $X"
#     fi
# else
#     echo "Warning: No pattern found in the string."
#     # Handle the case where no pattern is found as needed
# fi

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
 cyftools annotate - - -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 |
 $roicmd
 cyftools filter - - -a $TCELL_MARKER ${V} -M |
 cyftools annotate - - -f 0.50 -k 25 -d 100000 -t ${T} ${V} -F 16 | 
 cyftools island - - -n 5000 -T | cyftools island - - -S -n 5000 |
 cyftools margin -d 100 - - |
 $distann
 cyftools radialdens ${V} - - -t ${T} -f ${RAD} |
 cyftools delaunay -l 20 - ${output_file}"
    echo "$cmd" | tr '\n' ' '
    echo ""    
    eval "$cmd"
else
    parallel_echo "...running: cyftools chain on ${base}"
    cmd="cyftools pheno ${V} $input_file -t $pheno_file - |
 cyftools filter - - -a $TUMOR_MARKER ${V} -M |
 cyftools annotate - - -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 |
 $roicmd
 cyftools tls - - -b $BCELL_MARKER -i $IMMUNE_MARKER -m 300 -d 35 ${V} |
 cyftools island - - -n 5000 -T | cyftools island - - -S -n 5000 | cyftools margin -d 100 - - |
$distann
cyftools radialdens ${V} - - -t ${T} -f ${RAD} | cyftools delaunay -l 20 - ${output_file}"
    echo "$cmd" | tr '\n' ' '
    echo ""
    eval "$cmd"
fi

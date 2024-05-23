#!/usr/bin/env bash

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

declare -A maglut=(
    [LSP10388]=1.015625
    [LSP10353]=1.015625
    [LSP10375]=1.015625
    [LSP10364]=1.015625
    )

declare -A ymaxlut=(
    [LSP12601]=8818
    [LSP12603]=9482
    [LSP12605]=8129
    [LSP12607]=8142
    [LSP12609]=8146
    [LSP12611]=8140
    [LSP12613]=8135
    [LSP12615]=4750
    [LSP12617]=7459
    [LSP12619]=7457
    [LSP12621]=5832
    [LSP12623]=5427
    [LSP12625]=9216
    [LSP12627]=6784
    [LSP12629]=6789
    [LSP12631]=5434
    [LSP12633]=8816
    [LSP12635]=6106
    [LSP12637]=7457
    [LSP12639]=9485
    [LSP12641]=6781
    [LSP12643]=8811
    [LSP12645]=7465
    [LSP12647]=8132
    [LSP12649]=7459
    [LSP12651]=9485
    [LSP12653]=5429
    [LSP12655]=8817
    [LSP12657]=6784
)    

# Make the flip command 
if [[ $input_file =~ (LSP[0-9]+) ]]; then
    PATTERN=${BASH_REMATCH[1]}

    # Lookup the pattern in the LUT
    if [[ -n "${ymaxlut[$PATTERN]}" ]]; then
        ymax=${ymaxlut[$PATTERN]}
	flipcmd="cyftools flip - - -y 0 -Y $ymax |"
        echo "Found in LUT: $ymax"
    fi
else
    echo "Warning: No LSP pattern found in the string."
fi

# Make the magnify command
if [[ $input_file =~ (LSP[0-9]+) ]]; then
    PATTERN=${BASH_REMATCH[1]}
    
    # Lookup the pattern in the LUT
    if [[ -n "${maglut[$PATTERN]}" ]]; then
        mag=${maglut[$PATTERN]}
	magcmd="cyftools magnify - - -f $mag |"
        echo "Found in magnify LUT: $mag"
    fi
else
    echo "Warning: No LSP pattern found in the string."
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
else
    parallel_echo "...running: cyftools chain on ${base}"
    cmd="cyftools check $input_file - |
$magcmd
cyftools pheno ${V} - -t $pheno_file - |
cyftools filter - - -a $TUMOR_MARKER ${V} -M |
cyftools annotate - - -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 |
$flipcmd
$roicmd
cyftools tls - - -b $BCELL_MARKER -i $IMMUNE_MARKER -m 300 -d 35 ${V} |
cyftools island - - -n 5000 -T | cyftools island - - -S -n 5000 | cyftools margin -d 100 - - |
$distann
cyftools radialdens ${V} - - -t ${T} -f ${RAD} | cyftools delaunay -l 20 - ${output_file}"
    
    echo "$cmd" | tr '\n' ' '
    echo ""
    eval "$cmd"
fi

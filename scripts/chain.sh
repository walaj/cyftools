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

### PARAMS
# num threads
T=6
# verbose?
V="-v"

## extract the basename
base="${input_file%%.*}"

##########
# HOLD SPECIFIC FILE TRANSFORMATION PARAMS
##########
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

#######
## PROJECT SPECIFIC PARAMETERS
#######
# just some default "offsets"
xoffset=0;
yoffset=0;

## error if input file does not exist
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
## CRC Orion - samples 41-73    
elif contains_string "$orion41_73" "$input_file"; then
    echo "...chain.sh: detected Orion 41-73"
    RAD=${PROJ_HOME}/orion/radial.csv
    TUMOR_MARKER=131072
    TCELL_MARKER=4096
    roimag=0.325
    distann="cyftools filter - - -A 8 -M | cyftools dist -i tumor - - |
 cyftools filter - - -a 4096 -M | cyftools dist -i CD3 - - |
 cyftools filter - - -a 256 -M | cyftools dist -i CD8 - - |
 cyftools filter - - -a 4352 -M | cyftools dist -i CD3CD8 - - |
 cyftools filter - - -a 32768 -M | cyftools dist -i PD1 - - |
 cyftools filter - - -a 128 -M | cyftools dist -i FOXP3 - - |
 cyftools filter - - -a 131072 -M | cyftools dist -i PDL1 - - |
 cyftools filter - - -a 133120 -M | cyftools dist -i PanCKPDL1 - - |
 cyftools filter - - -a 139264 -M | cyftools dist -i CD163PDL1 - - |
 cyftools filter - - -A 32 -M | cyftools dist -i TLS - - |
 cyftools filter - - -a 1024 -M | cyftools dist -i CD20 - - |"
## CRC Orion - samples 1-40
elif contains_string "$orion1_40" "$input_file"; then
    echo "...chain.sh: detected Orion 1-40"
    RAD=${PROJ_HOME}/orion/radial.csv    
    TUMOR_MARKER=131072
    TCELL_MARKER=4096
    roimag=0.325
    distann="cyftools filter - - -A 8 -M | cyftools dist -i tumor - - |
 cyftools filter - - -a 4096 -M | cyftools dist -i CD3 - - |
 cyftools filter - - -a 256 -M | cyftools dist -i CD8 - - |
 cyftools filter - - -a 4352 -M | cyftools dist -i CD3CD8 - - |
 cyftools filter - - -a 32768 -M | cyftools dist -i PD1 - - |
 cyftools filter - - -a 128 -M | cyftools dist -i FOXP3 - - |
 cyftools filter - - -a 131072 -M | cyftools dist -i PDL1 - - |
 cyftools filter - - -a 133120 -M | cyftools dist -i PanCKPDL1 - - |
 cyftools filter - - -a 139264 -M | cyftools dist -i CD163PDL1 - - |
 cyftools filter - - -A 32 -M | cyftools dist -i TLS - - |
 cyftools filter - - -a 1024 -M | cyftools dist -i CD20 - - |"
## CRC CyCIF Immune panel
elif [[ "$input_file" == *"immune"* ]]; then
    echo "...chain.sh: detected CyCIF Immune. Need radial file"
    exit 1
## CRC CyCIF Tumor panel    
elif [[ "$input_file" == *"tumor"* ]]; then
    echo "...chain.sh: detected CyCIF Tumor. Need radial file"
    exit 1
## Prostate CYCIF
elif contains_string "$prostate" "$input_file"; then
    echo "...chain.sh: detected Prostate"
    RAD=${PROJ_HOME}/prostate/radial.csv
    TUMOR_MARKER=4
    TCELL_MARKER=2048
    BCELL_MARKER=64
    IMMUNE_MARKER=268435408
    roimag=0.325    
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
## JHU ORION    
elif contains_string "$jhu" "$input_file"; then
    echo "...chain.sh: detected JHU project -- ORION"
    RAD=${PROJ_HOME}/jhu/radial.csv
    TUMOR_MARKER=65536
    TCELL_MARKER=2048
    BCELL_MARKER=64
    IMMUNE_MARKER=35824
    yoffset=0
    roimag=1    
    distann="cyftools filter - - -A 8 -M | cyftools dist -i tumor - - |
 cyftools filter - - -a 2048 -M | cyftools dist -i CD3 - - |
 cyftools filter - - -a 256 -M | cyftools dist -i CD8 - - |
 cyftools filter - - -a 8192 -M | cyftools dist -i PD1 - - |
 cyftools filter - - -a 18432 -M | cyftools dist -i FOXP3 - - |
 cyftools filter - - -a 32 -M | cyftools dist -i CD163 - - |
 cyftools filter - - -a 16 -M | cyftools dist -i CD68 - - |
 cyftools filter - - -a 1024 -M | cyftools dist -i PDL1 - - |"
## JHU CYCIF    
elif contains_string "$jhu_cycif" "$input_file"; then
    echo "...chain.sh: detected JHU project --CYCIF"
    RAD=${PROJ_HOME}/jhu/cycif/radial.csv
    TUMOR_MARKER=1024
    roimag=1    
else
    echo "...header.sh: Warning: $input_file doesn't fit into cycif, prostate, orion, etc"
    exit 1
fi

##########
## FORM INDIVIDUAL COMMANDS
##########

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
    roicmd="cyftools roi - - -m $roimag -r $roi_file |"
else
    roicmd=""
fi

## make the tls commands
if [[ -n "$BCELL_MARKER" && -n "$IMMUNE_MARKER" ]]; then
    tlscmd="cyftools tls - - -b $BCELL_MARKER -i $IMMUNE_MARKER -m 300 -d 35 ${V} |"
fi

## make the radial density command
if [ -f "$RAD" ]; then
    radcmd="cyftools radialdens ${V} - - -t ${T} -f ${RAD} |"
fi

##########
# RUN THE CHAIN OF CYFTOOLS COMMANDS
##########

if [[ ! -f "$input_file" ]]; then
    echo "Error in chain.sh: File '$input_file' does not exist."
    exit 1
else
    echo "...running: cyftools chain on ${base}"
    cmd="cyftools check $input_file - |
$magcmd
cyftools pheno ${V} - -t $pheno_file - |
cyftools filter - - -a $TUMOR_MARKER ${V} -M |
cyftools annotate - - -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 |
$flipcmd
$roicmd
$tlscmd
cyftools island - - -n 5000 -T | cyftools island - - -S -n 5000 | cyftools margin -d 100 - - |
$distann
$radcmd
cyftools delaunay -l 20 - ${output_file}"
    
    echo "$cmd" | tr '\n' ' '
    echo ""
    eval "$cmd"
fi

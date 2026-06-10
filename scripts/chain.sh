#!/usr/bin/env bash

#SBATCH -c 4                               # Request four cores
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=12G                          # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)

if [[ "${USE_SLURM:-0}" -eq 1 ]]; then
    module load gcc/9.2.0
    module load cairo
    module load boost
fi

input_file=$1
output_file=$2
pheno_file=$3
roi_file=$4

source ~/git/cyftools/scripts/config.sh

### for chaining or debugging
IN="-"
OUT="-"
DELIM="|"

#IN="tmpin.cyf"
#OUT="tmpout.cyf"
#DELIM="&& mv tmpout.cyf tmpin.cyf"

IO="${IN} ${OUT}"
### PARAMS
# num threads
T=4
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
    [LSP20122]=0.66
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
    BCELL_MARKER=1024
    IMMUNE_MARKER=14296
    roimag=0.325
    distann="cyftools filter ${IO} -A 8 -M ${DELIM} cyftools dist ${IO} -i tumor ${DELIM}
 cyftools filter ${IO} -a 4096 -M ${DELIM} cyftools dist ${IO} -i CD3 ${DELIM}
 cyftools filter ${IO} -a 256 -M ${DELIM} cyftools dist ${IO} -i CD8 ${DELIM}
 cyftools filter ${IO} -a 4352 -M ${DELIM} cyftools dist ${IO} -i CD3CD8 ${DELIM}
 cyftools filter ${IO} -a 32768 -M ${DELIM} cyftools dist ${IO} -i PD1 ${DELIM}
 cyftools filter ${IO} -a 128 -M ${DELIM} cyftools dist ${IO} -i FOXP3 ${DELIM}
 cyftools filter ${IO} -a 131072 -M ${DELIM} cyftools dist ${IO} -i PDL1 ${DELIM}
 cyftools filter ${IO} -a 133120 -M ${DELIM} cyftools dist ${IO} -i PanCKPDL1 ${DELIM}
 cyftools filter ${IO} -a 139264 -M ${DELIM} cyftools dist ${IO} -i CD163PDL1 ${DELIM}
 cyftools filter ${IO} -A 32 -M ${DELIM} cyftools dist ${IO} -i TLS ${DELIM}
 cyftools filter ${IO} -a 1024 -M ${DELIM} cyftools dist ${IO} -i CD20 ${DELIM}"
## CRC Orion - samples 1-40
elif contains_string "$orion1_40" "$input_file"; then
    echo "...chain.sh: detected Orion 1-40"
    RAD=${PROJ_HOME}/orion/radial.csv    
    TUMOR_MARKER=131072
    TCELL_MARKER=4096
    BCELL_MARKER=1024
    IMMUNE_MARKER=14296
    roimag=0.325
    distann="cyftools filter ${IO} -A 8 -M ${DELIM} cyftools dist ${IO} -i tumor ${DELIM}
 cyftools filter ${IO} -a 4096 -M ${DELIM} cyftools dist ${IO} -i CD3 ${DELIM}
 cyftools filter ${IO} -a 256 -M ${DELIM} cyftools dist ${IO} -i CD8 ${DELIM}
 cyftools filter ${IO} -a 4352 -M ${DELIM} cyftools dist ${IO} -i CD3CD8 ${DELIM}
 cyftools filter ${IO} -a 32768 -M ${DELIM} cyftools dist ${IO} -i PD1 ${DELIM}
 cyftools filter ${IO} -a 128 -M ${DELIM} cyftools dist ${IO} -i FOXP3 ${DELIM}
 cyftools filter ${IO} -a 131072 -M ${DELIM} cyftools dist ${IO} -i PDL1 ${DELIM}
 cyftools filter ${IO} -a 133120 -M ${DELIM} cyftools dist ${IO} -i PanCKPDL1 ${DELIM}
 cyftools filter ${IO} -a 139264 -M ${DELIM} cyftools dist ${IO} -i CD163PDL1 ${DELIM}
 cyftools filter ${IO} -A 32 -M ${DELIM} cyftools dist ${IO} -i TLS ${DELIM}
 cyftools filter ${IO} -a 1024 -M ${DELIM} cyftools dist ${IO}  -i CD20 ${DELIM}"
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
distann="cyftools filter ${IO} -A 8 -M ${DELIM} cyftools dist ${IO} -i tumor ${DELIM}
cyftools filter ${IO} -a 2048 -M ${DELIM} cyftools dist ${IO} -i CD3 ${DELIM}
cyftools filter ${IO} -a 4096 -M ${DELIM} cyftools dist ${IO} -i CD8 ${DELIM}
cyftools filter ${IO} -a 32768 -M ${DELIM} cyftools dist ${IO} -i PD1 ${DELIM}
cyftools filter ${IO} -a 16384 -M ${DELIM} cyftools dist ${IO} -i FOXP3 ${DELIM}
cyftools filter ${IO} -a 4 -M ${DELIM} cyftools dist ${IO} -i AMCAR ${DELIM}
cyftools filter ${IO} -a 8 -M ${DELIM} cyftools dist ${IO} -i HMWCK ${DELIM}
cyftools filter ${IO} -a 65536 -M ${DELIM} cyftools dist ${IO} -i CD57 ${DELIM}
cyftools filter ${IO} -A 4096 -M ${DELIM} cyftools dist ${IO} -i GG1 ${DELIM}
cyftools filter ${IO} -A 8192 -M ${DELIM} cyftools dist ${IO} -i GG2 ${DELIM}
cyftools filter ${IO} -A 16384 -M ${DELIM} cyftools dist ${IO} -i GG3 ${DELIM}
cyftools filter ${IO} -A 32768 -M ${DELIM} cyftools dist ${IO} -i GG4 ${DELIM}
cyftools filter ${IO} -A 65536 -M ${DELIM} cyftools dist ${IO} -i GG5 ${DELIM}
cyftools filter ${IO} -A 131072 -M ${DELIM} cyftools dist ${IO} -i PNI ${DELIM}
cyftools filter ${IO} -A 262144 -M ${DELIM} cyftools dist ${IO} -i SV ${DELIM}
cyftools filter ${IO} -A 32 -M ${DELIM} cyftools dist ${IO} -i TLS ${DELIM}"

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
    distann="cyftools filter ${IO} -A 8 -M ${DELIM} cyftools dist ${IO} -i tumor ${DELIM}
cyftools filter ${IO} -a 2048 -M ${DELIM} cyftools dist ${IO} -i CD3 ${DELIM}
cyftools filter ${IO} -a 256 -M ${DELIM} cyftools dist ${IO} -i CD8 ${DELIM}
cyftools filter ${IO} -a 8192 -M ${DELIM} cyftools dist ${IO} -i PD1 ${DELIM}
cyftools filter ${IO} -a 18432 -M ${DELIM} cyftools dist ${IO} -i FOXP3 ${DELIM}
cyftools filter ${IO} -a 32 -M ${DELIM} cyftools dist ${IO} -i CD163 ${DELIM}
cyftools filter ${IO} -a 16 -M ${DELIM} cyftools dist ${IO} -i CD68 ${DELIM}
cyftools filter ${IO} -a 1024 -M ${DELIM} cyftools dist ${IO} -i PDL1 ${DELIM}"

## JHU REVISION
elif contains_string "$jhu_revision" "$input_file"; then
    echo "...chain.sh: detected JHU revision project"
    RAD="${PROJ_HOME}/jhu/revision/radial.csv"
    TUMOR_MARKER=4096
    TCELL_MARKER=1024
    BCELL_MARKER=0
    yoffset=0
    roimag=0.650
## JHU CYCIF    
elif contains_string "$jhu_cycif" "$input_file"; then
    echo "...chain.sh: detected JHU project --CYCIF"
    RAD=${PROJ_HOME}/jhu/cycif/radial.csv
    TUMOR_MARKER=1024
    roimag=1
    ## NEOADJUVANT PROSTATE
elif contains_string "$neo_prostate" "$input_file"; then
    echo "...chain.sh: detected neoadjuvant prostate project"
    RAD=${PROJ_HOME}/met/radial.csv
    TUMOR_MARKER=64
    roimag=1
    TCELL_MARKER=1024
    BCELL_MARKER=2097152
    IMMUNE_MARKER=34361852928 
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
	flipcmd="cyftools flip ${IO} -y 0 -Y $ymax ${DELIM}"
        echo "Found in LUT: $ymax"
    fi
else
    echo "Warning: No LSP pattern found in the string."
fi

# Make the magnif
y command
if [[ $input_file =~ (LSP[0-9]+) ]]; then
    PATTERN=${BASH_REMATCH[1]}
    
    # Lookup the pattern in the LUT
    if [[ -n "${maglut[$PATTERN]}" ]]; then
        mag=${maglut[$PATTERN]}
	magcmd="cyftools magnify ${IO} -f $mag ${DELIM}"
        echo "Found in magnify LUT: $mag"
    fi
else
    echo "Warning: No LSP pattern found in the string."
fi


## set the roi file
if [[ -f "$roi_file" ]]; then
    roicmd="cyftools roi ${IO} -m $roimag -r $roi_file ${DELIM}"
fi

## make the tls commands
if [[ -n "$BCELL_MARKER" && -n "$IMMUNE_MARKER" ]]; then
    tlscmd="cyftools tls ${IO} -b $BCELL_MARKER -i $IMMUNE_MARKER -m 300 -d 35 ${V} ${DELIM}"
fi

## make the radial density command
if [ -f "$RAD" ]; then
    radcmd="cyftools radialdens ${IO} -t ${T} -f ${RAD} ${V} ${DELIM}"
fi

##########
# RUN THE CHAIN OF CYFTOOLS COMMANDS
##########

if [[ ! -f "$input_file" ]]; then
    echo "Error in chain.sh: File '$input_file' does not exist."
    exit 1
else
    echo "...running: cyftools chain on ${base}"
cmd="cyftools check $input_file ${OUT} ${DELIM}
${magcmd:+$magcmd}
cyftools pheno ${IO} -t $pheno_file ${DELIM}
cyftools filter ${IO} -a $TUMOR_MARKER ${V} -M ${DELIM}
cyftools annotate ${IO} -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 ${DELIM}
${flipcmd:+$flipcmd}
${roicmd:+$roicmd}
${tlscmd:+$tlscmd}
${distann:+$distann}
${radcmd:+$radcmd}
cyftools delaunay ${IN} ${output_file} -l 20"

cmd_noisland="cyftools check $input_file ${OUT} ${DELIM}
${magcmd:+$magcmd}
cyftools pheno ${IO} -t $pheno_file ${DELIM}
cyftools filter ${IO} -a $TUMOR_MARKER ${V} -M ${DELIM}
cyftools annotate ${IO} -f 0.33 -k 25 -d 10000 -t ${T} ${V} -F 1 ${DELIM}
${flipcmd:+$flipcmd}
${roicmd:+$roicmd}
${distann:+$distann}
${radcmd:+$radcmd}
cyftools delaunay ${IN} ${output_file} -l 20"

short_cmd="cyftools check $input_file ${OUT} ${DELIM}
${magcmd:+$magcmd}
cyftools pheno ${IO} -t $pheno_file ${DELIM}
cyftools filter ${IO} -a $TUMOR_MARKER ${V} -M ${DELIM}
${roicmd:+$roicmd}
${tlscmd:+$tlscmd}
cyftools margin ${IO} ${V} -T 8 -M 64 -d 200 ${DELIM} cyftools margin ${IN} ${V} -T 128 -M 256 -d 200 ${output_file}"

echo "$cmd" | tr '\n' ' '
echo ""
eval "$cmd"

echo ${output_file}

fi

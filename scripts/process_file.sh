#!/bin/bash

# Check if USE_SLURM is set and non-zero
USE_SLURM=0

# Turn slurm off if on mac
if [ -d "/Users/" ]; then
    USE_SLURM=0
fi    

if [[ "${USE_SLURM:-0}" -eq 1 ]]; then
    echo "Submitting jobs to SLURM..."
else
    echo "Running jobs locally..."
fi

## configuration file
#### NEED TO EDIT TO LINK PROJECT NAME (e.g. JHU) to SAMPLES
source ${HOME}/git/cyftools/scripts/config.sh

############
### SET PROJECT SPECIFIC DIRECTORY STRUCTURE
############

#### ORION CRC
#HOMEBASE=${PROJ_DATA}/orion/orion_1_74
#ROIBASE=rois
#GATEBASE=phenotype
#CHAINDIR=chain_coy2

#### PROSTATE CYCIF
#PROJ_HOME=/Users/jeremiahwala/Dropbox_HMS/Primary_prostate_CyCIF_TME_manuscript_2023/primary_data
#PROJ_DATA=/Users/jeremiahwala/Dropbox_HMS/Primary_prostate_CyCIF_TME_manuscript_2023/primary_data
#HOMEBASE=${PROJ_DATA}
#ROIBASE=roisgu
#GATEBASE=gates_j
#CHAINDIR=chain

#### JHU ORION
#HOMEBASE=${PROJ_DATA}/jhu
#ROIBASE=roi
#GATEBASE=gates
#CHAINDIR=chain

#### JHU CYCIF
#HOMEBASE=${PROJ_DATA}/jhu/cycif
#ROIBASE=roi
#GATEBASE=phenotype
#CHAINDIR=chain

#### JHU CYCIF REVISION
#HOMEBASE=${PROJ_DATA}/jhu/revision
#ROIBASE=roi
#GATEBASE=phenotype
#CHAINDIR=chain

## NEOADJUVANT PROSTATE
HOMEBASE=${PROJ_DATA}/met
ROIBASE=${PROJ_DATA}=roi
GATEBASE=phenotype
CHAINDIR=chain

mkdir -p ${HOMEBASE}/${CHAINDIR}

echo "...getting file list from $HOMEBASE/clean"
for infile in $HOMEBASE/clean/*.cyf; do    

    ## uncomment to run just the one file
#if [[ ! "$base" =~ ^(LSP17711|LSP17708|LSP20071)$ ]]; then
#    continue
#fi    

    ## Extract the base LSP123 name
    base=$(basename "$infile" .cyf)
    echo $base
    
    echo "...process_file.sh: working on sample $base"
        
    ## Get the gates from the *p columns from csv's dumped from matlab files
    #if ! command -v Rscript &> /dev/null
    #then
    #    echo "Rscript could not be found"
    #    exit
    #fi
    #check_file_exists "$HOMEBASE/rawcsv/${base}.csv"
    #~/git/cysift/scripts/pheno.sh "$HOMEBASE/rawcsv/${base}.csv" "$HOMEBASE/$GATEBASE/${base}.csv"

    ## run the actual chain command
    if [[ "${USE_SLURM:-0}" -eq 1 ]]; then
        # Submit the job to SLURM
        sbatch ~/git/cyftools/scripts/chain.sh $HOMEBASE/clean/${base}.cyf\
					$HOMEBASE/$CHAINDIR/${base}.cyf\
					$HOMEBASE/$GATEBASE/${base}.phenotype.csv\
					${HOMEBASE}/${ROIBASE}/${base}.roi.csv
    else
        # Run interactively
        ~/git/cyftools/scripts/chain.sh $HOMEBASE/clean/${base}.cyf\
					$HOMEBASE/$CHAINDIR/${base}.cyf\
					$HOMEBASE/$GATEBASE/${base}.phenotype.csv\
					${HOMEBASE}/${ROIBASE}/${base}.roi.csv
    fi
    
    ## uncomment to run just one sample
    #exit 1
    
done

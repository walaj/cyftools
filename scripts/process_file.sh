#!/bin/bash

# Check if USE_SLURM is set and non-zero
USE_SLURM=1
if [[ "${USE_SLURM:-0}" -eq 1 ]]; then
    echo "Submitting jobs to SLURM..."
else
    echo "Running jobs interactively..."
fi

## configuration file
#### NEED TO EDIT TO LINK PROJECT NAME (e.g. JHU) to SAMPLES
source ${HOME}/git/cyftools/scripts/config.sh

############
### SET PROJECT SPECIFIC DIRECTORY STRUCTURE
############

#### ORION CRC
HOMEBASE=${PROJ_DATA}/orion/orion_1_74
ROIBASE=rois
GATEBASE=phenotype
CHAINDIR=chain_coy

#### PROSTATE CYCIF
#HOMEBASE=${PROJ_DATA}/prostate
#ROIBASE=roisgu
#GATEBASE=pheno
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
#CHAINDIR=chain_coy

echo "...getting file list from $HOMEBASE"
for infile in $HOMEBASE/clean/*.cyf; do    
    
    ## Extract the base LSP123 name
    base=$(basename "$infile" .cyf)
    echo "...process_file.sh: working on sample $base"
    
    ## uncomment to run just the one file
    ##if [[ $base != "LSP14483" ]]; then
    ##    continue
    ##fi
    
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

#!/bin/bash

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

#### PROSTATE CYCIF
#HOMEBASE=${PROJ_DATA}/prostate
#ROIBASE=roisgu
#GATEBASE=pheno

#### JHU ORION
#HOMEBASE=${PROJ_DATA}/jhu
#ROIBASE=roi
#GATEBASE=gates

#### JHU CYCIF
#HOMEBASE=${PROJ_DATA}/jhu/cycif
#ROIBASE=roi
#GATEBASE=phenotype

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
    ~/git/cyftools/scripts/chain.sh $HOMEBASE/clean/${base}.cyf\
				    $HOMEBASE/chain_coy/${base}.cyf\
				    $HOMEBASE/$GATEBASE/${base}.phenotype.csv\
				    ${HOMEBASE}/${ROIBASE}/${base}.roi.csv

    ## uncomment to run just one sample
    #exit 1
    
done

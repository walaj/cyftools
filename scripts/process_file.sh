#!/bin/bash

## set to 1 if want to only print run lines (create "parallel" compatible output)
#PCO=1

source ${HOME}/git/cyftools/scripts/config.sh

## download the data
#orion data
#wget https://www.dropbox.com/s/c88pb7ih3y6d94v/matlab-Orion_CRC_allbacthes-20220602.mat?dl=1 -O ~/projects/orion/matlab-Orion_CRC_allbacthes-20220602.mat

# prepare the matlab file
#matlab -nodisplay -r "run('/home/jaw34/git/cysift/matlab/jerry.m'); exit;"

#HOMEBASE=${PROJ_DATA}/orion/orion_1_74
#ROIBASE=rois

HOMEBASE=${PROJ_DATA}/prostate
ROIBASE=roisgu

parallel_echo "...getting file list from $HOMEBASE"
for infile in $HOMEBASE/rawcsv/*.csv; do

    if [[ ! $infile =~ rar ]]; then

	base=$(basename "$infile" .csv)
	
	# If the stripped name ends with .rar, skip to the next iteration
	if [[ $base == *.rar ]]; then
            continue
	fi

	# Extract the basename
	parallel_echo "...process_file.sh: working on sample $base"

	## uncomment to run just the one file
	##if [[ $base != "LSP14483" ]]; then
	##    continue
	##fi
	
	#check_file_exists $infile
	#~/git/cyftools/scripts/csv_rearrange.sh $infile "${HOMEBASE}/rawcsv/${base}.rar.csv"

	## Get the gates from the *p columns from csv's dumped from matlab files
	#if ! command -v Rscript &> /dev/null
	#then
	#    echo "Rscript could not be found"
	#    exit
	#fi
	#check_file_exists "$HOMEBASE/rawcsv/${base}.rar.csv"
	#~/git/cysift/scripts/pheno.sh $infile "$HOMEBASE/phenotype/${base}.phenotype.csv"
	
	## Put the cysift headers onto the csv files
	#check_file_exists "$HOMEBASE/phenotype/${base}.phenotype.csv"
	#~/git/cyftools/scripts/header.sh "$HOMEBASE/rawcsv/${base}.rar.csv" "$HOMEBASE/header/${base}.header.csv"

	## Convert the cysift csv files to cys files
	#check_file_exists "$HOMEBASE/header/${base}.header.csv"
	#~/git/cyftools/scripts/cerealed.sh "$HOMEBASE/header/${base}.header.csv" "$HOMEBASE/clean/${base}.cyf" 2>/dev/null

	#check_file_exists "$HOMEBASE/clean/${base}.cyf"
	~/git/cyftools/scripts/chain.sh "$HOMEBASE/clean/${base}.cyf" "$HOMEBASE/chain_flip/${base}.p.cyf" "$HOMEBASE/phenotype/${base}.phenotype.csv" "${HOMEBASE}/${ROIBASE}/${base}.roi.csv"
	exit 1
	#check_file_exists "${HOMEBASE}/chain/${base}.p.cyf"
	#sbatch ~/git/cyftools/scripts/margin_noisland.sh "${HOMEBASE}/chain/${base}.p.cyf" "${HOMEBASE}/margin_noisland/${base}.p.cyf"
    fi
done

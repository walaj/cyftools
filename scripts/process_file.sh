#!/bin/bash

source ~/git/cysift/scripts/config.sh


## download the data
#orion data
#wget https://www.dropbox.com/s/c88pb7ih3y6d94v/matlab-Orion_CRC_allbacthes-20220602.mat?dl=1 -O ~/projects/orion/matlab-Orion_CRC_allbacthes-20220602.mat

# prepare the matlab file
#matlab -nodisplay -r "run('/home/jaw34/git/cysift/matlab/jerry.m'); exit;"

#HOMEBASE=/n/scratch3/users/j/jaw34/projects/prostate/
HOMEBASE=/n/scratch3/users/j/jaw34/projects/orion/orion_1_74
HOMEBASE=/Users/jeremiahwala/Sorger/projects/orion/orion_1_74
for infile in $HOMEBASE/rawcsv/*.csv; do

    if [[ ! $infile =~ rar ]]; then

	base=$(basename "$infile" .csv)
	
	# If the stripped name ends with .rar, skip to the next iteration
	if [[ $base == *.rar ]]; then
            continue
	fi
	
	# Extract the basename
	echo "...process_file.sh: working on sample $base"

	#check_file_exists $infile
	#~/git/cysift/scripts/csv_rearrange.sh $infile "${base}.rar.csv"

	## Get the gates from the *p columns from csv's dumped from matlab files
	#if ! command -v Rscript &> /dev/null
	#then
	#    echo "Rscript could not be found"
	#    exit
	#fi
	#check_file_exists "$HOMEBASE/rawcsv/${base}.rar.csv"
	#~/git/cysift/scripts/pheno.sh $infile "$HOMEBASE/pheno/${base}.phenotype.csv"
	
	## Put the cysift headers onto the csv files
	#check_file_exists "$HOMEBASE/pheno/${base}.phenotype.csv"
	#~/git/cysift/scripts/header.sh "$HOMEBASE/rawcsv/${base}.rar.csv" "$HOMEBASE/header/${base}.header.csv"

	## Convert the cysift csv files to cys files
	#check_file_exists "$HOMEBASE/header/${base}.header.csv"
	#~/git/cysift/scripts/cerealed.sh "$HOMEBASE/header/${base}.header.csv" "$HOMEBASE/clean/${base}.cys" 2>/dev/null

	check_file_exists "$HOMEBASE/clean/${base}.cys"
	sbatch ~/git/cysift/scripts/chain.sh "$HOMEBASE/clean/${base}.cys" "$HOMEBASE/chain/${base}.ptrdim.cys" "$HOMEBASE/pheno/${base}.phenotype.csv" "${HOMEBASE}/roi/${base}.roi.csv"

	#check_file_exists "${HOMEBASE}/chain/${base}.ptrd.cys"
	#sbatch ~/git/cysift/scripts/margin_noisland.sh "${HOMEBASE}/chain/${base}.ptrd.cys" "${HOMEBASE}/margin_noisland/${base}.ptrdim.cys"
	
    fi
done

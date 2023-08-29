#!/bin/bash

if ! command -v Rscript &> /dev/null
then
    echo "Rscript could not be found"
    exit
fi

function check_file_exists {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "File not found: $file"
        exit 1
    fi
}

## download the data
#orion data
#wget https://www.dropbox.com/s/c88pb7ih3y6d94v/matlab-Orion_CRC_allbacthes-20220602.mat?dl=1 -O ~/projects/orion/matlab-Orion_CRC_allbacthes-20220602.mat

# prepare the matlab file
#matlab -nodisplay -r "run('/home/jaw34/git/cysift/matlab/jerry.m'); exit;"

HOMEBASE=/n/scratch3/users/j/jaw34/projects/orion/orion_1_74
for infile in $HOMEBASE/rawcsv/*.csv; do

    if [[ ! $infile =~ rar ]]; then

	#base="${infile%%.*}"
	base=$(basename "$infile" .csv)
	
	# Extract the basename
	echo "...working on sample $base"

	# check_file_exists $infile
	# ~/git/cysift/scripts/csv_rearrange.sh $infile "${base}.rar.csv"

	# check_file_exists "${base}.rar.csv"
	#~/git/cysift/scripts/pheno.sh $infile "${base}.phenotype.csv"

	# check_file_exists "${base}.phenotype.csv"
	#~/git/cysift/scripts/header.sh "${base}.rar.csv" "${base}.header.csv"

	#check_file_exists "$HOMEBASE/headered/${base}.header.csv"
	#~/git/cysift/scripts/cerealed.sh "$HOMEBASE/headered/${base}.header.csv" "$HOMEBASE/clean/${base}.cys" 2>/dev/null

	check_file_exists "$HOMEBASE/clean/${base}.cys"
	sbatch ~/git/cysift/scripts/chain.sh "$HOMEBASE/clean/${base}.cys" "$HOMEBASE/chain/${base}.ptrd.cys" "$HOMEBASE/phenotype/${base}.phenotype.csv"
	
	#check_file_exists "${base}.cys"
	#~/git/cysift/scripts/phenotype.sh "${base}.cys" "${base}.phenotype.cys" "${base}.phenotype.csv"

	#check_file_exists "${base}.phenotype.cys"
	#sbatch /home/jaw34/git/cysift/scripts/spatial.sh "${base}.phenotype.cys" "${base}.spat.cyz"

	#check_file_exists "${base}.phenotype.cys"
	#sbatch /home/jaw34/git/cysift/scripts/radial.sh "${base}.phenotype.cys" "${base}.rad.cys"

	#check_file_exists "${base}.rad.cys"
	#sbatch /home/jaw34/git/cysift/scripts/tumor.sh "${base}.rad.cys" "${base}.tumor.cys"

	#check_file_exists "${base}.tumor.cys"
	#cysift select -o 131072 -O 1 ${base}.tumor.cys ${base}.tselect.cys
	
	   #check_file_exists "${base}.tselect.cys"
	   #/home/jaw34/git/cysift/scripts/average.sh "${base}.tselect.cys" "${base}.avg.cys"
	
	   #check_file_exists "${base}.tselect.cys"
	   #sbatch /home/jaw34/git/cysift/scripts/frame.sh "${base}.tselect.cys" "${base}.tframe.csv"
	
    fi
done

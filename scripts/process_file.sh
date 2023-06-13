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

for infile in *.csv; do

    if [[ $(echo "$infile" | grep -o "\." | wc -l) -eq 1 ]]; then

	base="${infile%%.*}"
	
	# Extract the basename
	echo "...working on sample $base"

	# check_file_exists $infile
	# ~/git/cysift/scripts/csv_rearrange.sh $infile "${base}.rar.csv"

	# check_file_exists "${base}.rar.csv"
	#~/git/cysift/scripts/pheno.sh $infile "${base}.phenotype.csv"

	# check_file_exists "${base}.phenotype.csv"
	#~/git/cysift/scripts/header.sh "${base}.rar.csv" "${base}.header.csv"

	#check_file_exists "${base}.header.csv"
	#~/git/cysift/scripts/cerealed.sh "${base}.header.csv" "${base}.cys" 2>/dev/null

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
	
	   check_file_exists "${base}.tselect.cys"
	   /home/jaw34/git/cysift/scripts/average.sh "${base}.tselect.cys" "${base}.avg.cys"
	
	   #check_file_exists "${base}.tselect.cys"
	   #sbatch /home/jaw34/git/cysift/scripts/frame.sh "${base}.tselect.cys" "${base}.tframe.csv"
	
    fi
    
done

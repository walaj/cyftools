#!/bin/bash




for infile in *.csv; do

    if [[ $(echo "$infile" | grep -o "\." | wc -l) -eq 1 ]]; then

	base="${infile%%.*}"
	
	# Extract the basename
	echo "...working on sample $base"

	~/git/cysift/scripts/csv_rearrange.sh $infile "${base}.rar.csv"
	~/git/cysift/scripts/pheno.sh $infile "${base}.phenotype.csv"
	~/git/cysift/scripts/header.sh "${base}.rar.csv" "${base}.header.csv"
	~/git/cysift/scripts/cerealed.sh "${base}.header.csv" "${base}.cys" 2>/dev/null
	~/git/cysift/scripts/phenotype.sh "${base}.cys" "${base}.phenotype.cys" "${base}.phenotype.csv"
	sbatch /home/jaw34/git/cysift/scripts/spatial.sh "${base}.phenotype.cys" "${base}.spat.cyz"
    fi
    
done

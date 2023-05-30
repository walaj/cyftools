#!/bin/bash

O73="false"

#######
# orion 41-74
#######
if [[ "${O73}" == "true" ]]; then
    for file in data*.csv; do
	cp -- $file "${file#data}"
    done
    exit
fi


#######
# orion 1-40
#######

# read CSV file line by line
while IFS=, read -r specimen_id crc_id htma402id age gender primary histology grade location tnm stage ifm1 ifm2 ifm3 ifm4 lvi pni deposits border til mmrihc distant_mets recurrence location_of_recurrence death pfs_censor os_censor pfs_days os_days renato_tmb op_tmb hypermutant kras p53 nras braf pik3ca apc fbxw7 pten tcf7l2 sox9 ctnnb1 slide_id htan_participant_id cohort
do
    # skip the header line
    if [[ $specimen_id != "Specimen_ID" ]]; then
        # form old file name
        old_file_name="data${specimen_id}.csv"
        # form new file name
        new_file_name="${slide_id}.csv"
        # rename if old file exists
        if [[ -f $old_file_name ]]; then

	    echo "cp -- $old_file_name ${old_file_name#data}"
	    
	    # Change '-' to '_' in the first line of the file and write it to a temporary file
	    head -n 1 "$old_file_name" | sed 'y/-/_/' > tmp
	    
	    # Append the rest of the original file (from line 2 to the end) to the temporary file
	    tail -n +2 "$old_file_name" >> tmp

            echo "cp $old_file_name $new_file_name"
            cp tmp "$new_file_name"
	    rm tmp
	    
	else
	    echo "file not found $old_file_name"
        fi
    fi
done < ~/projects/orion/74OrionCasesDeidentifiedMasterList2023_scOPupdate_better.csv  

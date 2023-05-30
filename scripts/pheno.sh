#!/bin/bash

input_file=$1
out_file=$2

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    Rscript /home/jaw34/git/cysift/R/pheno_form.R $input_file $out_file
fi




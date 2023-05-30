#!/bin/bash

input_file=$1
out_file=$2
pheno_file=$3

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
fi

if [[ ! -f "$pheno_file" ]]; then
    echo "Error: Phenotype csv file $pheno_file does not exist"
    exit 1
fi

echo "...working on cysift pheno $input_file $output_file -t $pheno_file"
cysift pheno $input_file $out_file -t $pheno_file

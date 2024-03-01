#!/bin/bash

source ~/git/cysift/scripts/config.sh

input_file=$1
output_file=$2

base=$(basename "$input_file" .rar.csv)

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
elif [[ $orion41_73 == *"$base"* ]]; then
    echo "header.sh: detected Orion 41-73"
    header="${PROJ_HOME}/orion/header.txt"
elif [[ $orion1_40 == *"$base"* ]]; then
    echo "header.sh: detected Orion 1-40"
    header="${PROJ_HOME}/orion/header.txt"
elif [[ "$input_file" == *"immune"* ]]; then
    echo "header.sh: detected CyCIF Immune"
    header="${PROJ_HOME}/orion/header.immune.txt"
elif [[ "$input_file" == *"tumor"* ]]; then
    echo "header.sh: detected CyCIF Tumor"
    header="${PROJ_HOME}/orion/header.tumor.txt"
elif contains_string "$prostate" "$base"; then    
    echo "header.sh: detected Prostate"
    header="${PROJ_HOME}/prostate/header.txt"
else
    echo "header.sh: Warning: $input_file doesn't fit into cycif, prostate, orion, etc. Base: $base"
    exit 1
fi

# check that the header exists
if [[ ! -f "${header}" ]]; then
    echo "Error: Header file ${header} does not exist"
    exit 1
fi

# re-header the file
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "...re-headering on $input_file"
    tail -n +2 $input_file > tmpfile
    cat $header tmpfile > $output_file
    rm tmpfile
fi

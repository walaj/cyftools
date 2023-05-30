#!/bin/bash

input_file=$1
output_file=$2
header="/home/jaw34/projects/orion/header.txt"

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

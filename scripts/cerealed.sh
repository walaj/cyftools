#!/bin/bash

input_file=$1
output_file=$2

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "...running cereal on $input_file"
    cysift cereal ${input_file} ${output_file}
fi

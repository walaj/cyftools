#!/bin/bash

input_file=$1
output_file=$2

number="${input_file##*LSP}"
number="${number%%.*}"
echo "...cerealed.sh: detected sample number: $number"

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "cysift convert ${input_file} ${output_file} -s $number"
    cysift convert ${input_file} ${output_file} -s $number
fi

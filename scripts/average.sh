#!/bin/bash

input_file=$1
out_file=$2

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    echo "...cysift average $input_file $out_file"
    cysift average $input_file $out_file
fi




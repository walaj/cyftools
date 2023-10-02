#!/bin/bash

orion1_40="LSP10353,LSP10408,LSP10463,LSP10518,LSP10573,LSP10628,LSP10683,LSP10738,\
LSP10364,LSP10419,LSP10474,LSP10529,LSP10584,LSP10639,LSP10696,LSP10749,\
LSP10375,LSP10431,LSP10485,LSP10540,LSP10595,LSP10650,LSP10705,LSP10760,\
LSP10388,LSP10441,LSP10496,LSP10551,LSP10606,LSP10661,LSP10716,LSP10771,\
LSP10397,LSP10452,LSP10507,LSP10562,LSP10617,LSP10672,LSP10727,LSP10786"

orion41_73="LSP14383,LSP14388,LSP14363,LSP14373,LSP14393,LSP14398,LSP14438,\
LSP14468,LSP14503,LSP15304,LSP14403,LSP15280,LSP15308,LSP14408,LSP14443,LSP14473,LSP15284,\
LSP15312,LSP14413,LSP14448,LSP14483,LSP15288,LSP15316,LSP14418,LSP14453,LSP14493,LSP15292,\
LSP15320,LSP14423,LSP14458,LSP14498,LSP15300,LSP15324,LSP14463"

prostate="LSP12601,LSP12603,LSP12605,LSP12607,LSP12609,LSP12611,LSP12613,LSP12615,\
LSP12617,LSP12619,LSP12621,LSP12623,LSP12625,LSP12627,LSP12629,LSP12631,\
LSP12633,LSP12635,LSP12637,LSP12639,LSP12641,LSP12643,LSP12645,LSP12647,\
LSP12649,LSP12651,LSP12653,LSP12655,LSP12657"

input_file=$1
output_file=$2

contains_string() {
    local list="$1"
    local string_to_check="$2"
    local match_found=0

    IFS=',' read -ra ADDR <<< "$list"
    for i in "${ADDR[@]}"; do
        if [[ $string_to_check == *"$i"* ]]; then
            match_found=0
            break
        fi
    done

    return $match_found
}

base=$(basename "$input_file" .csv)

if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
elif [[ $orion41_73 == *"$base"* ]]; then
    echo "header.sh: detected Orion 41-73"
    header="/home/jaw34/projects/orion/header.txt"
elif [[ $orion1_40 == *"$base"* ]]; then
    echo "header.sh: detected Orion 1-40"
    header="/home/jaw34/projects/orion/header.txt"
elif [[ "$input_file" == *"immune"* ]]; then
    echo "header.sh: detected CyCIF Immune"
    header="/home/jaw34/projects/orion/header.immune.txt"
elif [[ "$input_file" == *"tumor"* ]]; then
    echo "header.sh: detected CyCIF Tumor"
    header="/home/jaw34/projects/orion/header.tumor.txt"
elif contains_string "$prostate" "$base"; then    
    echo "header.sh: detected Prostate"
    header="/home/jaw34/projects/prostate/header.txt"
else
    echo "header.sh: Warning: $input_file doesn't fit into cycif, prostate, orion, etc"
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

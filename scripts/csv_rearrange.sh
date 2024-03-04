#!/bin/bash

input_file=$1
rar_file=$2

orion1_40="LSP10353,LSP10408,LSP10463,LSP10518,LSP10573,LSP10628,LSP10683,LSP10738,\
LSP10364,LSP10419,LSP10474,LSP10529,LSP10584,LSP10639,LSP10696,LSP10749,\
LSP10375,LSP10431,LSP10485,LSP10540,LSP10595,LSP10650,LSP10705,LSP10760,\
LSP10388,LSP10441,LSP10496,LSP10551,LSP10606,LSP10661,LSP10716,LSP10771,\
LSP10397,LSP10452,LSP10507,LSP10562,LSP10617,LSP10672,LSP10727,LSP10786"

orion41_73="LSP14383,LSP14388,LSP14363,LSP14373,LSP14393,LSP14398,LSP14438,\
LSP14468,LSP14503,LSP15304,LSP14403,LSP15280,LSP15308,LSP14408,LSP14443,LSP14473,LSP15284,\
LSP15312,LSP14413,LSP14448,LSP14483,LSP15288,LSP15316,LSP14418,LSP14453,LSP14493,LSP15292,\
LSP15320,LSP14423,LSP14458,LSP14498,LSP15300,LSP15324,LSP14463"

prostate="dataLSP12601.csv,dataLSP12603.csv,dataLSP12605.csv,dataLSP12607.csv,dataLSP12609.csv,dataLSP12611.csv,dataLSP12613.csv,dataLSP12615.csv,\
dataLSP12617.csv,dataLSP12619.csv,dataLSP12621.csv,dataLSP12623.csv,dataLSP12625.csv,dataLSP12627.csv,dataLSP12629.csv,dataLSP12631.csv,\
dataLSP12633.csv,dataLSP12635.csv,dataLSP12637.csv,dataLSP12639.csv,dataLSP12641.csv,dataLSP12643.csv,dataLSP12645.csv,dataLSP12647.csv,\
dataLSP12649.csv,dataLSP12651.csv,dataLSP12653.csv,dataLSP12655.csv,dataLSP12657.csv"

base=$(basename "$input_file" .csv)

## install with "pip install csvkit"
## then arrange the columns below how you want in the output
if [[ ! -f "$input_file" ]]; then
    echo "...csv_rearrange.sh: Error: File '$input_file' does not exist."
    exit 1
elif [[ $orion41_73 == *"$base"* ]]; then
    echo "...csv_rearrange.sh: Orion 41-73 detected"
    csvcut -c "Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Blank,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD3e,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,Region" ${input_file} > ${rar_file}
elif [[ $orion1_40 == *"$base"* ]]; then
    echo "...csv_rearrange.sh: Orion 1-40 detected"
    csvcut -c "Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Argo550,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD3e,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,Region" ${input_file} > ${rar_file}
elif [[ "$input_file" == *"immune"* ]]; then
    echo "...csv_rearrange.sh: CyCIF immune detected"
    csvcut -c "Xt,Yt,Hoechst1,A488,A555,A647,CD3,CD11c,GranzymeB,Ki67,panCK,CD45,CD11b,CD68,CD14,CD163,FOXP3,CD8a,CD15,CD44,PD_L1,CD4,p_TBK1,PD_1,CD57,,Ki67,AMCAR,HMWCK,CD19,SMA,CD20,CD11b,CD68,CD163,CD4,CD3d,CD8a,TCF1,FOXP3,PD1,CD57,CD11c,GranzymeB,CD15,HLADR,CD103,CD31,pTBK1,HLAA,CD24,CD44,CD206" ${input_file} > ${rar_file}
elif [[ $prostate == *"$base"* ]]; then
    echo "...csv_rearrange.sh: Prostate detected"
    csvcut -c "Xt,Yt,Hoechst1,Ki67,AMCAR,HMWCK,CD19,SMA,CD20,CD11b,CD68,CD163,CD4,CD3d,CD8a,TCF1,FOXP3,PD1,CD57,CD11c,GranzymeB,CD15,HLADR,CD103,CD31,pTBK1,HLAA,CD24,CD44,CD206" $input_file > $rar_file
else
    echo "Warning: Unknown file"
fi

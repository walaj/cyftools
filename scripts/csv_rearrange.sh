#!/bin/bash

input_file=$1
rar_file=$2

## install with "pip install csvkit"
## then arrange the columns below how you want in the output
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else
    csvcut -c "CellID,Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Blank,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation,X,Y,Pan_CKp,CD31p,Ki67p,CD68p,CD163p,CD20p,CD4p,CD8ap,CD45ROp,PD_L1p,CD3ep,E_cadherinp,PD_1p,FOXP3p,CD45p,SMAp,col,row,frame,TDist,Region,topics" ${input_file} > ${rar_file}
fi




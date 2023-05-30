#!/bin/bash

input_file=$1
rar_file=$2

orion41_73="LSP14383,LSP14388,LSP14363,LSP14373,LSP14393,LSP14398,LSP14438"
orion41_73="${orion41_73},LSP14468,LSP14503,LSP15304,LSP14403,LSP15280"
orion41_73="${orion41_73},LSP15308,LSP14408,LSP14443,LSP14473,LSP15284"
orion41_73="${orion41_73},LSP15312,LSP14413,LSP14448,LSP14483,LSP15288"
orion41_73="${orion41_73},LSP15316,LSP14418,LSP14453,LSP14493,LSP15292"
orion41_73="${orion41_73},LSP15320,LSP14423,LSP14458,LSP14498,LSP15300,LSP15324,LSP14463"

## install with "pip install csvkit"
## then arrange the columns below how you want in the output
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
else

    base="${input_file%%.*}"
    
    if [[ $orion41_73 == *"$base"* ]]; then
	echo "Orion 41-73"
	csvcut -c "CellID,Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Blank,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD3e,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation,X,Y,Pan_CKp,CD31p,Ki67p,CD68p,CD163p,CD20p,CD4p,CD8ap,CD45ROp,PD_L1p,CD3ep,E_cadherinp,PD_1p,FOXP3p,CD45p,SMAp,col,row,frame,TDist,Region,topics" ${input_file} > ${rar_file}
    else
	echo "Orion 1-40"
	csvcut -c "CellID,Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Argo550,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation" ${input_file} > ${rar_file}
    fi

fi

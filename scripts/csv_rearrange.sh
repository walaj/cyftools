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

base=$(basename "$input_file" .csv)

## install with "pip install csvkit"
## then arrange the columns below how you want in the output
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' does not exist."
    exit 1
elif [[ $orion41_73 == *"$base"* ]]; then
    echo "Orion 41-73"
    csvcut -c "Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Blank,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD3e,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation,X,Y,Pan_CKp,CD31p,Ki67p,CD68p,CD163p,CD20p,CD4p,CD8ap,CD45ROp,PD_L1p,CD3ep,E_cadherinp,PD_1p,FOXP3p,CD45p,SMAp,col,row,frame,TDist,Region,topics" ${input_file} > ${rar_file}
elif [[ $orion1_40 == *"$base"* ]]; then
    echo "Orion 1-40"
    csvcut -c "Xt,Yt,Hoechst,AF1,CD31,CD45,CD68,Argo550,CD4,FOXP3,CD8a,CD45RO,CD20,PD_L1,CD3e,CD163,E_cadherin,PD_1,Ki67,Pan_CK,SMA,Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation" ${input_file} > ${rar_file}
elif [[ "$input_file" == *"immune"* ]]; then
    echo "CyCIF immune"
    csvcut -c "Xt,Yt,Hoechst1,A488,A555,A647,CD3,CD11c,GranzymeB,Ki67,panCK,CD45,CD11b,CD68,CD14,CD163,FOXP3,CD8a,CD15,CD44,PD_L1,CD4,p_TBK1,PD_1,CD57,HLA_DR,LAG3,CD20,TIM3,HLA_A,CD16,pSRC,pSTAT1,CD24,Desmin,PDPN,STING,pRB,pTyr" ${input_file} > ${rar_file}
elif [[ "$input_file" == *"tumor"* ]]; then
    echo "CyCIF Tumor"
    csvcut -c "Xt,Yt,Hoechst1,A488,A555,A647,NaKATPase,MCM6,p53,PCNA,panCK,SMA_1,E_Cad,SMA_2,Vimentin,EGFR,CDX2,CD45,pERK,H3K27me3,H3K27ac,beta_Catenin,CD44,PD_L1,N_Cadherin,H2ax,ZEB1,Nestin,CD68,CD31,S100a,GFAP,IRF_1,pS6_235,cPARP,p21,Var49,KAP1,Ki67,FN3,VEGFR2,NGFR" ${input_file} > ${rar_file}
else
    echo "Warning: Unknown file"
fi

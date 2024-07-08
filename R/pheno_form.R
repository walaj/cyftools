library(data.table)

options(scipen = 999)

contains_any <- function(substrings, target_string) {
  results <- sapply(substrings, function(sub) grepl(sub, target_string, fixed = TRUE))
  any(results)
}

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# check if command line arguments are provided
if (length(args) < 2) {
  stop("Two arguments are required: input and output file paths.")
}

# check if input file exists and is readable
if (!file.exists(args[1]) || !file.access(args[1], 4) == 0) {
  stop(paste("Input file does not exist or is not readable:", args[1]))
}

# check if output directory is writable
output_directory <- dirname(args[2])
if (!dir.exists(output_directory) || !file.access(output_directory, 2) == 0) {
  stop(paste("Output directory does not exist or is not writable:", output_directory))
}

# read input file
dt <- data.table::fread(args[1])

jhu_cycif <- c("LSP19083", "LSP19118", "LSP19193", "LSP19233", "LSP19273", "LSP19329", "LSP19378", "LSP20006", "LSP20051", "LSP20071", "LSP20091", "LSP20126", "LSP20146", "LSP19088", "LSP19158", "LSP19203", "LSP19238", "LSP19298", "LSP19339", "LSP19383", "LSP20041", "LSP20056", "LSP20081", "LSP20096", "LSP20136", "LSP20171", "LSP19093", "LSP19163", "LSP19228", "LSP19243", "LSP19324", "LSP19358", "LSP20001", "LSP20046", "LSP20066", "LSP20086", "LSP20121", "LSP20141")

# get column namesv
cols <- colnames(dt)

##
jhu_cycif <- c("LSP19083", "LSP19118", "LSP19193", "LSP19233", "LSP19273", "LSP19329", "LSP19378", "LSP20006", "LSP20051", "LSP20071", "LSP20091", "LSP20126", "LSP20146",                                                                    
               "LSP19088", "LSP19158", "LSP19203", "LSP19238", "LSP19298", "LSP19339", "LSP19383", "LSP20041", "LSP20056", "LSP20081", "LSP20096", "LSP20136", "LSP20171",                                                                    
               "LSP19093", "LSP19163", "LSP19228", "LSP19243", "LSP19324", "LSP19358", "LSP20001", "LSP20046", "LSP20066", "LSP20086", "LSP20121", "LSP20141") 

## link the file name to the project
if (grepl("immune",args[1])) {
    marker_cols <- c("panCKp","PDPNp","Desminp","Ki67p","CD11bp","CD11cp","CD14p","CD15p","CD16p","CD163p","CD20p","CD24p","CD4p","CD44p","CD45p","CD57p","CD68p","CD8ap","FOXP3p","GranzymeBp","HLA_Ap","HLA_DRp","LAG3p","p_TBK1p","PD_1p","PD_L1p","pSRCp","pSTAT1p","pTyrp","STINGp","pRBp")
    cat("...pheno_form.R: cycif crc immune\n")        
} else if (grepl("tumor",args[1])) {
    marker_cols <- c("panCKp","SMA_1p","ZEB1p","Vimentinp","E_Cadp","IRF_1p","Ki67p","CD31p","CD44p","CD45p","CD68p","N_Cadherinp","Nestinp","NGFRp","GFAPp","PD_L1p","H2axp","cPARPp")
    cat("...pheno_form.R: cycif crc tumor\n")    
} else if (grepl("LSP126", args[1])) { # prostate
    cat("...pheno_form.R: prostate\n")
    marker_cols <- c("AMCARp","HMWCKp","SMAp","CD20p","CD68p","CD163p","CD4p","CD3dp","CD8ap","FOXP3p","PD1p","CD57p","CD11cp","CD15p","HLADRp","CD103p","CD31p","pTBK1p","HLAAp","CD44p","CD206p")
} else if (contain_any(jhu_cycif, args[1]))
    cat("...pheno_form.R: jhu cycif\n")
    marker_cols <- c("panCKp","pAKTp","LINE1p","pERKp","pS6_240p","SMAp","Cateninp","CDX2p","gH2axp","pS6_235p","PCNAp","MLH1p","HER2p","Ki67p","H3K27me3p","TROP2p","CD44p","pNDRG1p","pRBp","p27p")
} else {
    cat("...pheno_form.R: cycif orion\n")    
    marker_cols <- grep("p$", cols, value = TRUE)
}

stopifnot(all(marker_cols %in% cols))

# remove trailing p
marker_cols <- sub("p$", "", marker_cols)
stopifnot(all(marker_cols %in% cols))


# for each marker, find the min value where corresponding p column equals 1
result <- lapply(marker_cols, function(marker) {
    p_col <- paste0(marker, "p")
    min_val <- dt[get(p_col) == 1, min(get(marker), na.rm=TRUE)]
  return(data.table::data.table(marker, min_val, 100000))
})

# convert the result to a data frame and write to output file
result_df <- data.table::rbindlist(result)
result_df[, marker := gsub("-","_",marker)]
write.table(result_df, args[2], row.names = FALSE,
            quote = FALSE, col.names = FALSE, sep = ",")

cat("...pheno_form.R successful\n")

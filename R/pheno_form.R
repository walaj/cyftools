library(data.table)

options(scipen = 999)

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

# get column namesv
cols <- colnames(dt)

# find the columns that match the pattern "AF\d+"
marker_cols <- grep("p$", cols, value = TRUE)
marker_cols <- sub("p$", "", marker_cols)

# for each marker, find the min value where corresponding p column equals 1
result <- lapply(marker_cols, function(marker) {
    p_col <- paste0(marker, "p")
    min_val <- dt[get(p_col) == 1, min(get(marker), na.rm=TRUE)]
  return(data.table::data.table(marker, min_val, 100000))
})

# convert the result to a data frame and write to output file
result_df <- data.table::rbindlist(result)
result_df[, marker := gsub("_","-",marker)]
write.table(result_df, args[2], row.names = FALSE,
            quote = FALSE, col.names = FALSE, sep = ",")

library(data.table)

options(scipen = 999)

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

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

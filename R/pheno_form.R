# load necessary libraries
library(dplyr)
library(tidyr)

options(scipen = 999)

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# read input file
df <- read.csv(args[1])

# get column names
cols <- names(df)

# find the columns that match the pattern "AF\d+"
marker_cols <- grep("\\d+$", cols, value = TRUE)

# for each marker, find the min value where corresponding p column equals 1
result <- lapply(marker_cols, function(marker) {
  p_col <- paste0(marker, "p")
  min_val <- df[df[[p_col]] == 1, marker] %>% min(na.rm = TRUE)
  return(c(marker, min_val, 100000))
})

# convert the result to a data frame and write to output file
result_df <- do.call(rbind, result)
colnames(result_df) <- c("Marker", "MinVal", "FixedVal")
write.table(result_df, args[2], row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ",")

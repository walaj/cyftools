library(data.table)

options(scipen=999)

DEBUG=FALSE
CELL_COUNT <- 10
radcols=c("CD31_100r","CD45_100r","CD68_100r","CD4_100r","FOXP3_100r","CD8_100r","CD20_100r","PD_L1_100r","CD3_100r","CD163_100r","Ecad_100r","PD1_100r","PanCK_100r","SMA_100r")
frame_size <- 300

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
cysfile <- args[1] 

# extract the sapmle
split_str <- strsplit(basename(cysfile), split = "\\.")
sample <- split_str[[1]][1]

if (DEBUG) cysfile <- "~/Sorger/orion/orion_1_40/LSP10353.rad.cys"

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

##########
### DATA READ
##########

# header read
dth <- fread(cmd=paste("cysift view -H", cysfile, "| grep '^@\\(MA\\|CA\\)'"), header=FALSE)
dth[, V1 := gsub("^@","",V1)]
dth[, V2 := gsub("^ID:","",V2)]
dth[, V3 := NULL]
dth[, V4 := NULL]
setnames(dth, c("V1","V2"),c("tag","id"))

# data read
cat("...reading data\n")
dt  <- fread(cmd=paste("cysift view", cysfile),header=FALSE)
setnames(dt, paste0("V",seq(ncol(dt))), c("cellid","cell_flag","pheno_flag","x","y",dth$id))
#radcols <- dth$id[seq(from=which(dth$id=="Orientation")+1, to=ncol(dt)-4)]
dt_rad <- dt[rowSums(dt != 0) > 0]

##########
### FRAMES
##########
cat("...setting up frames\n")
# frame function
get_frame <- function(coordinate) {
  return(floor(coordinate / frame_size) + 1)
}

# setup the frames
dt[, c("frame_x", "frame_y") := .(get_frame(x), get_frame(y))]
dt[, frame_id := frame_x * 1e6 + frame_y]

# count the cells per frame
dt.n <- dt[, .N, by=frame_id]
setnames(dt.n, "N","cellcount")

# find the frame identity of each cell
dt.framed <- dt[, lapply(.SD, mean), by = frame_id, .SDcols = radcols]

# round the counts
for (col in radcols) {
    dt.framed[[col]] <- round(dt.framed[[col]])
}

# remove rows
dt.framed <- dt.framed[rowSums(dt.framed[,..radcols] != 0) > 0]
dt.framed[, sample := sample]

# what are the frames themselves
frames <- dt[, .(centroid_x = (as.numeric(frame_x) - 0.5) * frame_size, 
                 centroid_y = (as.numeric(frame_y) - 0.5) * frame_size), 
             by = frame_id]
frames <- frames[!duplicated(frame_id)]
dt.framed <- merge(dt.framed, frames, all.x=TRUE, by="frame_id")
dt.framed <- merge(dt.framed, dt.n, all.x=TRUE, by="frame_id")
setcolorder(dt.framed, c("frame_id", "centroid_x", "centroid_y", "sample", "cellcount", setdiff(names(dt.framed), c("frame_id", "centroid_x", "centroid_y", "sample","cellcount"))))

write.table(dt.framed[cellcount > CELL_COUNT], args[2], row.names = FALSE,
            quote = FALSE, col.names = FALSE, sep = ",")

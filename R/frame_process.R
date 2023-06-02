library(data.table)

DEBUG=FALSE
radcols=c("CD31_200r","CD45_200r","CD68_200r","CD4_200r","FOXP3_200r","CD8_200r","CD20_200r","PD_L1_200r","CD3_200r","CD163_200r","Ecad_200r","PD1_200r","PanCK_200r","SMA_200r")
frame_size <- 400

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
setnames(dt, paste0("V",seq(ncol(dt))), c("cellid","flag","x","y",dth$id))
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

# find the frame identity of each cell
dt.framed <- dt[, lapply(.SD, mean), by = frame_id, .SDcols = radcols]
dt.framed <- dt.framed[, lapply(.SD, round), .SDcols = radcols]
dt.framed <- dt.framed[rowSums(dt.framed[,..radcols] != 0) > 0]
dt.framed[, sample := sample]

# what are the frames themselves
frames <- dt[, .(centroid_x = (as.numeric(frame_x) - 0.5) * frame_size, 
                 centroid_y = (as.numeric(frame_y) - 0.5) * frame_size), 
             by = frame_id]
frames <- frames[!duplicated(frame_id)]
dt.framed <- merge(dt.framed, frames, all.x=TRUE, by="frame_id")

setcolorder(dt.framed, c("frame_id", "centroid_x", "centroid_y", "sample", setdiff(names(dt.framed), c("frame_id", "centroid_x", "centroid_y", "sample"))))


write.table(dt.framed, args[2], row.names = FALSE,
            quote = FALSE, col.names = FALSE, sep = ",")
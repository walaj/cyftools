library(data.table)
library(ldatuning)
library(topicmodels)
library(parallel)
source("~/Sorger/orion/orion_helpers.R")

xfile <- file.path("~/Sorger/orion/orion_1_74/ptrdvl/all74.2M.cys")

# load the header
hd <- header_read(xfile)

# load the data
dt <- data.table::fread(cmd=paste("cysift subsample", xfile, "-r 0.1 - | cysift view -"))
colnames(dt) <- c("frameid","x","y","sample","count",hd[tag_type %in% c("MA","CA")]$id)
radcols=c("CD31_200r","CD45_200r","CD68_200r","CD4_200r","FOXP3_200r","CD8_200r","CD20_200r","PD_L1_200r","CD3_200r","CD163_200r","Ecad_200r","PD1_200r","PanCK_200r","SMA_200r")
dtm <- round(as(as.matrix(dt[sample(.N, 1e3),..radcols]), "CsparseMatrix"))

# remove empty rows
dtm <- dtm[apply(dtm, 1, function(row) sum(row) != 0), ]

system.time(model <- topicmodels::LDA(dtm, k = 8, method = "VEM", 
                                      control = list(seed = 42, verbose=1)))
system.time(model <- topicmodels::CTM(dtm, k = 3, method = "VEM", 
                                      control = list(seed = 42, verbose=1)))

system.time(result <- FindTopicsNumber(
 dtm,
  topics = seq(from = 8, to = 12, by = 1),
  metrics = "CaoJuan2009", # "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 831),
  mc.cores = 1L,
  verbose = TRUE))
FindTopicsNumber_plot(result)

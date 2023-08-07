library(data.table)
library(ldatuning)
library(topicmodels)
library(parallel)
library(Matrix)
source("~/Sorger/orion/orion_helpers.R")

xfile <- file.path("~/Sorger/orion/orion_1_74/ptrdvl/all74.2M.cys")

# load the header
hd <- header_read(xfile)

# load the data
dt <- data.table::fread(cmd=paste("cysift subsample", xfile, "-r 0.1 - | cysift view -"))
colnames(dt) <- c("cellid","pflag","cflag","x","y",hd[tag_type %in% c("MA","CA")]$id)
radcols=c("CD31_200r","CD45_200r","CD68_200r","CD4_200r","FOXP3_200r","CD8_200r","CD20_200r","PD_L1_200r","CD3_200r","CD163_200r","Ecad_200r","PD1_200r","PanCK_200r","SMA_200r")
dt[, (radcols) := lapply(.SD, function(x) x/mean(x)), .SDcols = radcols]

# allocate subsample tables
set.seed(42)
powers <- seq(2,5)
dtm <- vector("list", length(powers))
models <- vector("list", length(powers))
for (i in seq_along(powers)) {
    dtm[[i]] <- round(as(as.matrix(dt[sample(.N, 10^powers[i]),..radcols]), "CsparseMatrix"))
    dtm[[i]] <- dtm[[i]][apply(dtm[[i]], 1, function(row) sum(row) != 0), ] # remove empty rows
}

# scramble and allocate
dtscramble = data.table::copy(dt)
dtscramble[, (radcols) := lapply(.SD, sample), .SDcols = radcols]

# test training corpus size
for (i in seq_along(dtm)) {
    cat("working on ", powers[i], "\n")
    system.time(models[[i]] <- topicmodels::LDA(dtm[[i]], k = 10, method = "VEM", 
                                                control = list(seed = 42, verbose=0)))
}

# plot model comparison
markers_name <- sub("_200r$", "", radcols)
beta_dt <- vector("list", length(models))
for (i in seq_along(models))  {
    beta_dt[[i]] <- as.data.table(t(exp(models[[i]]@beta)))
    setnames(beta_dt[[i]], colnames(beta_dt[[i]]), paste("Topic", seq(ncol(beta_dt[[i]]))))
    beta_dt[[i]][, marker :=  markers_name]
    beta_dt[[i]][, model := powers[i]]
}
beta_long <- as.data.table(melt(rbindlist(beta_dt), id.vars = c("marker","model"), variable.name = "Topic", value.name = "Distribution"))
beta_long$Topic <- factor(beta_long$Topic, levels = unique(beta_long$Topic))
beta_long[, model := factor(model)]

# Create bar plot
g <- ggplot(beta_long, aes(x = factor(marker), y = Distribution, fill=model)) +
  #geom_bar(stat = "identity", position = "dodge", fill="gray30", color="black") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank()) +
  facet_wrap(~ Topic, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        legend.text = element_text(size = 6)) 

pdf("~/Sorger/figs/lda_iter_10.pdf", width=10, height=6, useDingbats = FALSE)
print(g)
dev.off()


#####
                          # SCRAMBLE
#####

                                        # scramble and allocate
set.seed(72)
dtscramble = data.table::copy(dt)
dtscramble[, (radcols) := lapply(.SD, sample), .SDcols = radcols]

# allocate subsample tables
dtm <- vector("list", length(powers))
models <- vector("list", length(powers))
for (i in seq(4)) {
    dtm[[i]] <- round(as(as.matrix(dtscramble[sample(.N, 10^4),..radcols]), "CsparseMatrix")*1000)
    dtm[[i]] <- dtm[[i]][apply(dtm[[i]], 1, function(row) sum(row) != 0), ] # remove empty rows
}

# test training corpus size
for (i in seq_along(dtm)) {
    cat("working on ", i, "\n")
    system.time(models[[i]] <- topicmodels::LDA(dtm[[i]], k = 10, method = "VEM", 
                                                control = list(seed = 42, verbose=0)))
}

# plot model comparison
markers_name <- sub("_200r$", "", radcols)
beta_dt <- vector("list", length(models))
for (i in seq_along(models))  {
    beta_dt[[i]] <- as.data.table(t(exp(models[[i]]@beta)))
    setnames(beta_dt[[i]], colnames(beta_dt[[i]]), paste("Topic", seq(ncol(beta_dt[[i]]))))
    beta_dt[[i]][, marker :=  markers_name]
    beta_dt[[i]][, model := i]
}
beta_long <- as.data.table(melt(rbindlist(beta_dt), id.vars = c("marker","model"), variable.name = "Topic", value.name = "Distribution"))
beta_long$Topic <- factor(beta_long$Topic, levels = unique(beta_long$Topic))
beta_long[, model := factor(model)]

# Create bar plot
g <- ggplot(beta_long, aes(x = factor(marker), y = Distribution, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank()) +
  facet_wrap(~ Topic, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        legend.text = element_text(size = 6)) 

pdf("~/Sorger/figs/lda_iter_scramble_10.pdf", width=10, height=6, useDingbats = FALSE)
print(g)
dev.off()


# test training corpus size
for (i in seq_along(dtm)) {
    cat("working on ", powers[i], "\n")

    system.time(models[[i]] <- topicmodels::LDA(dtm[[i]], k = 10, method = "VEM", 
                                                control = list(seed = 42, verbose=0)))
}



system.time(model <- topicmodels::CTM(dtm, k = 4, method = "VEM", 
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

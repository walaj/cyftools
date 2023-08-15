library(data.table)
library(jsonlite)
library(reshape2)
library(ggplot2)

model_file <- "~/Desktop/lda7M.200r.50its.12t.tumor.json"
model_file <- c("~/Sorger/orion/tmp/scramble.s0.l42.json",
                "~/Sorger/orion/tmp/scramble.s0.l100.json",
                "~/Sorger/orion/tmp/scramble.s0.l200.json")

beta_dt <- rbindlist(lapply(model_file, function(x) {
  # Replace "path_to_your_file.json" with the actual path to your JSON file
  json_content <- jsonlite::fromJSON(x)
  
  # Convert to data.table
  alpha_dt <- data.table(json_content$LDAModel$alpha)
  beta_dt <-
    data.table(t(json_content$LDAModel$beta))  # transpose to get one row per array
  setnames(beta_dt, colnames(beta_dt), paste("Topic", seq(ncol(beta_dt))))
  markers <- json_content$LDAModel$markers
  markers_name <- sub("_200r$", "", markers)
  beta_dt[, marker :=  markers_name]
  beta_dt[, model := basename(x)]
  return (beta_dt)
}))

##### PLOT
# Reshape data from wide to long format
beta_long <- as.data.table(melt(beta_dt, id.vars = c("marker","model"), variable.name = "Topic", value.name = "Distribution"))

# Convert Topic to factor for ordering in plot
beta_long$Topic <- factor(beta_long$Topic, levels = unique(beta_long$Topic))

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

pdf("~/Sorger/figs/lda_scramble_3seedsonscramble.pdf", width=8, height=6, useDingbats = FALSE)
print(g)
dev.off()


library(lda)
library(Matrix)
library(text2vec)
library(textTinyR)
library(lsa)
library(ggplot2)

num_topics <- c(5, 10, 15, 20)  # replace with your list of number of topics
num_keywords <- 15

# Assuming model_list is a list of lda models and models were trained on a DTM named dtm
LDA_models <- model_list
LDA_topics <- lapply(LDA_models, function(x) {
  top_words <- apply(x$beta, 1, function(y) order(y, decreasing = TRUE)[1:num_keywords])
  return(top_words)
})

jaccard_similarity <- function(topic_1, topic_2){
  intersection = length(intersect(topic_1, topic_2))
  union = length(union(topic_1, topic_2))
  return (intersection / union)
}

LDA_stability <- list()
for (i in 1:(length(num_topics) - 1)) {
  jaccard_sims <- list()
  for (topic1 in LDA_topics[[num_topics[i]]]) {
    sims <- list()
    for (topic2 in LDA_topics[[num_topics[i + 1]]]) {
      sims <- c(sims, jaccard_similarity(topic1, topic2))
    }
    jaccard_sims <- c(jaccard_sims, list(unlist(sims)))
  }
  LDA_stability[[num_topics[i]]] <- unlist(jaccard_sims)
}

mean_stabilities <- sapply(LDA_stability, mean)

# Coherence score calculation
coherences <- list()
for (num_topic in num_topics) {
  model <- LDA_models[[num_topic]]
  top_terms <- LDA_topics[[num_topic]]
  coherences[[num_topic]] <- coherence(x = model, dtm = dtm, top_terms = top_terms, measure = "UMass")
}

# Optimal number of topics
coh_sta_diffs <- mapply('-', coherences, mean_stabilities)
ideal_topic_num <- as.numeric(names(coh_sta_diffs)[which.max(coh_sta_diffs)])

# Plot
df <- data.frame(Number_of_Topics = num_topics, 
                 Mean_Stabilities = mean_stabilities, 
                 Coherences = coherences, 
                 stringsAsFactors = FALSE)

ggplot(df, aes(Number_of_Topics)) +
  geom_line(aes(y = Mean_Stabilities), color = 'blue') +
  geom_line(aes(y = Coherences), color = 'red') +
  geom_vline(xintercept = ideal_topic_num, color = "black", linetype = "dashed") +
  labs(x = "Number of Topics", y = "Metric Level",
       title = "Model Metrics per Number of Topics")

d
lapply(

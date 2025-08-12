library(tidyverse)
library(vegan)
library(corrplot)

setwd("")
data <- read.csv("KoSC.csv", stringsAsFactors = FALSE)
data[data == "Shared"] <- "Species"
data_clean <- data[!(apply(data, 1, function(row) all(row == "NoHit"))), ]
binary_data <- as.data.frame(lapply(data_clean, function(col) ifelse(col == "Species", 1, 0)))
binary_t <- t(binary_data)
jaccard_dist <- vegdist(binary_t, method = "jaccard")
jaccard_sim <- 1 - as.matrix(jaccard_dist)
jaccard_sim_full <- matrix(1, nrow = nrow(jaccard_sim), ncol = ncol(jaccard_sim))
rownames(jaccard_sim_full) <- rownames(jaccard_sim)
colnames(jaccard_sim_full) <- colnames(jaccard_sim)
jaccard_sim_full[lower.tri(jaccard_sim_full)] <- jaccard_sim[lower.tri(jaccard_sim)]
jaccard_sim_full[upper.tri(jaccard_sim_full)] <- jaccard_sim[upper.tri(jaccard_sim)]
pdf("Fig_4-A.pdf", width = 8, height = 5)
# Plot
corrplot(
  jaccard_sim_full,
  method  = "color",
  type    = "lower",
  is.corr = FALSE,                             
  col     = colorRampPalette(c("blue","white","red"))(200),
  col.lim = c(0, 1),                         
  cl.lim  = c(0, 1),                       
  addCoef.col = "black",
  tl.pos  = "n",
  cl.pos  = "r",
  cl.cex  = 0.8,
  number.cex = 2,
  tl.cex = 0.9,
  number.font = 1,
  mar = c(1, 3, 2, 1),
  addgrid.col = "black"
)
custom_labels <- c(
  expression("QIIME2"),
  expression("KoSCAPE"[bio]),
  expression("Metaphlan (1%)"),
  expression("SRST2(" * italic("bla"[OXY]) * ")")
)
n <- length(custom_labels)
x_pos <- 1:n
y_pos <- n:1
text(x = x_pos, y = n + 0.55, labels = custom_labels, xpd = TRUE, cex = 0.9)
text(x = 0.3, y = y_pos, labels = custom_labels, xpd = TRUE, cex = 0.9)
dev.off()
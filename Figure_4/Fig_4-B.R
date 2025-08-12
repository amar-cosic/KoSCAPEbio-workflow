library(tidyverse)
library(vegan)
library(corrplot)

setwd("")
data <- read.csv("kosc_plus.csv", stringsAsFactors = FALSE)
data <- data[!(apply(data, 1, function(row) all(row == "NoHit"))), ]
data <- data.frame(lapply(data, function(col) {
  ifelse(col == "NoHit", NA, trimws(tolower(col)))
}), stringsAsFactors = FALSE)
methods <- colnames(data)
jaccard_sim <- matrix(NA, nrow = length(methods), ncol = length(methods),
                      dimnames = list(methods, methods))
for (i in seq_along(methods)) {
  for (j in seq_along(methods)) {
    method_i <- data[[methods[i]]]
    method_j <- data[[methods[j]]]
    valid_rows <- !(is.na(method_i) & is.na(method_j))
    agree <- sum(method_i[valid_rows] == method_j[valid_rows], na.rm = TRUE)
    total <- sum(valid_rows)
    jaccard_sim[i, j] <- if (total == 0) NA else agree / total
  }
}

jaccard_sim_full <- jaccard_sim
diag(jaccard_sim_full) <- 1
col_pal <- colorRampPalette(c("blue", "white", "red"))(200)
pdf("Fig_4-B.pdf", width = 8, height = 5)
# --- plot ------------------------------------------------------------------
corrplot(
  jaccard_sim_full,
  method  = "color",
  type    = "lower",
  is.corr = FALSE,        
  col     = col_pal,
  cl.lim  = c(0, 1),       
  addCoef.col = "black",
  tl.pos  = "n",
  cl.pos  = "r",
  cl.cex  = 0.8,
  number.cex = 1.5,
  number.font = 1,  
  tl.cex = 0.9,
  mar = c(1, 3, 2, 1),
  addgrid.col = "black"
)
custom_labels <- c(
  expression("QIIME2"),
  expression("KoSCAPE"[bio]),
  expression("Metaphlan (1%)"),
  expression("SRST2 (" * italic("bla")[OXY] * ")")
)
n <- length(custom_labels)
text(1:n, y = n + 0.55, labels = custom_labels, xpd = TRUE, cex = 0.9)
text(x = 0.3, y = n:1, labels = custom_labels, xpd = TRUE, cex = 0.9)
dev.off()

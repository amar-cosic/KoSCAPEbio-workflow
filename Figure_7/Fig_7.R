library(tidyverse)
library(epitools)

setwd("")
df <- read.csv("raw_data.csv", stringsAsFactors=FALSE)
df$NEC <- factor(df$NEC)
df$Definition <- factor(df$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))
incidence_table <- df %>%
  group_by(Definition, NEC) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = NEC, values_from = n, values_fill = 0)

custom_colors <- c(
  "KoSC" = "tomato",
  "KpSC" = "mediumseagreen",
  "co-occurance" = "skyblue",
  "Empty" = "gray70"
)
counts_df <- df %>%
  group_by(NEC, Definition) %>%
  summarise(n = n(), .groups = "drop")
percent_df <- counts_df %>%
  group_by(NEC) %>%
  mutate(percent = n / sum(n) * 100)
raw_bar <- ggplot(df, aes(x = NEC, fill = Definition)) +
  geom_bar(position = "stack") +
  labs(title = "Raw Occurrences of Definition by NEC Status",
       y = "Count",
       x = "NEC") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()
raw_bar_numbers <- raw_bar +
  geom_text(
    data = counts_df,
    aes(x = NEC, y = n, label = n, group = Definition),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  )
norm_bar <- ggplot(df, aes(x = NEC, fill = Definition)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Relative Occurrence (Normalized to 100%) of Definition by NEC Status",
       y = "Percentage",
       x = "NEC") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()
norm_bar_numbers <- norm_bar +
  geom_text(
    data = percent_df,
    aes(x = NEC, y = percent/100, label = paste0(round(percent, 1), "%"), group = Definition),
    position = position_fill(vjust = 0.5),
    color = "black",
    size = 4
  )
outcomes <- levels(df$Definition)
or_results <- lapply(outcomes, function(def) {
  tbl <- table(df$NEC, df$Definition == def)
  res <- oddsratio.wald(tbl)
  pval <- fisher.test(tbl)$p.value
  data.frame(
    Definition = def,
    Odds_Ratio = res$measure[2,1],
    Lower_CI = res$measure[2,2],
    Upper_CI = res$measure[2,3],
    P_value = pval
  )
})
or_df <- do.call(rbind, or_results)
or_df$Definition <- factor(or_df$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))

forest_plot <- ggplot(or_df, aes(x = Odds_Ratio, y = Definition)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  labs(title = "Odds Ratio Forest Plot by Definition",
       x = "Odds Ratio (95% CI)",
       y = NULL) +
  theme_minimal()
summary_counts <- data.frame(
  Definition = c("co-occurance", "KoSC", "KpSC", "Empty"),
  NEC1 = c(5, 11, 5, 1),
  NEC0 = c(21, 14, 5, 5),
  stringsAsFactors = FALSE
)
summary_counts$Definition <- factor(summary_counts$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))




library(epitools)
or_counts_results <- lapply(1:nrow(summary_counts), function(i) {
  this_row <- summary_counts[i, ]
  nec1_in  <- this_row$NEC1
  nec1_out <- sum(summary_counts$NEC1) - nec1_in
  nec0_in  <- this_row$NEC0
  nec0_out <- sum(summary_counts$NEC0) - nec0_in
  tbl <- matrix(c(nec1_in, nec1_out, nec0_in, nec0_out), nrow = 2)
  res <- oddsratio.wald(tbl)
  pval <- fisher.test(tbl)$p.value
  data.frame(
    Definition = as.character(this_row$Definition),
    Odds_Ratio = res$measure[2, 1],
    Lower_CI = res$measure[2, 2],
    Upper_CI = res$measure[2, 3],
    P_value = pval
  )
})
or_counts_df <- do.call(rbind, or_counts_results)
or_counts_df$Definition <- factor(or_counts_df$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))
print(or_counts_df)


library(ggplot2)
forest_plot_counts <- ggplot(or_counts_df, aes(x = Odds_Ratio, y = Definition)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  labs(title = "Odds Ratio Forest Plot Sim",
       x = "Odds Ratio (95% CI)",
       y = NULL) +
  theme_minimal()
forest_plot_counts
library(dplyr)
raw_summary <- df %>%
  filter(Definition != "Empty") %>%   # <---- Exclude the "false" Empty rows
  group_by(Definition, NEC) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = NEC, values_from = n, values_fill = 0) %>%
  rename(NEC1 = `1`, NEC0 = `0`)
raw_summary$Definition <- as.character(raw_summary$Definition)
summary_counts$Definition <- as.character(summary_counts$Definition)
combined <- full_join(
  raw_summary, summary_counts,
  by = "Definition",
  suffix = c("_raw", "_summary")
)
combined[is.na(combined)] <- 0

combined <- combined %>%
  mutate(
    NEC1 = NEC1_raw + NEC1_summary,
    NEC0 = NEC0_raw + NEC0_summary
  ) %>%
  select(Definition, NEC1, NEC0)

combined$Definition <- factor(combined$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))
library(epitools)
library(ggplot2)

or_combined_results <- lapply(1:nrow(combined), function(i) {
  this_row <- combined[i, ]
  nec1_in  <- this_row$NEC1
  nec1_out <- sum(combined$NEC1) - nec1_in
  nec0_in  <- this_row$NEC0
  nec0_out <- sum(combined$NEC0) - nec0_in
  tbl <- matrix(c(nec1_in, nec1_out, nec0_in, nec0_out), nrow = 2)
  res <- oddsratio.wald(tbl)
  pval <- fisher.test(tbl)$p.value
  data.frame(
    Definition = as.character(this_row$Definition),
    Odds_Ratio = res$measure[2, 1],
    Lower_CI = res$measure[2, 2],
    Upper_CI = res$measure[2, 3],
    P_value = pval
  )
})
or_combined_df <- do.call(rbind, or_combined_results)
or_combined_df$Definition <- factor(or_combined_df$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))
print(or_combined_df)
or_combined_df_dummy <- do.call(rbind, or_combined_results)
or_combined_df_dummy$Definition <- factor(or_combined_df_dummy$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))
or_combined_df_dummy$Definition <- recode(or_combined_df_dummy$Definition, "Empty" = "Unclassified")
print(or_combined_df_dummy)
write.csv(or_combined_df_dummy, "or_combined_df_dummy.csv", row.names = FALSE)

or_combined_df_dummy$Definition <- factor(
  or_combined_df_dummy$Definition,
  levels = c( "Unclassified", "co-occurance", "KpSC","KoSC")
)
forest_colors <- c(
  "KoSC" = "#ff0000",
  "KpSC" = "#ffd36f",
  "co-occurance" = "#009e73",
  "Unclassified" = "gray70"
)
# Forest plot
forest_plot_combined <- ggplot(or_combined_df_dummy, aes(x = Odds_Ratio, y = Definition)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, color = "black") +
  geom_point(
    aes(fill = Definition),  
    color = "black",         
    size = 4,
    shape = 21,              
    stroke = 1
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(
    values = forest_colors,
    drop = FALSE,
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_x_log10() +
  labs(
    x = expression("Odds Ratio (95% CI, log"[10]*" scale)"),
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_blank()
  )

library(tidyverse)
combined$Definition <- factor(combined$Definition, levels = c("KoSC", "KpSC", "co-occurance", "Empty"))
combined_long <- combined %>%
  pivot_longer(cols = c("NEC1", "NEC0"), names_to = "NEC", values_to = "n") %>%
  mutate(NEC = recode(NEC, "NEC1" = "1", "NEC0" = "0"))

combined_long$Definition <- recode(combined_long$Definition, "Empty" = "Unclassified")
combined_long$Definition <- factor(
  combined_long$Definition,
  levels = c("KoSC", "KpSC", "co-occurance", "Unclassified")
)

custom_colors <- c(
  "KoSC" = "#e11c1c",           
  "KpSC" = "#f1cd94",           
  "co-occurance" = "#159a85",   
  "Unclassified" = "#8f8f8f"    
)

combined_long <- combined_long[!is.na(combined_long$Definition), ]
bar_combined_raw_numbers <- ggplot(combined_long, aes(x = NEC, y = n, fill = Definition)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    geom_text(
      aes(label = n),
      position = position_stack(vjust = 0.5),
      color = "black", size = 4
    ) +
    scale_fill_manual(values = custom_colors, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 200, by = 25)
    ) +
    scale_x_discrete(
      labels = c(
        "0" = "Control\n(n=175)",
        "1" = "NEC\n(n=45)"
      )
    ) +
    labs(
      y = "Number of samples",
      x = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_blank(),
      plot.title = element_blank()
    )
bar_combined_raw_numbers

combined_percent <- combined_long %>%
  group_by(NEC) %>%
  mutate(percent = n / sum(n) * 100)

combined_percent$Definition <- recode(combined_percent$Definition, "Empty" = "Unclassified")
combined_percent$Definition <- factor(
  combined_percent$Definition,
  levels = c("KoSC", "KpSC", "co-occurance", "Unclassified")
)
combined_percent <- combined_percent[!is.na(combined_percent$Definition), ]
bar_combined_norm_numbers <- ggplot(combined_percent, aes(x = NEC, y = percent, fill = Definition)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_text(
    aes(label = paste0(round(percent, 1), "%")),
    position = position_fill(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  scale_fill_manual(values = c(
    "KoSC" = "#ff0000",
    "KpSC" = "#ffd36f",
    "co-occurance" = "#009e73",
    "Unclassified" = "gray70"
  )) +
  scale_y_continuous(
    labels = function(x) scales::number_format(accuracy = 1)(x * 100),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_x_discrete(
    labels = c(
      "0" = "Control\n(n=175)",
      "1" = "NEC\n(n=45)"
    )
  ) +
  labs(
    y = "Distribution population subtype (%)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank()
  )
bar_combined_norm_numbers


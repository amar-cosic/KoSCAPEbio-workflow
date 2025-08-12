library(ape)
library(ggtree)
library(ggplot2)
library(stringr)
library(dplyr)

tree <- read.tree("tree_Shortest.nwk")

tip_labels <- tree$tip.label
tip_group <- case_when(
  str_starts(tip_labels, "Ko_") ~ "Ko",
  str_starts(tip_labels, "Kg_") ~ "Kg",
  str_starts(tip_labels, "Km_") ~ "Km",
  str_starts(tip_labels, "Kp_") ~ "Kp",
  str_starts(tip_labels, "Kh_") ~ "Kh",
  str_starts(tip_labels, "Ks_") ~ "Ks",
  TRUE ~ "Other"
)

tip_data <- data.frame(label = tip_labels, group = tip_group)

group_colors <- c(
  "Ko" = "red",
  "Kg" = "#56B4E9",   # Light blue
  "Km" = "limegreen",
  "Kp" = "purple",
  "Kh" = "black",
  "Ks" = "gold",
  "Other" = "gray"
)

p <- ggtree(tree, layout="daylight") %<+% tip_data +
  geom_tippoint(aes(color = group), size = 1.2) +
  scale_color_manual(values = group_colors) +
  theme(legend.position = "right") +
  theme_void() +
  theme(legend.text = element_text(size = 12), legend.title = element_blank())


p + geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = group_colors) +
  theme(legend.position = "right") +
  theme_void() +   # Clean background, no grid/axes
  geom_treescale(x = -7.5, y = 0.1, width = 0.02, offset = 0.001, fontsize = 4) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank())

p
p_flipped <- p + scale_x_reverse()

p_flipped_2 <- p_flipped +  scale_y_reverse()
p_flipped_2

ggsave("radial_tree_final.pdf", plot = p, width = 8, height = 8)
ggsave("radial_tree_final.png", plot = p, width = 8, height = 8, dpi = 600)

ggsave("radial_tree_final_flipped.pdf", plot = p_flipped_2, width = 8, height = 8)
ggsave("radial_tree_final_flipped.png", plot = p_flipped_2, width = 8, height = 8, dpi = 600)

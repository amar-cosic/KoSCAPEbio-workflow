# Load required libraries
library(ggplot2)
library(reshape2)
setwd("C:\\Users\\amar_\\OneDrive\\Desktop\\KoSCAPE\\GitHUB\\fig5")
a=1
if (a>2) {
  # Read the data
  main_pre_abs <- read.csv("main_pressence_absence.csv", header = TRUE, sep = ";", row.names = 1)
  
  # Convert '+' to 8 (Red) and '-' to 9 (Blue)
  main_pre_abs[main_pre_abs == "+"] <- 8
  main_pre_abs[main_pre_abs == "-"] <- 9
  
  # Ensure Risk probability stays numeric
  main_pre_abs["Risk probability", ] <- as.numeric(main_pre_abs["Risk probability", ])
  
  # Convert the entire data frame to numeric for plotting
  main_pre_abs_numeric <- as.data.frame(lapply(main_pre_abs, as.numeric))
  rownames(main_pre_abs_numeric) <- rownames(main_pre_abs)
  
  # Add row names as a separate column for 'melt'
  main_pre_abs_numeric$Species <- rownames(main_pre_abs_numeric)
  
  # Melt the data specifying 'Species' as the ID variable
  data_long <- melt(main_pre_abs_numeric, id.vars = "Species", variable.name = "Marker", value.name = "Value")
  
  # Create a color scale
  custom_colors <- c("9" = "blue", "8" = "red", "2" = "black", "1" = "gray", "0" = "white")
  
  # Define the exact order and set as factor levels
  species_order <- c("Klebsiella spp.",
                     "Klebsiella_oxytoca", "Klebsiella_grimontii", "Klebsiella_michiganensis",
                     "42_Shared_G/P/S (++)", "47_Shared_M/P (+)", "44_Shared_M/G/P_V2 (+)", "46_Shared_M/G (+)", "48_Shared_M/EC/SE (-)",
                     "Klebsiella_aerogenes",
                     "QIIME2 level 7 - Ko",
                     "Klebsiella spp. (rest)",
                     "Risk probability")
  
  species_order <- rev(species_order)
  # Update the Species column to the factor with the specified levels
  data_long$Species <- factor(data_long$Species, levels = species_order)
  
  # Ensure the Value is a factor for coloring
  data_long$Value <- as.factor(data_long$Value)
  
  # Plot the heatmap with borders around the tiles
  heatmap_pressence <- ggplot(data_long, aes(x = Marker, y = Species, fill = Value)) +
      geom_tile(colour = "black", size = 0.5) +  # White borders, adjust 'size' for thickness
      scale_fill_manual(values = custom_colors) +
      labs(title = "Presence/Absence Heatmap", x = "Markers", y = "Bacterial Species") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(color = "black"))
  ggsave("main_pressence.pdf", plot = heatmap_pressence, width = 25, height = 8, unit = "in")
}
# Load required libraries
library(ggplot2)
library(reshape2)

# Read the data
main_abundance <- read.csv("main_abundance.csv", header = TRUE, sep = ";", row.names = 1)

# Convert the data frame to numeric for plotting, except row names
main_abundance_numeric <- as.data.frame(lapply(main_abundance, as.numeric))
rownames(main_abundance_numeric) <- rownames(main_abundance)

# Add row names as a separate column for 'melt'
main_abundance_numeric$Species <- rownames(main_abundance_numeric)

# Melt the data specifying 'Species' as the ID variable
abundance_long <- melt(main_abundance_numeric, id.vars = "Species", variable.name = "Marker", value.name = "Value")

# Define the same species order as before and set as factor levels
species_order <- c("Klebsiella spp.",
                   "Klebsiella_oxytoca", "Klebsiella_grimontii", "Klebsiella_michiganensis",
                   "42_Shared_G/P/S (++)", "47_Shared_M/P (+)", "44_Shared_M/G/P_V2 (+)", "46_Shared_M/G (+)", "48_Shared_M/EC/SE (-)",
                   "Klebsiella_aerogenes",
                   "QIIME2 level 7 - Ko",
                   "Klebsiella spp. (rest)",
                   "Risk probability")

species_order <- rev(species_order)
abundance_long$Species <- factor(abundance_long$Species, levels = species_order)

# Define the color scale with enhanced contrast for low values
colors <- c("0" = "#5a8ac6", "low" = "#FFCC99", "high" = "red")
min_positive <- round(min(abundance_long$Value[abundance_long$Value > 0], na.rm = TRUE), digits=3)
max_value <- 100

# Adjust the values parameter to increase the visibility of low values
color_gradient <- scale_fill_gradientn(
  colors = colors,
  values = scales::rescale(c(0, 0.0001, 100)),  # Adjust this line to change color distribution
  breaks = c(0, min_positive, max_value),
  labels = c("0", as.character(min_positive), "100"),
  limits = c(0, max_value)
)

# Plot the heatmap with the adjusted color scale
heatmap_plot <- ggplot(abundance_long, aes(x = Marker, y = Species, fill = Value)) +
  geom_tile(color = "black", size = 0.5) +  # Correct use of 'size' for tile border width
  color_gradient +
  labs(title = "Relative Abundance Heatmap", x = "Markers", y = "Bacterial Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black"))
heatmap_plot
ggsave("main_abundance.pdf", plot = heatmap_plot, width = 25, height = 8, unit = "in")

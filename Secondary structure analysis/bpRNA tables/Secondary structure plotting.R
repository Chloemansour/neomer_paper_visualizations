install.packages("ggplot2")
install.packages("reshape2")
install.packages("pheatmap")
install.packages("patchwork")
install.packages("scales")
install.packages("RColorBrewer")
install.packages("gridExtra")

library(ggplot2)
library(reshape2)
library(pheatmap)
library(patchwork)
library(scales)
library(RColorBrewer)
library(gridExtra)

#################### Custom theme, palettes and functions ######################

# Define the theme for the plots
theme_custom <- function(base_size = 12) {
  theme_grey(base_size = base_size) + 
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Bold title
      plot.subtitle = element_text(size = 14, hjust = 0.5), # Subtitle
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_text(size = 12), # y-axis title
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), # X axis labels at 45 degrees
      axis.text.y = element_text(size = 10), # Text size for y-axis labels
      axis.line = element_line(size = 0.5, color = "black"), # Axis lines
      legend.position = "none", # No legend for individual plots
      
      panel.grid.major = element_line(color = "grey80"), # Lighter grid lines
      panel.grid.minor = element_blank(), # No minor grid lines
      panel.background = element_rect(fill = "white", color = NA), # White panel background
      plot.background = element_rect(fill = "white", color = NA), # White plot background
      strip.background = element_blank(), # No background for facet labels
      strip.text = element_text(face = "bold", size = 14) # Bold facet labels
    )
}

# Define the custom publication theme
publication_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.line = element_line(size = 1, color = "black"), # Darker and thicker axis lines
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    )
}

# Example of using a more professional color palette from RColorBrewer
palette_qualitative <- RColorBrewer::brewer.pal(8, "Set2")  # Qualitative palette for discrete data
palette_sequential <- RColorBrewer::brewer.pal(9, "Blues")   # Sequential palette for continuous data

# Adjust the plotting functions to use the new color palettes
plot_character_counts <- function(data, dataset_name, max_value) {
  gg <- ggplot(data, aes(x=Structure, y=Frequency, fill=Structure)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=palette_qualitative) +
    labs(y=paste("Structure Frequency -", dataset_name)) +
    ylim(0, max_value) + # Standardize the y-axis scale
    theme_custom()
  return(gg)
}

plot_avg_contiguous_lengths <- function(data, dataset_name, max_value) {
  gg <- ggplot(data, aes(x=Structure, y=AvgLength_Nucleotides, fill=Structure)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=palette_qualitative) +
    labs(y=paste("Nucleotide Length -", dataset_name)) +
    ylim(0, max_value) + # Standardize the y-axis scale
    theme_custom()
  return(gg)
}

# Adjust the plotting function for the heatmap
plot_heatmap_by_position <- function(data, dataset_name, title_label) {
  # Convert Positions to numeric for sorting
  data$Position <- as.numeric(data$Position)
  data <- data[order(data$Position),]
  
  # Melt the data for plotting
  data_melted <- melt(data, id.vars = 'Position')
  
  # Create the ggplot heatmap
  gg <- ggplot(data_melted, aes(x = Position, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "Reds")) +
    labs(title = title_label, y = dataset_name) +  # Set the y-axis label to the dataset name
    theme_custom() +
    theme(axis.title.x = element_text(size = 12), legend.position = "right") +
    guides(fill = guide_legend(title = "Relative Abundance"))
  
  return(gg)
}

#################### Read in data and format ######################

# Reading the average contiguous lengths data into variables
selex_lengths_data <- read.csv('selex_avg_contiguous_lengths_nucleotides.csv')[-8,]
selex_structural_data <- read.csv('selex_frequencies.csv')[-8,]
selex_heatmap_by_position_data <- read.csv('selex_heatmap_by_position.csv')[,c(-1,-9)]

neomer_lengths_data <- read.csv('neomer_avg_contiguous_lengths_nucleotides.csv')[-8,]
neomer_structural_data <- read.csv('neomer_frequencies.csv')[-8,]
neomer_heatmap_by_position_data <- read.csv('neomer_heatmap_by_position.csv')[,c(-1,-9)]

# Find the maximum values for setting the same scale
max_frequency <- max(c(selex_structural_data$Frequency, neomer_structural_data$Frequency))
max_length <- max(c(selex_lengths_data$AvgLength_Nucleotides, neomer_lengths_data$AvgLength_Nucleotides))

#Labels and position
colnames(selex_heatmap_by_position_data) = Labels
selex_heatmap_by_position_data['Position'] = rownames(selex_heatmap_by_position_data)

colnames(neomer_heatmap_by_position_data) = Labels
neomer_heatmap_by_position_data['Position'] = rownames(neomer_heatmap_by_position_data)

# Specify labels for the plot
Labels = c("Dangling end", "Stem", "Hairpin",
           "Internal loop", "Bulge", "External loop",
           "Multi loop")

selex_lengths_data$Structure <- Labels
neomer_lengths_data$Structure <- Labels
selex_structural_data$Structure <- Labels
neomer_structural_data$Structure <- Labels

# Assuming the Structure column is factor and ordering it for the plots
selex_lengths_data$Structure <- factor(selex_lengths_data$Structure, levels=Labels)
neomer_lengths_data$Structure <- factor(neomer_lengths_data$Structure, levels=Labels)

#################### Plotting - secondary structure analysis ######################

# Create the plots without printing them
p1 <- plot_character_counts(selex_structural_data, 'SELEX', max_frequency) + 
  labs(tag = "(A)") # Add tag for the plot letter
p2 <- plot_avg_contiguous_lengths(selex_lengths_data, 'SELEX', max_length) + 
  labs(tag = "(B)") # Add tag for the plot letter
p3 <- plot_character_counts(neomer_structural_data, 'Neomer', max_frequency) + 
  labs(tag = "(C)") # Add tag for the plot letter
p4 <- plot_avg_contiguous_lengths(neomer_lengths_data, 'Neomer', max_length) + 
  labs(tag = "(D)") # Add tag for the plot letter

# Combine the plots into a single panel with patchwork
plot_panel <- (p1 | p2) / 
  (p3 | p4) + 
  plot_layout(guides = "collect")

# Add the main title using plot_annotation
plot_panel <- plot_panel + 
  plot_annotation(
    title = "Secondary Structure Analysis",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
    )
  )

# Print the combined plot panel with the main title
plot_panel

#################### Plotting - secondary structure distribution analysis ######################



# Create the heatmaps with the updated function using the common color scale
p5 <- plot_heatmap_by_position(selex_heatmap_by_position_data, 'SELEX', NULL)
p6 <- plot_heatmap_by_position(neomer_heatmap_by_position_data, 'Neomer', NULL)

# Combine the heatmaps into a two-panel figure with a main title using gridExtra
final_plot <- grid.arrange(p5, p6, ncol = 1, 
                           top = textGrob("Secondary Structure Distribution Analysis", 
                                          gp = gpar(fontface = "bold", fontsize = 20)))

# Print the final plot
print(final_plot)


# Add the main title using plot_annotation
plot_panel <- plot_panel + 
  plot_annotation(
    title = "Secondary Structure Distribution Analysis",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
    )
  )

# Combine the heatmaps into a two-panel figure with a main title using gridExtra
final_plot <- grid.arrange(p5, p6, ncol=1, 
                           top=textGrob("Secondary Structure Distribution Analysis", 
                                        gp=gpar(fontface="bold", fontsize=20)))

# Print the final plot
print(final_plot)

#################### Plotting - diversity analysis ######################

# Read SDI data
sdi_data <- read.csv('sdi_bpRNA_data.csv')

# Plot SDI for SELEX and NEOMER with defined borders along the axes
sdi_plot <- ggplot(sdi_data, aes(x=Position)) +
  geom_line(aes(y=SELEX_SDI, colour="SELEX"), size=1, alpha = 0.5) +
  geom_line(aes(y=NEOMER_SDI, colour="NEOMER"), size=1, alpha = 0.5) +
  scale_colour_manual(values = c("SELEX" = "blue", "NEOMER" = "red")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.5) # Add border here
  ) +
  labs(
    title="Shannon Diversity Index across Positions",
    colour="Dataset",
    x="Position",
    y="Shannon Diversity Index"
  )

# Print the SDI plot
print(sdi_plot)

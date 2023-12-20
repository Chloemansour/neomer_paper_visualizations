# Install required libraries
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("grid")
install.packages("dplyr")

# Load the required libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# Read and prepare the data for Figure 2
Final_all <- read.delim('Final_All.high.txt', sep='\t')[c(1:10000),]
Final_all$Fold_change <- Final_all$Target.average/Final_all$Counter.target.average
Final_all <- Final_all %>% filter(!is.na(Z.score), !is.na(Fold_change))

Folds <- read.delim('Folds.txt', sep='\t')[c(1:10000),]
Full_folds <- read.delim('Folds_frequencies_IL6.txt', sep='\t')

# Define the custom theme for improved aesthetics for Figure 5
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

# Panel A
p1 <- ggplot(Final_all, aes(x = Z.score, y = Fold_change)) +
  geom_point(alpha = 0.6, size = 2, color = "black") +
  scale_x_continuous(limits = c(min(Final_all$Z.score, na.rm = TRUE), max(Final_all$Z.score, na.rm = TRUE))) +
  scale_y_continuous(limits = c(min(Final_all$Fold_change, na.rm = TRUE), max(Final_all$Fold_change, na.rm = TRUE))) +
  labs(title = "(A) IL-6 top 10,000 selected aptamers", x = "Z-score", y = "Fold change") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none")

# Panel D
p4 <- ggplot(Folds, aes(x = Fold.IL.6.Naive, y = Fold.IL.6.HSA)) +
  geom_point() +
  labs(title = "(D) IL-6/Naive vs IL-6/HSA", x = "IL-6/Naive", y = "IL-6/HSA") +
  custom_theme

# Panel E
p5 <- ggplot(Folds, aes(x = Fold.IL.6.UL, y = Fold.IL.6.HSA)) +
  geom_point() +
  labs(title = "(E) IL-6/Ultralink vs IL-6/HSA", x = "IL-6/Ultralink", y = "IL-6/HSA") +
  custom_theme

# Panel F
p6 <- ggplot(Folds, aes(x = Fold.IL.6.Nickel, y = Fold.IL.6.UL)) +
  geom_point() +
  labs(title = "(F) IL-6/Nickel resin-1 vs IL-6/Ultralink-1", x = "IL-6/Nickel resin-1", y = "IL-6/Ultralink-1") +
  custom_theme

# Panel B
p2 <- ggplot(Full_folds, aes(x = Frequency, fill = Condition)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "(B) Density Plot of Frequency", x = "Frequency", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Panel C 
p3 <- ggplot(Full_folds, aes(x = Fold_Change, fill = Condition)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "(C) Density Plot of Fold Change", x = "Fold Change", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine the plots into a 3x2 grid, aligning the titles and ensuring consistent aesthetics
combined_plot <- grid.arrange(p1, p4, p5, p6, p2, p3, ncol = 3, nrow = 2, 
                              top = textGrob("IL-6 Aptamer Analysis", 
                                             gp = gpar(fontface = "bold", fontsize = 20)))

# Save the combined plot with high resolution for publication
ggsave("publication_ready_plot.png", combined_plot, width = 18, height = 12, dpi = 300)


# Define a custom theme function to set dark axis lines and position the legend at the top right
publication_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size * 1.2, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size * 0.8),
      legend.position = "top",
      legend.justification = "right",
      legend.text = element_text(size = base_size * 0.8),
      legend.title = element_text(size = base_size),
      axis.line = element_line(size = 0.5, color = "black")
    )
}

# Apply the theme to all plots and adjust the legends
p1 <- p1 + publication_theme()
p2 <- p2 + publication_theme() + guides(fill=guide_legend(title="IL-6 vs HSA")) # Remove "Fold" from legend title
p3 <- p3 + publication_theme() + guides(fill=guide_legend(title="IL-6 vs Ultralink")) # Remove "Fold" from legend title
p4 <- p4 + publication_theme()
p5 <- p5 + publication_theme()
p6 <- p6 + publication_theme()

# Adjust axis labels for plots D, E, and F to Z-score comparisons
p4 <- p4 + labs(x = "Fold change IL-6 vs Naïve", y = "Fold change IL-6 vs HSA")
p5 <- p5 + labs(x = "Fold change IL-6 vs Ultralink", y = "Fold change IL-6 vs HSA")
p6 <- p6 + labs(x = "Fold change IL-6 vs Nickel resin", y = "Fold change IL-6 vs Ultralink")

# Remove subfigure titles and update the axis titles
p1 <- p1 + labs(title = "(A)", x = "Z-score", y = "Fold change")
p2 <- p2 + labs(title = "(B)", x = "Frequency")
p3 <- p3 + labs(title = "(C)", x = "Fold Change")
p4 <- p4 + labs(title = "(D)")
p5 <- p5 + labs(title = "(E)")
p6 <- p6 + labs(title = "(F)")

# Correct the legend position by specifying the right theme elements
correct_legend_theme <- function() {
  theme(
    legend.justification = c(1, 1), # Justify legend position at the top-right
    legend.position = c(1, 1) # Position the legend at the top-right corner
  )
}

# Modify the scale_fill_brewer function to have no title
p2 <- p2 + scale_fill_brewer(palette = "Set1", name = NULL)
p3 <- p3 + scale_fill_brewer(palette = "Set1", name = NULL)

# Ensure the guides function also specifies no title for the legend
p2 <- p2 + guides(fill=guide_legend(title=NULL))
p3 <- p3 + guides(fill=guide_legend(title=NULL))

# Apply the correct_legend_theme only to plots B and C
p2 <- p2 + correct_legend_theme()
p3 <- p3 + correct_legend_theme()

# Recombine the plots
combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, 
                              top = textGrob("IL-6 Aptamer Analysis", 
                                             gp = gpar(fontsize = 20, fontface = "bold")))

# Save the final combined plot
ggsave("publication_ready_plot.png", combined_plot, width = 18, height = 12, dpi = 300)
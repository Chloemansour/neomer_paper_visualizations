################### Figure 3 #################
install.packages("ggplot2")
install.packages("viridis")
install.packages("tidyr")

library(tidyr)
library(gridExtra)
library(ggplot2)
library(grid)
library(viridis)

setwd("~/NGS Related Files/Naive publication plots")

# Function to create a line plot
create_line_plot <- function(aptamer, column1, column2, column3, column4, aptamer_name, max_y = NULL) {
  # Reorder factor levels to control the order in the legend
  legend_order <- c("IL6 500 nM", "IL6 750 nM", "IL6 1000 nM", "HSA 1000 nM")
  
  # Create line plot
  plot <- ggplot(aptamer, aes(x = Time..sec.)) +
    geom_line(aes(y = !!column1, linetype = "IL6 500 nM", color = "IL6 500 nM"), size = 1) +
    geom_line(aes(y = !!column2, linetype = "IL6 750 nM", color = "IL6 750 nM"), size = 1) +
    geom_line(aes(y = !!column3, linetype = "IL6 1000 nM", color = "IL6 1000 nM"), size = 1) +
    geom_line(aes(y = !!column4, linetype = "HSA 1000 nM", color = "HSA 1000 nM"), size = 1) + 
    scale_linetype_manual(values = c("IL6 500 nM" = "dotted", 
                                     "IL6 750 nM" = "dashed", 
                                     "IL6 1000 nM" = "solid", 
                                     "HSA 1000 nM" = "solid"),
                          name = "Condition",
                          breaks = legend_order,
                          labels = c("IL6 500 nM", 
                                     "IL6 750 nM", 
                                     "IL6 1000 nM", 
                                     "HSA 1000 nM")) +
    scale_color_manual(values = c("IL6 500 nM" = "black", 
                                  "IL6 750 nM" = "black", 
                                  "IL6 1000 nM" = "black", 
                                  "HSA 1000 nM" = "maroon"),
                       name = "Condition",
                       breaks = legend_order,
                       labels = c("IL6 500 nM", "IL6 750 nM", "IL6 1000 nM", "HSA 1000 nM")) +
    labs(x = "Time (seconds)", y = paste("Binding Resonance:", aptamer_name, sep = " ")) +
    theme_minimal() + 
    theme(
      axis.line = element_line(color = "black", linewidth = 0.5), 
      axis.text = element_text(size = 10, color = "black"), 
      axis.title = element_text(size = 11.), 
      legend.position = "bottom", 
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) 
  # Set maximum y-axis value if specified
  if (!is.null(max_y)) {
    plot <- plot + ylim(-0.3, max_y)
  }
  
  return(plot)
}

# Read in binding assay values
IL6_4202_1 <- read.csv("IL6_4202_1_binding_assay.csv")
IL6_7326_1 <- read.csv("IL6_7326_1_binding_assay.csv")
IL6_6449 <- read.csv("IL6_6449_binding_assay.csv")
IL6_9805 <- read.csv("IL6_9805_binding_assay.csv")

# Plot binding resonance for each IL-6 candidate aptamer 
p1 <- create_line_plot(IL6_4202_1, quote(X500.nM.IL.6),quote(X750.nM.IL.6), quote(X1.uM.IL.6),quote(X1.uM.HSA), "IL6-4202.1", max_y = 1)
p2 <- create_line_plot(IL6_7326_1, quote(X500.nM.IL.6),quote(X750.nM.IL.6), quote(X1.uM.IL.6), quote(X1.uM.HSA), "IL6-7326.1", max_y = 1)
p3 <- create_line_plot(IL6_6449, quote(X500.nM.IL.6),quote(X750.nM.IL.6), quote(X1.uM.IL.6), quote(X1.uM.HSA), "IL6-6649", max_y = 1)
p4 <- create_line_plot(IL6_9805, quote(X500.nM.IL.6),quote(X750.nM.IL.6), quote(X1.uM.IL.6), quote(X1.uM.HSA), "IL6-9805", max_y = 1)

# Create 2x2 figure panel
combined_plot <- grid.arrange(p1, p3, p2, p4, ncol = 2, nrow = 2)

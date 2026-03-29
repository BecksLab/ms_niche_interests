library(ggplot2)
library(ggrepel)
library(ggtext)
library(patchwork)

## 2. General Figure Theme ----

figure_theme = 
  theme_classic() +
  theme(panel.border = element_rect(colour = '#00181F',
                                    fill = "white"), # Transparent fill
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::lighten("#00181F", 0.7),
                                  linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = "#00181F"), # Soft gray text for better readability
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

## 3. Model Color Palette ----
# Assigns specific hex codes to each modelling framework
# This ensures a model is always the same colour across all plots

# 1. DISCRETE MODEL COLORS
model_colours <- c(
  # Group 1: The Anchor
  "ATN"      = "#FFB81C",
  "LTM"      = "#FB8D6C", 
  "ADBM"     = "#EA7200", 
  "Random"   = "#002B36", 
  "Niche"    = "#7A74C2", 
  "Cascade"  = "#006D75",
  "MaxEnt"   = "#9CDBD9"
)

pal_df <- data.frame(l = names(model_colours), c = model_colours)

# 2. CONTINUOUS RAMP

col_cont <- c("#CCF5FF", "#83C5BE", "#002B36")

# 3. DIVERGING RAMP

col_div <- c("#006D75", "#F4F1DE", "#EA7200")
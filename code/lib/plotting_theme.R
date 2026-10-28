library(ggplot2)
library(ggrepel)
library(ggtext)
library(patchwork)

# for accents
shark_silver <- "#A2AAAD"
shark_black <- "#212322"
shark_white <- "#F4F1DE"

## 2. General Figure Theme ----

figure_theme = 
  theme_classic(16) +
  theme(panel.border = element_rect(colour = shark_black,
                                    fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::lighten(shark_black, 0.7),
                                  linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = shark_black),
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

# this set the theme once so dont need to call every time when plotting
set_theme(figure_theme)

## 3. Model Color Palette ----
# Assigns specific hex codes to each modelling framework
# This ensures a model is always the same colour across all plots

# 1. DISCRETE MODEL COLORS
model_colours <- c(
  "ATN"      = "#9CDBD9",
  "LTM"      = "#002B36",
  "ADBM"     = "#006D75",
  "Random"   = "#FFB81C",
  "Niche"    = "#E6BCD7",
  "Cascade"  = "#C085C2",
  "MaxEnt"   = "#EA7200"
)

pal_df <- data.frame(l = names(model_colours), c = model_colours)

# 2. CONTINUOUS RAMP

col_cont <- c("#CCF5FF", "#83C5BE", "#002B36")

# 3. DIVERGING RAMP

col_div <- c("#EA7200", shark_white, "#006D75")

# Lookup dictionaries for variables

# models

model_order <- c("ADBM", "ATN", "LTM", "Niche", "Cascade", "MaxEnt", "Random")
model_colours <- model_colours[model_order]

# model groups

model_groups <- c(
  ADBM = "Body-size",
  ATN = "Body-size",
  LTM = "Body-size",
  Niche = "Theoretical",
  Cascade = "Theoretical",
  MaxEnt = "Structural",
  Random = "Structural"
)

mod_shapes <- c(
  ADBM = 19,
  ATN = 19,
  LTM = 19,
  Niche = 17,
  Cascade = 17,
  MaxEnt = 15,
  Random = 15
)

# stablity

stability_labels <- c(
  persistence = "Species persistence",
  biomass_shannon = "Shannon diversity of biomasses",
  gini_fluxes = "Gini coeffient of consumption fluxes",
  skewness_IS = "Skewness of interaction strengths",
  resilience = "Resilience",
  reactivity = "Reactivity"
)

stability_groups <- c(
  persistence = "Network stability",
  biomass_shannon = "Network stability",
  gini_fluxes = "Network stability",
  skewness_IS = "Network stability",
  resilience = "Local stability",
  reactivity = "Local stability"
)


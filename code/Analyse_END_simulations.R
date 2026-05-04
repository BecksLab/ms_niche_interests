rm(list = ls())

library(tidyverse)
library(here)
library(MASS)

setwd(here())

source("lib/plotting_theme.R")

#------1. Load data ------------------------------------------------------------
dynamic_metric <- read_csv(
  "data/outputs/END_simulation_summary_03_05_2026.csv",
  show_col_types = FALSE)

dynamic_metric <- dynamic_metric %>%
  mutate(
    # use log10(x + 1) because some food webs collapse at equilibrium
    log_maxTL = log10(max_trophic_level + 1),
    log_total_biomass = log10(total_biomass + 1),
    model = factor(model),
    model = relevel(model, ref = "Niche")) # set niche model as reference

eq_topology <- read_csv(
  "data/outputs/END_topology_03_05_2026.csv",
  show_col_types = FALSE) %>%
  dplyr::select(-fw_ID) %>%
  mutate(
    model = factor(model),
    model = relevel(model, ref = "Niche"))

#------2. Difference checks across models --------------------------------------
# MANOVA for dynamic outputs
dyn_metrics <- c(
  "persistence",
  "log_maxTL",
  "log_total_biomass",
  "shannon",
  "evenness",
  "resilience",
  "reactivity"
)

dynamic_df <- dynamic_metric %>%
  select(model, all_of(dyn_metrics)) %>%
  na.omit()

dynamic_dep <- dynamic_df %>%
  select(all_of(dyn_metrics)) %>%
  as.matrix()

dynamic_manova <- manova(dynamic_dep ~ model, data = dynamic_df)
summary(dynamic_manova, test = "Pillai")


# MANOVA for equilibrium topology outputs
topology_metrics <- eq_topology %>%
  select(-model) %>%
  names()

topology_df <- eq_topology %>%
  select(model, all_of(topology_metrics)) %>%
  na.omit()

topology_dep <- topology_df %>%
  select(all_of(topology_metrics)) %>%
  as.matrix()

fit_topology <- manova(topology_dep ~ model, data = topology_df)
summary(fit_topology, test = "Pillai")

#------3. Boxplot of all simulation outputs across models ---------------------
dyn_metric_labels <- c(
  persistence       = "Species persistence",
  log_maxTL         = "Log10(max trophic level + 1)",
  log_total_biomass = "Log10(community biomass + 1)",
  shannon           = "Shannon diversity",
  evenness          = "Community evenness",
  resilience        = "Resilience",
  reactivity        = "Reactivity"
)

dynamic_out_long <- dynamic_df %>%
  mutate(across(any_of(dyn_metrics), ~ ifelse(is.nan(.x), NA_real_, .x))) %>%
  pivot_longer(
    cols = all_of(dyn_metrics), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    # Keep facet order following dyn_metrics
    metric = factor(metric, levels = dyn_metrics),
    
    # Rename metrics using the same vector
    metric_label = recode(as.character(metric), !!!dyn_metric_labels),
    
    # Keep renamed facet labels in the same order
    metric_label = factor(metric_label, levels = dyn_metric_labels[dyn_metrics])
  )

boxplot_sim <- ggplot(dynamic_out_long, aes(x = model, y = value, colour = model)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~ metric_label, scales = "free_y") +  
  labs(x = "", y = "") +
  scale_colour_manual(values = model_colours)  +
  figure_theme +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

boxplot_sim

ggsave("../figures/boxplots_dynamic_metrics.png", boxplot_sim, 
       width = 12, height = 7, dpi = 600)

#------4. LDA on dynamic metrics -----------------------------------------------
# Fit LDA using dynamic metrics
lda_dynamic_fit <- MASS::lda(model ~ ., data = dynamic_df)

# Extract LDA scores
lda_dynamic_scores <- as.data.frame(predict(lda_dynamic_fit)$x)

# Use original model labels, not predicted classes
lda_dynamic_scores$model <- dynamic_df$model

# Compute Niche centroid in LDA space
niche_centroid_dynamic_lda <- lda_dynamic_scores %>%
  filter(model == "Niche") %>%
  summarise(across(starts_with("LD"), mean))

# Centre all LDA scores relative to the Niche centroid
lda_dynamic_scores <- lda_dynamic_scores %>%
  mutate(
    LD1 = LD1 - niche_centroid_dynamic_lda$LD1,
    LD2 = LD2 - niche_centroid_dynamic_lda$LD2
  )

# Plot dynamic LDA
lda_dynamic_plot <- ggplot(
  lda_dynamic_scores,
  aes(x = LD1, y = LD2, colour = model, fill = model)) +
  geom_hline(
    yintercept = 0,
    colour = "#A5ACAF") +
  geom_vline(
    xintercept = 0,
    colour = "#A5ACAF") +
  stat_ellipse(
    level = 0.95,
    linetype = 2,
    show.legend = FALSE) +
  geom_point(
    alpha = 0.7,
    size = 2.5,
    shape = 21) +
  scale_fill_manual(values = model_colours) +
  scale_colour_manual(values = model_colours) +
  labs(
    x = "LD1, centred on Niche",
    y = "LD2, centred on Niche",
    fill = "Model",
    colour = "Model") +
  figure_theme +
  guides(colour = "none") +
  theme(
    legend.position = "right")

lda_dynamic_plot

ggsave("../figures/lda_dynamic_metrics.png",
  lda_dynamic_plot,
  width = 8, height = 6, dpi = 600)


#------5. LDA on topology at equilibrium ---------------------------------------
# Fit LDA using equilibrium topology metrics
lda_topology_fit <- MASS::lda(model ~ ., data = topology_df)

# Extract LDA scores
lda_topology_scores <- as.data.frame(predict(lda_topology_fit)$x)

# Use original model labels, not predicted classes
lda_topology_scores$model <- topology_df$model

# Compute Niche centroid in LDA space
niche_centroid_topology_lda <- lda_topology_scores %>%
  filter(model == "Niche") %>%
  summarise(across(starts_with("LD"), mean))

# Centre all LDA scores relative to the Niche centroid
lda_topology_scores <- lda_topology_scores %>%
  mutate(
    LD1 = LD1 - niche_centroid_topology_lda$LD1,
    LD2 = LD2 - niche_centroid_topology_lda$LD2)

# Plot topology LDA
lda_topology_plot <- ggplot(
  lda_topology_scores,
  aes(x = LD1, y = LD2, colour = model, fill = model)) +
  geom_hline(
    yintercept = 0,
    colour = "#A5ACAF") +
  geom_vline(
    xintercept = 0,
    colour = "#A5ACAF") +
  stat_ellipse(
    level = 0.95,
    linetype = 2,
    show.legend = FALSE) +
  geom_point(
    alpha = 0.7,
    size = 2.5,
    shape = 21 ) +
  scale_fill_manual(values = model_colours) +
  scale_colour_manual(values = model_colours) +
  labs(
    x = "LD1, centred on Niche",
    y = "LD2, centred on Niche",
    fill = "Model",
    colour = "Model"
  ) +
  figure_theme +
  guides(colour = "none") +
  theme(
    legend.position = "right")

lda_topology_plot

ggsave("../figures/lda_equilibrium_topology.png",
  lda_topology_plot,
  width = 8, height = 6, dpi = 600)



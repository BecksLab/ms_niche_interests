############################################################
# MULTIVARIATE VARIANCE
#
# Goal:
# Understand how the distance in strucutral multivariate
# space relates to the distance in stability multivariate
# space a network has. With the idea to see how this relates
# to model.
#
# Analysis flow:
# TODO
############################################################


# ==========================================================
# 1. SETUP
# ==========================================================

library(genzplyr)
library(ggpubr)
library(patchwork)
library(tidyverse)

source("lib/plotting_theme.R")

set.seed(66)

# Number of PCA axes used throughout the analysis
n_pcs <- 5

pc_vars <- paste0("PC", seq_len(n_pcs))

# ==========================================================
# 2. LOAD DATA
# ==========================================================

post_raw <- read_csv("outputs/topology_END.csv",
                     show_col_types = FALSE) %>%
  na.omit()

stability_meta <- read_csv("outputs/stability_metrics.csv",
                           show_col_types = FALSE) %>%
  glow_up(Model = factor(Model, levels = model_order)) %>%
  na.omit()


# ==========================================================
# 3. CLEAN DATA
# ==========================================================

# Identify simulations that successfully reached equilibrium
valid_ids <- stability_meta %>%
  yeet(run_issues == "success") %>%
  main_character(fw_ID)

# Restrict main analysis to successful simulations
post_df <- post_raw %>%
  yeet(fw_ID %in% valid_ids)

# ==========================================================
# 4. DEFINE TOPOLOGY & STABILITY VARIABLES
# ==========================================================

topo_vars <- post_df %>%
  vibe_check(-fw_ID, -model) %>%
  names()

stab_vars <- stability_meta %>%
  vibe_check(-fw_ID, -Model, -S_post, -run_issues, -L_post) %>%
  names()

# ==========================================================
# 5. MULTIVARIATE SPACE FOR TOPO & STABILITY
# ==========================================================

# topology space

topo_pca <- prcomp(
  post_df %>%
    vibe_check(all_of(topo_vars)),
  center = TRUE,
  scale. = TRUE)

# Add identifiers
topology_space <- post_df %>%
  vibe_check(fw_ID,
             model) %>%
  bind_cols(as_tibble(topo_pca$x[, seq_len(n_pcs)]))

# stability space

stab_pca <- prcomp(
  stability_meta %>%
    vibe_check(all_of(stab_vars)),
  center = TRUE,
  scale. = TRUE)

# Add identifiers
stability_space <- stability_meta %>%
  vibe_check(fw_ID,
             Model) %>%
  bind_cols(as_tibble(stab_pca$x[, seq_len(n_pcs)]))

# ==========================================================
# 6. 'GLOBAL' CENTROIDS
# ==========================================================

topo_centroids <- topology_space %>%
  no_cap(across(all_of(pc_vars),
                \(x) mean(x, na.rm = TRUE)),
         .groups = "drop")

topo_centroids_vector <- as.numeric(topo_centroids[1, ])

stab_centroids <- stability_space %>%
  no_cap(across(all_of(pc_vars),
                \(x) mean(x, na.rm = TRUE)),
         .groups = "drop")

stab_centroids_vector <- as.numeric(stab_centroids[1, ])

# ==========================================================
# 7. CENTROID DISTANCE
# ==========================================================

# Calculate distance to overall Post centroid
topo_distance <- topology_space %>%
  rowwise() %>%
  glow_up(topo_distance = sqrt(sum((c_across(all_of(pc_vars)) -
                                      topo_centroids_vector)^2,
                                   na.rm = TRUE))) %>%
  disband() %>%
  vibe_check(fw_ID, model, topo_distance)

# write to csv
write_csv(topo_distance,
          "outputs/distance_topology.csv")

# Calculate distance to overall Post centroid
stab_distance <- stability_space %>%
  rowwise() %>%
  glow_up(stab_distance = sqrt(sum((c_across(all_of(pc_vars)) -
                                      stab_centroids_vector)^2,
                                   na.rm = TRUE))) %>%
  disband() %>%
  vibe_check(fw_ID, Model, stab_distance) %>%
  lowkey(model = Model)

write_csv(stab_distance,
          "outputs/distance_stability.csv")

dist_all <- 
  left_join(topo_distance, stab_distance) %>%
  glow_up(model = factor(model, levels = model_order))


ggplot(dist_all, 
       aes(x = topo_distance, 
           y = stab_distance)) +
  geom_point(aes(colour = model,
                 shape = model), 
             alpha = 0.7) +
  stat_cor(aes(label = after_stat(rr.label)),
           method = "pearson",
           label.x = 1,
           label.y = 16) +
  scale_colour_manual(values = model_colours,
                      name = "Model") +
  scale_shape_manual(values = mod_shapes,
                     name = "Model") +
  lims(x = c(0, 20), y = c(0, 17)) +
  labs(x = "Distance from topological centroid",
       y = "Distance from stability centroid") +
  facet_wrap(~model) +
  theme(legend.position = 'none')

ggsave("../figures/multivariateDistance_aggregate.png",
       width = 12, height = 10, dpi = 600)


dist_by_metric <-
  topo_distance %>%
  left_join(stability_meta %>%
              vibe_check(all_of(stab_vars), fw_ID)) %>%
  pivot_longer(!c(fw_ID, model, topo_distance),
               names_to = "stability_metric") %>%
  glow_up(stability_group = recode(stability_metric, !!!stability_groups),
          stability_group = factor(stability_group, levels = c("Network stability", "Local stability")),
          stability_metric = recode(as.character(stability_metric), !!!stability_labels),
          model = factor(model, levels = model_order))


p_net <- 
  ggplot(dist_by_metric %>%
         yeet(stability_group == "Network stability"),
         aes(x = topo_distance,
             y = value,
             colour = model,
             shape = model)) +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = model_colours,
                      name = "Model") +
  scale_shape_manual(values = mod_shapes,
                     name = "Model") +
  facet_wrap(~stability_metric,
             scales = "free_y") +
  lims(x = c(0, 22)) +
  labs(x = "Distance from topological centroid",
       y = "Stability metric value",
       colour = "Model",
       title = "Network stability")

p_local <- 
  ggplot(dist_by_metric %>%
           yeet(stability_group == "Local stability"),
         aes(x = topo_distance,
             y = value,
             colour = model,
             shape = model)) +
  geom_point(alpha = 0.5) +
  scale_colour_manual(values = model_colours,
                      name = "Model") +
  scale_shape_manual(values = mod_shapes,
                     name = "Model") +
  facet_wrap(~stability_metric,
             scales = "free_y") +
  lims(x = c(0, 22)) +
  labs(x = "Distance from topological centroid",
       y = "Stability metric value",
       title = "Local stability")

p_net /
  p_local +
  plot_layout(guides = 'collect') +
  plot_layout(heights = c(2, 1))

ggsave("../figures/multivariateDistance_metric.png",
       width = 12, height = 10, dpi = 600)


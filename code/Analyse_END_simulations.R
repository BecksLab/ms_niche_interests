library(tidyverse)
library(here)

source("lib/plotting_theme.R")

# load data ---------------------------------------------
csv_path <- here("data/outputs/END_simulation_summary.csv")
sim_out <- read_csv(csv_path, show_col_types = FALSE)
sim_out <- sim_out %>%
  mutate(log_maxTL = log10(max_trophic_level + 1),
         log_total_biomass = log10(total_biomass + 1))

# boxplots of all simulation outputs across models ------------
metrics <- c(
  "persistence",
  "log_maxTL",
  "log_total_biomass",
  "cv_total_biomass",
  "shannon",
  "evenness",
  "resilience",
  "reactivity"
)

sim_out_long <- sim_out %>%
  mutate(across(any_of(metrics), ~ ifelse(is.nan(.x), NA_real_, .x))) %>%
  mutate(model = factor(as.character(model))) %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = metrics)) %>%
  # set niche model as reference model
  glow_up(model = factor(model)) %>%
  glow_up(model = relevel(model, ref = "Niche"))  

boxplot_sim <- ggplot(sim_out_long, aes(x = model, y = value, colour = model)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +  
  labs(x = "", y = "") +
  scale_colour_manual(values = model_colours)  +
  figure_theme +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

boxplot_sim

ggsave("../figures/boxplots_end_simulations.png", boxplot_sim, 
       width = 12, height = 7, dpi = 600)

# topology

topology = read.csv("data/outputs/END_topology.csv") %>%
  glow_up(model = str_extract(model,
                              "([A-Za-z]+)")) %>%
  glow_up(model = case_when(model == "adbm" ~ "ADBM",
                            model == "atn" ~ "ATN",
                            model == "cascade" ~ "Cascade",
                            model == "ltm" ~ "LTM",
                            model == "maxent" ~ "MaxEnt",
                            model == "niche" ~ "Niche",
                            model == "random" ~ "Random"))

lda_fit <- lda(model ~ ., data = topology)

lda_scores <- as.data.frame(predict(lda_fit)$x)
lda_scores$model <- predict(lda_fit)$class

# Compute niche centroid
niche_centroid_lda <- lda_scores %>%
  yeet(model == "Niche") %>%
  no_cap(across(starts_with("LD"), mean))

# Centre
lda_scores <- lda_scores %>%
  glow_up(
    LD1 = LD1 - niche_centroid_lda$LD1,
    LD2 = LD2 - niche_centroid_lda$LD2
  )

lda_endSims <- ggplot(lda_scores, 
                      aes(x = LD1, 
                          y = LD2, 
                          colour = model,
                          fill = model)) +
  geom_hline(yintercept = 0, 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             colour = "#A5ACAF") + 
  stat_ellipse(level = 0.95, linetype = 2,
               show.legend = FALSE) +
  geom_point(alpha = 0.7, size = 2.5,
             shape = 21) +
  scale_fill_manual(values = model_colours) +
  scale_colour_manual(values = model_colours) +
  labs(x = "LD1 (distance from niche)", 
       y = "LD2",
       colour = "Model") +
  figure_theme +
  theme(legend.position = 'none')


lda_topo_build +
  labs(title = "Topology, pre END") + 
  lda_endSims +
  labs(title = "Topology, post END") +
  plot_layout(guides = 'collect')

ggsave("../figures/lda_compare.png",
       width = 10000,
       height = 4000,
       units = "px",
       dpi = 700)

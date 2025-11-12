
rm(list=ls())

library(tidyverse)
library(here)

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
  mutate(metric = factor(metric, levels = metrics))  

boxplot_sim <- ggplot(sim_out_long, aes(x = model, y = value, color = model)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +  
  labs(x = "", y = "") +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )
boxplot_sim

ggsave("../figures/boxplots_end_simulations.png", boxplot_sim, 
       width = 14, height = 9, dpi = 300)





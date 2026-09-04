
library(tidyverse)
library(here)
library(MASS)
library(patchwork)
library(broom)
library(emmeans)
library(multcomp)
library(effectsize)
library(genzplyr)

setwd(here())

source("lib/plotting_theme.R")

#------1. Load data ------------------------------------------------------------
model_order <- c("ADBM", "ATN", "LTM", "Niche", "Cascade", "MaxEnt", "Random")
model_colours <- model_colours[model_order]

stability_meta <- read_csv("outputs/stability_metrics.csv", show_col_types = FALSE) %>%
  mutate(Model = factor(Model, levels = model_order))

stability_metrics <- c("persistence", "biomass_shannon", "gini_fluxes", "skewness_IS",
                 "resilience","reactivity")
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
#------2. Difference checks across models --------------------------------------
# MANOVA check
stability_df <- stability_meta %>%
  vibe_check(Model, all_of(stability_metrics)) %>%
  na.omit()

stability_dep <- stability_df %>%
  vibe_check(all_of(stability_metrics)) %>%
  as.matrix()

manova <- manova(stability_dep ~ Model, data = stability_df)
summary(manova, test = "Pillai")

#------3. Boxplot of all simulation outputs across models ---------------------
stability_long <- stability_df %>%
  pivot_longer(
    cols = all_of(stability_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    metric = factor(metric, levels = stability_metrics),
    
    metric_label = recode(as.character(metric), !!!stability_labels),
    metric_label = factor(metric_label, levels = stability_labels[stability_metrics]),
    
    stability_group = recode(as.character(metric), !!!stability_groups),
    stability_group = factor(
      stability_group,
      levels = c("Network stability", "Local stability")
    )
  )

p_network <- stability_long %>%
  filter(stability_group == "Network stability") %>%
  ggplot(aes(x = Model, y = value, colour = Model)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~ metric_label, scales = "free_y", nrow = 2, ncol = 2) +
  labs(x = "", y = "Network stability") +
  scale_colour_manual(values = model_colours) +
  figure_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

p_local <- stability_long %>%
  filter(stability_group == "Local stability") %>%
  ggplot(aes(x = Model, y = value, colour = Model)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~ metric_label, scales = "free_y", nrow = 1, ncol = 2) +
  labs(x = "", y = "Local stability") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = scales::breaks_extended(n = 4),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_colour_manual(values = model_colours) +
  figure_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

stability_boxplot <- p_network / p_local +
  plot_layout(heights = c(2, 1))

stability_boxplot

ggsave("../figures/stability_boxplot.png", stability_boxplot, 
       width = 12, height = 10, dpi = 600)

#------4. Difference in group means and variance contribution ------------------
effect_list <- list()
effect_size_list <- list()
emm_list <- list()

for (m in stability_metrics) {
  
  df_m <- stability_meta %>%
    filter(!is.na(.data[[m]]), !is.na(S_post))
  
  mod <- lm(as.formula(paste(m, "~ Model + S_post")), data = df_m)
  anova_m <- car::Anova(mod, type = 2)
  emm <- emmeans(mod, specs = "Model") # at the same S_post value (mean)
  
  effect_list[[m]] <- broom::tidy(anova_m) %>% mutate(metric = m)
  
  effect_size_list[[m]] <- effectsize::eta_squared(anova_m, partial = TRUE) %>%
    as.data.frame() %>% mutate(metric = m)
  
  emm_list[[m]] <- multcomp::cld(emm, Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>% mutate(metric = m)
}

effect_df <- bind_rows(effect_list)
effect_size_df <- bind_rows(effect_size_list)
emm_df <- bind_rows(emm_list)

effect_size_summary <- effect_size_df %>%
  mutate(metric_label = recode(metric, !!!stability_labels),
         metric_label = factor(metric_label, levels = stability_labels[stability_metrics]),
         Parameter = recode(Parameter,
                            Model = "Model",
                            S_post = "Equilibrium richness")) %>%
  dplyr::select(metric, metric_label, Parameter, Eta2_partial, CI_low, CI_high) %>%
  arrange(metric, desc(Eta2_partial))

effect_size_wide <- effect_size_summary %>%
  dplyr::select(metric_label, Parameter, Eta2_partial) %>%
  pivot_wider(names_from = Parameter, values_from = Eta2_partial)

emm_summary <- emm_df %>%
  mutate(metric_label = recode(metric, !!!stability_labels),
         metric_label = factor(metric_label, levels = stability_labels[stability_metrics]),
         .group = str_trim(.group)) %>%
  dplyr::select(metric, metric_label, Model, emmean, SE, df, lower.CL, upper.CL, .group)

effect_df
effect_size_summary
effect_size_wide
emm_summary

#------5. Plot effect contribution and adjusted model means ------------
effect_size_plot <- effect_size_summary %>%
  filter(Parameter %in% c("Model", "Equilibrium richness")) %>%
  ggplot(aes(x = metric_label, y = Eta2_partial, fill = Parameter)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  coord_flip() +
  labs(x = "", y = "Partial eta-squared", fill = "") +
  scale_fill_manual(values = c(col_div[1], col_div[3])) +
  figure_theme +
  theme(legend.position = "top")

effect_size_plot

ggsave("../figures/stability_effect_size_model_vs_richness.png",
       effect_size_plot, width = 8, height = 5, dpi = 600)


model_colours <- model_colours[model_order]

emm_plot_df <- emm_summary %>%
  mutate(
    metric_label = factor(metric_label, levels = stability_labels[stability_metrics]),
    stability_group = recode(metric, !!!stability_groups),
    stability_group = factor(stability_group, levels = c("Network stability", "Local stability"))
  ) %>%
  group_by(metric_label) %>%
  mutate(
    Model = factor(Model, levels = model_order),
    letter_y = upper.CL + 0.08 * (max(upper.CL, na.rm = TRUE) - min(lower.CL, na.rm = TRUE))
  ) %>%
  ungroup()

emm_network_plot <- emm_plot_df %>%
  filter(stability_group == "Network stability") %>%
  ggplot(aes(x = Model, y = emmean, colour = Model)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                  linewidth = 0.8, size = 1.1) +
  geom_text(aes(y = letter_y, label = .group),
            colour = "black", size = 4.5, show.legend = FALSE) +
  facet_wrap(~ metric_label, scales = "free_y", 
             nrow = 2, ncol = 2) +
  scale_colour_manual(values = model_colours, limits = model_order) +
  labs(x = "", y = "Adjusted mean", title = "Network stability") +
  figure_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

emm_local_plot <- emm_plot_df %>%
  filter(stability_group == "Local stability") %>%
  ggplot(aes(x = Model, y = emmean, colour = Model)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                  linewidth = 0.8, size = 1.1) +
  geom_text(aes(y = letter_y, label = .group),
            colour = "black", size = 4.5, show.legend = FALSE) +
  facet_wrap(~ metric_label, scales = "free_y", 
             nrow = 1, ncol = 2) +
  scale_colour_manual(values = model_colours, 
                      limits = model_order) +
  labs(x = "", y = "Adjusted mean", title = "Local stability") +
  figure_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

adjusted_mean_plot <- emm_network_plot / emm_local_plot +
  plot_layout(heights = c(2, 1))

adjusted_mean_plot

ggsave("../figures/stability_adjusted_model_means.png",
       adjusted_mean_plot, width = 12, height = 10, dpi = 600)



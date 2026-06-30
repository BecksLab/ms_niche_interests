library(emmeans)
library(tidyverse)
library(here)
library(MASS)

setwd(here())

source("lib/plotting_theme.R")

#------1. Load data ------------------------------------------------------------
dynamic_metric <- read_csv(
  "outputs/dynamic_metrics.csv",
  show_col_types = FALSE)

#------2. Difference checks across models --------------------------------------

# MANOVA for dynamic outputs
dyn_metrics <- dynamic_metric %>%
  vibe_check(-fw_ID, -model, -S_post, -L_post) %>%
  names()

dynamic_df <- dynamic_metric %>%
  vibe_check(model, all_of(dyn_metrics)) %>%
  na.omit()

dynamic_dep <- dynamic_df %>%
  vibe_check(all_of(dyn_metrics)) %>%
  as.matrix()

dynamic_manova <- manova(dynamic_dep ~ model, data = dynamic_df)
summary(dynamic_manova, test = "Pillai")

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
  pivot_longer(
    cols = all_of(dyn_metrics), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    # Keep facet order following dyn_metrics
    metric = factor(metric, levels = dyn_metrics)#,
    
    # Rename metrics using the same vector
    #metric_label = recode(as.character(metric), !!!dyn_metric_labels),
    
    # Keep renamed facet labels in the same order
    #metric_label = factor(metric_label, levels = dyn_metric_labels[dyn_metrics])
  ) 

boxplot_sim <- ggplot(dynamic_out_long, aes(x = model, y = value, colour = model)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y") +  
  labs(x = "", y = "") +
  scale_colour_manual(values = model_colours)  +
  figure_theme +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

boxplot_sim

ggsave("../figures/dynamics_boxplot.png", boxplot_sim, 
       width = 12, height = 7, dpi = 600)

#------4. Difference in group means --------------------------------------------

# Fit model and calculate EMMs for each metric
emm_list <- list()

for (m in dyn_metrics) {
  
  formula <- as.formula(paste(m, "~ model"))
  mod <- lm(formula, data = dynamic_df)
  
  # Calculate Estimated Marginal Means and compact letter display
  emm <- emmeans(mod, specs = "model")
  cld_res <- multcomp::cld(emm, Letters = letters)
  
  emm_list[[m]] <- as.data.frame(cld_res) %>% mutate(metric = m)
}

# Convert the EMMs list to a df
emm_df <- bind_rows(emm_list)


# plot
plot_local <- 
  ggplot(emm_df %>%
           yeet(metric %in% c("resilience", "reactivity")), 
       aes(x = model, 
           y = emmean, 
           color = model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.2, 
                linewidth = 0.8) +
  geom_text(aes(label = .group, 
                y = upper.CL),
            vjust = -0.6, size = 3.5, 
            colour = shark_black) +
  facet_wrap(~ metric,
             scales = "free_y") +
  labs(x = NULL,
       y = "Predicted value",
       title = 'Local Stability')  +
  scale_colour_manual(values = model_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# plot
plot_net <- ggplot(emm_df %>%
         yeet(metric %in% c("gini_fluxes_formula", "persistence", "skewness_IS")), 
       aes(x = model, 
           y = emmean, 
           color = model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.2, 
                linewidth = 0.8) +
  geom_text(aes(label = .group, 
                y = upper.CL),
            vjust = -0.6, size = 3.5, 
            colour = shark_black) +
  facet_wrap(~ metric,
             scales = "free_y",
             ncol = 2) +
  labs(x = NULL,
       y = "Predicted value",
       title = "Network Stability")  +
  scale_colour_manual(values = model_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_local /
  plot_net +
  plot_layout(guides='collect',
              heights = c(1, 2))

ggsave("../figures/dynamics_emmeans.png", 
       width = 9, 
       height = 10, 
       dpi = 600)

#------5. LDA on dynamic metrics -----------------------------------------------
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

ggsave("../figures/dynamics_lda.png",
       lda_dynamic_plot,
       width = 8, height = 6, dpi = 600)



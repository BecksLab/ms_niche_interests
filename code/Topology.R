# libraries
library(ggpubr)
library(here)
library(effectsize)
library(genzplyr)
library(ggrepel)
library(ggtext)
library(MASS)
library(patchwork)
library(rstatix)
library(tidyverse)
library(vegan)
library(MVN)        # For Henze-Zirkler test
library(biotools)   # For Box's M test
library(candisc)    # For Canonical Discriminant Analysis
library(emmeans)    # For pairwise comparisons
library(lmerTest)   # For mixed models

# set path to code sub dir
setwd(here())

source("lib/plotting_theme.R")

# import data
df <- read_csv("data/outputs/topology.csv") %>%
  vibe_check(-c(richness)) %>%
  na.omit() %>%
  # set niche model as reference model
  glow_up(model = factor(model)) %>%
  glow_up(model = relevel(model, ref = "Niche"))

# Dependent variable matrix for multivariate tests
# Assuming columns 2 onwards are your topological metrics
dep_vars <- as.matrix(df[2:ncol(df)])

# =========================
# 1. MANOVA + Assumption Checks
# =========================
fit <- manova(dep_vars ~ model, data = df)

# --- Multivariate test
summary(fit, test = "Pillai")

# --- Assumption Checks
# Henze-Zirkler Multivariate Normality
result_hz <- hz(data = dep_vars)
print(result_hz)

# Box's M test for Homogeneity of Covariance
boxM(dep_vars, df$model) 

# =========================
# 2. Canonical Discriminant Analysis (CDA)
# =========================
# This replaces the simple LDA with a more robust canonical variate approach
cda <- candisc(fit)
summary(cda)

# shift everything relative to niche

# Get canonical scores
scores <- as.data.frame(cda$scores)
scores$model <- df$model

# Compute niche centroid
niche_centroid <- scores %>%
  yeet(model == "Niche") %>%
  no_cap(across(starts_with("Can"), mean))

# Centre all scores on niche
scores_centered <- scores %>%
  glow_up(
    Can1 = Can1 - niche_centroid$Can1,
    Can2 = Can2 - niche_centroid$Can2
  )

# Extract canonical loadings for plotting
loadings_df <- as.data.frame(cda$structure[, 1:2]) %>%
  rownames_to_column("Metric") %>%
  rename(CV1 = Can1, CV2 = Can2) %>%
  glow_up(Level = case_when(
    Metric %in% c("complexity", "connectance", "trophic_level") ~ "Macro",
    Metric %in% c("generality", "vulnerability") ~ "Micro",
    TRUE ~ "Meso"
  ))

# CDA Loadings Plot
ggplot(loadings_df, aes(x = CV1, y = CV2)) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = CV1, yend = CV2, color = Level),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_text_repel(aes(label = Metric)) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  scale_colour_manual(values = c("#006D75", "#2F2F2F", "#EA7200")) +
  figure_theme

ggsave("../figures/lda_corr.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

# =========================
# 3. LDA Visualization
# =========================
# Using LDA for the group separation plot
lda_fit <- lda(model ~ ., data = df)
lda_scores <- predict(lda_fit)$x

plot_lda <- data.frame(
  model = df$model,
  lda = lda_scores
)

ggplot(scores_centered, aes(x = Can1, y = Can2, 
                            fill = model, colour = model)) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#A5ACAF") +
  geom_point(alpha = 0.7, size = 2.5,
             colour = "white",
             shape = 21) +
  stat_ellipse(level = 0.95, linetype = 4,
               show.legend = FALSE) +
  scale_colour_manual(values = model_colours) +
  scale_fill_manual(values = model_colours) +
  labs(x = "CV1 (distance from niche)",
       y = "CV2 (distance from niche)",
       fill = "Model") +
  figure_theme

ggsave("../figures/cv.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

# niche centred LDA

lda_fit <- lda(model ~ ., data = df)
lda_scores <- as.data.frame(predict(lda_fit)$x)
lda_scores$model <- df$model

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

lda_topo_build <-
  ggplot(lda_scores, 
         aes(x = LD1, 
             y = LD2, 
             colour = model)) +
  geom_hline(yintercept = 0, 
             colour = "#A5ACAF") +
  geom_vline(xintercept = 0, 
             colour = "#A5ACAF") + 
  stat_ellipse(level = 0.95, linetype = 2) +
  geom_point(alpha = 0.6, size = 2) +
  scale_colour_manual(values = model_colours) +
  labs(x = "LD1 (distance from niche)", 
       y = "LD2",
       colour = "Model") +
  figure_theme

ggsave("../figures/lda.png",
       lda_topo_build,
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

# =========================
# 4. Pairwise Comparisons (EMMeans)
# =========================
# This performs univariate comparisons between models for specific metrics
metrics_to_test <- df %>% select(!model) %>% names()

emm_list <- list()

for (m in metrics_to_test) {
  # We use a linear model here. If you have a 'time' or 'replicate' random effect, 
  # use lmer(as.formula(paste(m, "~ model + (1|replicate)")), data = df)
  formula <- as.formula(paste(m, "~ model"))
  mod <- lm(formula, data = df)
  
  # Calculate Estimated Marginal Means and compact letter display
  emm <- emmeans(mod, specs = "model")
  cld_res <- multcomp::cld(emm, Letters = letters)
  
  emm_list[[m]] <- as.data.frame(cld_res) %>% mutate(metric = m)
}

emm_df <- bind_rows(emm_list)

# --- 1. Define Metric Categories ---
# This matches the levels defined in your 06_structural_differences logic
emm_df <- emm_df %>%
  mutate(stat_label = case_when(
    metric == "connectance" ~ "Connectance",
    metric == "trophic_level" ~ "Max Trophic Level",
    metric == "generality" ~ "Generality",
    metric == "vulnerability" ~ "Vulnerability",
    metric == "S1" ~ "No. of linear chains",
    metric == "S2" ~ "No. of omnivory motifs",
    metric == "S5" ~ "No. of apparent competition motifs",
    metric == "S4" ~ "No. of direct competition motifs",
    TRUE ~ as.character(metric)
  )) %>%
  mutate(level = case_when(
    metric %in% c("complexity", "connectance", "trophic_level") ~ "Macro",
    metric %in% c("generality", "vulnerability") ~ "Micro",
    TRUE ~ "Meso"
  ))

# --- 2. Create Categorized Plot List ---
levs <- c("Macro", "Meso", "Micro")
plot_list_emm <- vector("list", length = 3)

for (i in seq_along(levs)) {
  
  plot_list_emm[[i]] <- ggplot(
    emm_df %>% filter(level == levs[i]),
    aes(x = model, 
        y = emmean, 
        colour = model)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower.CL, 
                      ymax = upper.CL), 
                  width = 0.2, linewidth = 0.6) +
    # Add Compact Letter Display (CLD) for significance
    geom_text(aes(label = .group, 
                  y = upper.CL),
              vjust = -0.6, size = 3.5, colour = "#00181F") +
    facet_wrap(vars(stat_label), scales = "free_y", ncol = 2) +
    scale_colour_manual(values = model_colours) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = NULL, 
         y = "Estimated Mean", 
         title = levs[i]) +
    figure_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

# --- 3. Patchwork Assembly ---
# Adjusting heights to accommodate the number of facets in each level
(plot_list_emm[[1]] / plot_list_emm[[2]] / plot_list_emm[[3]]) +
  plot_layout(heights = c(2, 2, 1), guides = "collect")

# Save
ggsave("../figures/emm_summary.png", width = 9, height = 12, dpi = 400)

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
df <- read_csv("data/outputs/topology_initial_05-06-2026.csv") %>%
  vibe_check(-c(richness, S_initial, connectance_initial, prop_basal_initial,
                mx_tl_initial, mean_FCL_initial, closed_loops_initial, isolated_initial,
                illogical_initial, fw_ID)) %>%
  yeet(model != "LTM_Laura") %>%
  yeet(max_trophic_level < 16) %>%
  na.omit() %>%
  # set niche model as reference model
  glow_up(model = factor(model)) %>%
  glow_up(model = relevel(model, ref = "Niche"))

df %>%
  squad_up(model) %>%
  tally()

dynamic_topo <- read_csv("data/outputs/END_topology_03_05_2026.csv",
                          show_col_types = FALSE) %>%
  vibe_check(-richness) %>%
  glow_up(model = factor(model),
          model = relevel(model, ref = "Niche")) %>%
  na.omit()

# Dependent variable matrix for multivariate tests
# Assuming columns 2 onwards are your topological metrics
dep_vars <- as.matrix(df[2:ncol(df)])

# general boxplot
ggplot(df %>%
         pivot_longer(-model)) +
  geom_boxplot(aes(x = model,
                   y = value,
                   colour = model),
               outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(vars(name), scales = "free_y", ncol = 2) +
  scale_colour_manual(values = model_colours) +
  figure_theme


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

# Center all scores on niche
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
    Metric %in% c("complexity", "connectance", "max_trophic_level", "ChLen") ~ "Macro",
    Metric %in% c("generality") ~ "Role",
    Metric %in% c("vulnerability", "top") ~ "Heterogeneity",
    Metric %in% c("distance") ~ "Path",
    TRUE ~ "Scaling"
  ))

# CDA Loadings Plot
ggplot(loadings_df, aes(x = CV1, y = CV2)) +
  geom_hline(yintercept = 0, 
             colour = shark_silver) +
  geom_vline(xintercept = 0, 
             colour = shark_silver) +
  geom_segment(aes(x = 0, y = 0, xend = CV1, yend = CV2, color = Level),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_text_repel(aes(label = Metric)) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  scale_colour_manual(values = c("#006D75", "#2F2F2F", "#EA7200", "#B2B4B2", "#FFB81C")) +
  figure_theme

ggsave("../figures/cda_corr.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

# =========================
# 3. CDA Visualization
# =========================
# Using the canonical variates (from candisc) for group separation plot
topology_plot <- ggplot(scores_centered, aes(x = Can1, y = Can2, colour = model)) +
  geom_hline(yintercept = 0, colour = shark_silver) +
  geom_vline(xintercept = 0, colour = shark_silver) +
  geom_point(alpha = 0.7, size = 2.5) +
  stat_ellipse(level = 0.95, linetype = 4, show.legend = FALSE) +
  scale_colour_manual(values = model_colours) +
  labs(x = "CV1, centred on Niche",
       y = "CV2, centred on Niche",
       colour = "Model") +
  figure_theme

# =========================
# 4. CDA on topology at equilibrium
# =========================
# Fit CDA using equilibrium topology metrics
dep_vars <- as.matrix(dynamic_topo[2:ncol(dynamic_topo)])
dynam_fit <- manova(dep_vars ~ model, data = dynamic_topo)

dynam_cda <- candisc(dynam_fit)
summary(dynam_cda)

# shift everything relative to niche

# Get canonical scores
dynam_scores <- as.data.frame(dynam_cda$scores)
dynam_scores$model <- dynamic_topo$model

# Compute niche centroid
dynam_niche_centroid <- dynam_scores %>%
  yeet(model == "Niche") %>%
  no_cap(across(starts_with("Can"), mean))

# Center all scores on niche
dynam_scores_centered <- dynam_scores %>%
  glow_up(
    Can1 = Can1 - dynam_niche_centroid$Can1,
    Can2 = Can2 - dynam_niche_centroid$Can2
  )

# Plot topology LDA
dynam_topology_plot <- ggplot(
  dynam_scores_centered, aes(x = Can1, y = Can2, colour = model)) +
  geom_hline(yintercept = 0, colour = "#A5ACAF") +
  geom_vline(xintercept = 0, colour = "#A5ACAF") +
  geom_point(alpha = 0.7, size = 2.5) +
  stat_ellipse(level = 0.95, linetype = 4, show.legend = FALSE) +
  scale_colour_manual(values = model_colours) +
  labs(x = "CV1, centred on Niche",
       y = "CV2, centred on Niche",
       colour = "Model") +
  figure_theme

# =========================
# 5. Pairwise Comparisons (EMMeans)
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
    metric %in% c("complexity", "connectance", "max_trophic_level", "S2", "S1", "ChLen") ~ "Macro",
    metric %in% c("generality") ~ "Role",
    metric %in% c("vulnerability", "top") ~ "Heterogeneity",
    metric %in% c("distance") ~ "Path",
    metric %in% c("centrality") ~ "Hubs",
    TRUE ~ "Scaling"
  ))

# --- 2. Create Categorized Plot List ---
levs <- c("Macro", "Role", "Heterogeneity", "Path", "Scaling")
plot_list_emm <- vector("list", length = length(levs))

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
(plot_list_emm[[1]] / (plot_list_emm[[2]] + plot_list_emm[[4]]) / 
    plot_list_emm[[3]] / plot_list_emm[[5]])+
  plot_layout(heights = c(3, 1, 1, 1), guides = "collect")

# Save
ggsave("../figures/emm_summary.png", width = 9, height = 12, dpi = 400)

# =========================
# 6. Standard Errors
# =========================

sims <- 
  read_csv("data/outputs/topology.csv") %>%
  na.omit() %>%
  pivot_longer(-c(model),
               names_to = "metric") %>%
  squad_up(model, metric)

# use Niche model as reference point
ref_stats <- sims %>%
  yeet(model == "Niche") %>%
  squad_up(metric) %>%
  no_cap(mu_ref = mean(value),
         sd_ref = sd(value))

comparison <- sims %>%
  squad_up(model, metric) %>%
  no_cap(mu_model = mean(value, na.rm = TRUE)) %>%
  left_join(ref_stats, by = "metric") %>%
  glow_up(z = (mu_model - mu_ref) / sd_ref)

ggplot(comparison %>%
         yeet(metric != "richness")) +
  geom_hline(yintercept = 2, 
             linetype = "dashed", 
             colour = shark_silver) +
  geom_hline(yintercept = -2, 
             linetype = "dashed", 
             colour = shark_silver) +
  geom_point(aes(x = metric,
                 y = z,
                 colour = model),
             alpha = 0.7) +
  scale_colour_manual(values = model_colours) +
  labs(x = "Model",
       y = "Normalised Error (niche as reference)") +
  figure_theme +
  theme(panel.grid.major = element_blank(),
        strip.text = ggtext::element_markdown())

# =========================
# 7. Combine LDAs
# =========================

topology_plot + labs(title = "Topology") + 
  dynam_topology_plot + labs(title = "Topology at equilibrium") + 
  plot_layout(guides='collect')

ggsave("../figures/cda_compare.png",
       width = 10000,
       height = 4000,
       units = "px",
       dpi = 700)


# =========================
# 8. Mahalanobis distance
# =========================

pre_cov <- scores_centered %>%
  filter(model == "Niche") %>%
  select(Can1, Can2) %>%
  cov()

eq_cov <- dynam_scores_centered %>%
  filter(model == "Niche") %>%
  select(Can1, Can2) %>%
  cov()

centroids_maha <-
  centroids %>%
  left_join(niche_ref, by = "type") %>%
  rowwise() %>%
  glow_up(
    maha = sqrt(
      mahalanobis(
        x = matrix(c(Can1, Can2), nrow = 1),
        center = c(niche_Can1, niche_Can2),
        cov = if (type == "Pre") pre_cov else eq_cov
      )
    )
  ) %>%
  disband() %>%
  glow_up(type = factor(type, levels = c("Pre", "Equilibrium")))

centroids_maha %>%
  filter(model != "Niche") %>%
  ggplot() +
  geom_segment(
    aes(
      x = 0,
      xend = maha,
      y = model,
      yend = model,
      colour = model
    )
  ) +
  geom_point(
    aes(
      x = maha,
      y = model,
      colour = model
    ),
    size = 3
  ) +
  facet_wrap(vars(type), ncol = 1) +
  scale_colour_manual(values = model_colours) +
  labs(
    x = "Mahalanobis distance from Niche centroid",
    y = NULL
  ) +
  figure_theme +
  theme(legend.position = "none")

# Which models are significantly further from Niche than Niche itself

pre_mu <- scores_centered %>%
  filter(model == "Niche") %>%
  select(Can1, Can2) %>%
  colMeans()

pre_cov <- scores_centered %>%
  filter(model == "Niche") %>%
  select(Can1, Can2) %>%
  cov()

pre_dist <- scores_centered %>%
  mutate(
    dist_niche =
      sqrt(
        mahalanobis(
          select(., Can1, Can2),
          center = pre_mu,
          cov = pre_cov
        )
      ),
    type = "Pre"
  )

eq_mu <- dynam_scores_centered %>%
  filter(model == "Niche") %>%
  select(Can1, Can2) %>%
  colMeans()

eq_cov <- dynam_scores_centered %>%
  filter(model == "Niche") %>%
  select(Can1, Can2) %>%
  cov()

eq_dist <- dynam_scores_centered %>%
  mutate(
    dist_niche =
      sqrt(
        mahalanobis(
          select(., Can1, Can2),
          center = eq_mu,
          cov = eq_cov
        )
      ),
    type = "Equilibrium"
  )

dist_df <- bind_rows(
  pre_dist,
  eq_dist
) %>%
  mutate(
    type = factor(
      type,
      levels = c("Pre", "Equilibrium")
    )
  )

ggplot(
  dist_df,
  aes(
    x = model,
    y = dist_niche,
    colour = model
  )
) +
  geom_boxplot(
    outlier.alpha = 0.4
  ) +
  facet_wrap(
    vars(type)
  ) +
  scale_colour_manual(
    values = model_colours
  ) +
  labs(
    x = NULL,
    y = "Distance from Niche centroid"
  ) +
  figure_theme

# =========================
# PRE
# =========================

pre_mod <- lm(
  dist_niche ~ model,
  data = filter(dist_df, type == "Pre")
)

pre_emm <- emmeans(
  pre_mod,
  ~ model
)

pre_contrasts <- contrast(
  pre_emm,
  method = "trt.vs.ctrl",
  ref = "Niche",
  adjust = "holm"
)

pre_results <- summary(pre_contrasts) %>%
  as_tibble() %>%
  mutate(type = "Pre")

pre_cld <- multcomp::cld(
  pre_emm,
  Letters = letters
) %>%
  as.data.frame() %>%
  mutate(type = "Pre")

# =========================
# EQUILIBRIUM
# =========================

eq_mod <- lm(
  dist_niche ~ model,
  data = filter(dist_df, type == "Equilibrium")
)

eq_emm <- emmeans(
  eq_mod,
  ~ model
)

eq_contrasts <- contrast(
  eq_emm,
  method = "trt.vs.ctrl",
  ref = "Niche",
  adjust = "holm"
)

eq_results <- summary(eq_contrasts) %>%
  as_tibble() %>%
  mutate(type = "Equilibrium")

eq_cld <- multcomp::cld(
  eq_emm,
  Letters = letters
) %>%
  as.data.frame() %>%
  mutate(type = "Equilibrium")

dist_results <- bind_rows(
  pre_results,
  eq_results
)

cld_df <- bind_rows(
  pre_cld,
  eq_cld
) %>%
  glow_up(type = factor(type, levels = c("Pre", "Equilibrium")))

ggplot(
  cld_df,
  aes(
    model,
    emmean,
    colour = model
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    width = 0.2
  ) +
  geom_text(
    aes(
      label = .group,
      y = upper.CL
    ),
    vjust = -0.6
  ) +
  facet_wrap(vars(type)) +
  scale_colour_manual(values = model_colours) +
  labs(
    y = "Mahalanobis distance from Niche",
    x = NULL
  ) +
  figure_theme +
  theme(legend.position = 'none')

ggsave("../figures/cda_emm_distance.png",
       width = 7000,
       height = 3000,
       units = "px",
       dpi = 700)

pre_results <- confint(pre_contrasts) %>%
  as_tibble() %>%
  mutate(type = "Pre")

eq_results <- confint(eq_contrasts) %>%
  as_tibble() %>%
  mutate(type = "Equilibrium")

bind_rows(
  pre_results,
  eq_results
) %>%
  separate(
    contrast,
    into = c("model", "ref"),
    sep = " - "
  ) %>%
  squad_up(type) %>%
  glow_up(model = forcats::fct_reorder(model, estimate)) %>%
  disband() %>%
  glow_up(type = factor(type, levels = c("Pre", "Equilibrium"))) %>%
  ggplot(
    aes(
      x = estimate,
      y = model,
      colour = model
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = model_colours["Niche"]
  ) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      xmin = lower.CL,
      xmax = upper.CL
    ),
    width = 0.15
  ) +
  facet_wrap(
    vars(type),
    ncol = 1
  ) +
  scale_colour_manual(values = model_colours) +
  labs(
    x = "Additional Mahalanobis distance relative to Niche",
    y = NULL
  ) +
  figure_theme +
  theme(
    legend.position = "none"
  )

# 9. LOADINGS

# Extract canonical loadings for plotting
loadings <- as.data.frame(dynam_cda$structure[, 1:2]) %>%
  glow_up(type = "Equilibrium") %>%
  rownames_to_column("Metric") %>%
  rbind(as.data.frame(cda$structure[, 1:2]) %>%
          glow_up(type = "Pre") %>%
          rownames_to_column("Metric")) %>%
  rename(CV1 = Can1, CV2 = Can2) %>%
  glow_up(Level = case_when(Metric %in% c("complexity", "connectance", "max_trophic_level", "ChLen") ~ "Macro",
                            Metric %in% c("generality") ~ "Role",
                            Metric %in% c("vulnerability", "top") ~ "Heterogeneity",
                            Metric %in% c("distance") ~ "Path",
                            TRUE ~ "Scaling"),
          type = factor(type, levels = c("Pre", "Equilibrium")))

# CDA Loadings Plot
ggplot(loadings, aes(x = CV1, y = CV2)) +
  geom_hline(yintercept = 0, 
             colour = shark_silver) +
  geom_vline(xintercept = 0, 
             colour = shark_silver) +
  geom_segment(aes(x = 0, y = 0, xend = CV1, yend = CV2, color = Level),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_text_repel(aes(label = Metric)) +
  facet_wrap(vars(type),
             ncol = 2) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  scale_colour_manual(values = c("#006D75", "#2F2F2F", "#EA7200", "#B2B4B2", "#FFB81C")) +
  figure_theme


ggsave("../figures/cda_loadings_compare.png",
       width = 10000,
       height = 5000,
       units = "px",
       dpi = 700)

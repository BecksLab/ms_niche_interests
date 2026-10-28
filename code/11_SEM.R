rm(list = ls())

library(tidyverse)
library(here)
library(lavaan)
library(patchwork)
library(diy.sem.plot)

source("lib/plotting_theme.R")

#------1. Load and prepare data ------------------------------------------------
distance_topology <- read.csv(here("outputs", "distance_topology.csv"))

stability_metrics <- read.csv(here("outputs", "stability_metrics.csv")) %>%
  rename(model = Model)

model_order <- c("ADBM", "ATN", "LTM", "Niche", "Cascade", "MaxEnt", "Random")

sem_df <- distance_topology %>%
  left_join(stability_metrics, by = c("fw_ID", "model")) %>%
  mutate(model = factor(model, levels = model_order))

sem_vars <- c("topo_distance", "S_post", "biomass_shannon", "gini_fluxes",
              "skewness_IS", "resilience", "reactivity")

# Standardise across all generative models using common means and SDs
sem_data <- sem_df %>%
  mutate(across(all_of(sem_vars), ~ as.numeric(scale(.x))))


#------2. Multi-group SEM and comparisons of path coefficients ------------------
sem_structure <- '
  topo_distance ~ S_post
  biomass_shannon ~ topo_distance + S_post
  gini_fluxes ~ biomass_shannon + topo_distance + S_post
  skewness_IS ~ biomass_shannon + topo_distance + S_post
  resilience ~ gini_fluxes + skewness_IS + biomass_shannon + topo_distance + S_post
  reactivity ~ gini_fluxes + skewness_IS + biomass_shannon + topo_distance + S_post

  gini_fluxes ~~ skewness_IS
  resilience ~~ reactivity
'

# Same structure, separate coefficients for each generative model
fit_sem <- sem(
  sem_structure, data = sem_data, group = "model", group.label = model_order,
  estimator = "MLR", meanstructure = TRUE, fixed.x = FALSE, auto.cov.y = FALSE
)

summary(fit_sem, standardized = TRUE, rsquare = TRUE)

model_names <- lavInspect(fit_sem, "group.label")
coef_vcov <- vcov(fit_sem)

param_table <- parameterTable(fit_sem) %>%
  filter(op == "~") %>%
  select(lhs, rhs, group, free)

path_df <- parameterEstimates(fit_sem, ci = TRUE) %>%
  filter(op == "~") %>%
  left_join(param_table, by = c("lhs", "rhs", "group")) %>%
  transmute(model = model_names[group], group, response = lhs, predictor = rhs,
            beta = est, se, ci_lower = ci.lower, ci_upper = ci.upper,
            p_value = pvalue, free,
            significance = case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**",
                                     pvalue < 0.05 ~ "*", TRUE ~ ""))

r2_list <- lavInspect(fit_sem, "r2")

r2_df <- map_dfr(seq_along(model_names), function(g) {
  tibble(model = model_names[g], response = names(r2_list[[g]]),
         R2 = unname(r2_list[[g]]))
})

# Omnibus Wald test: does each path differ among generative models?
# Significant: evidence that at least one model has a different coefficient for that path.
# Nonsignificant: insufficient evidence of a difference.
path_names <- path_df %>% distinct(response, predictor)
comparison_list <- list()

for (j in seq_len(nrow(path_names))) {
  df_j <- path_df %>%
    filter(response == path_names$response[j],
           predictor == path_names$predictor[j]) %>%
    arrange(group)
  
  contrast <- cbind(-1, diag(nrow(df_j) - 1))
  differences <- contrast %*% df_j$beta
  contrast_vcov <- contrast %*% coef_vcov[df_j$free, df_j$free] %*% t(contrast)
  wald_stat <- as.numeric(t(differences) %*% solve(contrast_vcov, differences)) #Calculate the Wald statistic
  
  comparison_list[[j]] <- tibble(
    response = path_names$response[j], predictor = path_names$predictor[j],
    chisq = wald_stat, df = nrow(contrast),
    p_value = pchisq(wald_stat, df = nrow(contrast), lower.tail = FALSE)
  )
}

comparison_df <- bind_rows(comparison_list) %>%
  mutate(p_holm = p.adjust(p_value, method = "holm"))


#------3. Decompose stability variance through the SEM pathways -----------------
stability_responses <- c("resilience", "reactivity")

source_blocks <- list(
  Richness = "S_post",
  Topology = "topo_distance",
  BiomassInteractions = c("biomass_shannon", "gini_fluxes", "skewness_IS")
)

component_labels <- c(
  Full = "Full",
  Richness = "Richness",
  Topology = "Topology\ndistance",
  BiomassInteractions = "B-F-IS"
)

sem_parameters <- parameterEstimates(fit_sem)
n_vars <- length(sem_vars)
partition_list <- list()

for (g in seq_along(model_names)) {
  m <- model_names[g]
  param_g <- sem_parameters %>% filter(group == g)
  
  # All fitted paths, including nonsignificant paths
  reg_g <- param_g %>% filter(op == "~")
  cov_g <- param_g %>% filter(op == "~~")
  
  B <- matrix(0, n_vars, n_vars, dimnames = list(sem_vars, sem_vars))
  Psi <- matrix(0, n_vars, n_vars, dimnames = list(sem_vars, sem_vars))
  
  reg_positions <- cbind(match(reg_g$lhs, sem_vars),
                         match(reg_g$rhs, sem_vars))
  cov_positions <- cbind(match(cov_g$lhs, sem_vars),
                         match(cov_g$rhs, sem_vars))
  
  B[reg_positions] <- reg_g$est
  Psi[cov_positions] <- cov_g$est
  Psi[cov_positions[, c(2, 1), drop = FALSE]] <- cov_g$est
  
  # Propagate each source through all direct and indirect pathways
  total_paths <- solve(diag(n_vars) - B)
  dimnames(total_paths) <- list(sem_vars, sem_vars)
  
  Sigma <- total_paths %*% Psi %*% t(total_paths)
  
  for (resp in stability_responses) {
    
    variance_parts <- vapply(source_blocks, function(vars) {
      weights <- total_paths[resp, vars]
      as.numeric(t(weights) %*% Psi[vars, vars, drop = FALSE] %*% weights)
    }, numeric(1))
    
    fractions <- variance_parts / Sigma[resp, resp]
    full_r2 <- 1 - Psi[resp, resp] / Sigma[resp, resp]
    
    partition_list[[paste(m, resp, sep = "_")]] <- tibble(
      model = m, response = resp, component = names(component_labels),
      R2 = c(full_r2, unname(fractions))
    )
  }
}

partition_summary <- bind_rows(partition_list)


#------4. Combined figure: SEM paths and variation partitioning -----------------
# SEM path plots follow:
# https://cran.r-project.org/web/packages/diy.sem.plot/vignettes/diyPaths-vignette.html

node_list <- list(
  node("topo_distance",   x = -0.2, y =  2, label = "Topology\ndistance"),
  node("S_post",          x = -0.2, y = -2, label = "Richness"),
  node("biomass_shannon", x = 3.2, y =  0,   label = "Biomass\ndiversity"),
  node("gini_fluxes",     x = 7, y =  3.8,   label = "Flux Gini"),
  node("skewness_IS",     x = 7, y = -3.8,   label = "IS\nskewness"),
  node("resilience",      x = 10, y =  1.7, label = "Resilience"),
  node("reactivity",      x = 10, y = -1.7, label = "Reactivity")
)

path_list <- list(
  path("S_post", "topo_distance", side_from = "top", side_to = "bottom"),
  path("topo_distance", "biomass_shannon", side_from = "bottom", side_to = "left"),
  path("S_post", "biomass_shannon", side_from = "top", side_to = "left"),
  
  path("biomass_shannon", "gini_fluxes", side_from = "right", side_to = "left"),
  path("topo_distance", "gini_fluxes", side_from = "right", side_to = "left"),
  path("S_post", "gini_fluxes", side_from = "top", side_to = "left"),
  
  path("biomass_shannon", "skewness_IS", side_from = "right", side_to = "left"),
  path("topo_distance", "skewness_IS", side_from = "bottom", side_to = "left"),
  path("S_post", "skewness_IS", side_from = "right", side_to = "left"),
  
  path("gini_fluxes", "resilience", side_from = "right", side_to = "top"),
  path("skewness_IS", "resilience", side_from = "top", side_to = "left",
       nudge_text_x = 0.3, nudge_text_y = 0.55),
  path("biomass_shannon", "resilience", side_from = "right", side_to = "bottom"),
  path("topo_distance", "resilience", side_from = "right", side_to = "left",
       nudge_text_x = -0.9),
  path("S_post", "resilience", side_from = "right", side_to = "bottom"),
  
  path("gini_fluxes", "reactivity", side_from = "bottom", side_to = "left",
       nudge_text_x = -0.3, nudge_text_y = 0.55),
  path("skewness_IS", "reactivity", side_from = "right", side_to = "bottom"),
  path("biomass_shannon", "reactivity", side_from = "right", side_to = "top"),
  path("topo_distance", "reactivity", side_from = "right", side_to = "top"),
  path("S_post", "reactivity", side_from = "right", side_to = "left",
       nudge_text_x = -0.9)
)

title_list <- lapply(seq_along(model_names), function(g) {
  panel_title(g, paste0("(", letters[g], ") ", model_names[g]))
})

# Use coefficients from the globally standardised data
sem_panels <- diyPaths(
  fit = fit_sem,
  node_positions = node_list,
  path_positions = path_list,
  panel_titles = title_list,
  panel_cols = 3,
  standardised = FALSE,
  digits = 2,
  est_stars = FALSE, show_variances = FALSE, show_group_labels = FALSE,
  node_width = 2.4, node_height = 1.1,
  observed_node_text_size = 3.8,
  path_text_size = 3.8,
  line_thickness = 0.65, arrow_size = 0.2,
  margin_x = 0.2, margin_y = 0.2
)

sem_plots <- list()

# change path color and line thickness
for (g in seq_along(model_names)) {
  p <- sem_panels[[g]]
  
  # Keep significant regression paths and label their direction
  for (k in seq_along(p$layers)) {
    d <- p$layers[[k]]$data
    
    # Thickness of node-box borders
    if (inherits(p$layers[[k]]$geom, "GeomRect")) {
      p$layers[[k]]$aes_params$linewidth <- 0.4
    }
    
    if ("pvalue" %in% names(d)) {
      d <- d %>%
        filter(type == "~", pvalue < 0.05) %>%
        mutate(direction = if_else(est > 0, "Positive", "Negative"),
               width_group = case_when(
                   abs(est) < 0.2 ~ "0-0.2",
                   abs(est) < 0.5 ~ "0.2-0.5",
                   abs(est) < 1 ~ "0.5-1.0",
                   TRUE ~ "1.0+"
                 )
               )
      
      p$layers[[k]]$data <- d
      p$layers[[k]]$mapping$colour <- aes(colour = direction)$colour
      p$layers[[k]]$show.legend <- FALSE
      
      if (inherits(p$layers[[k]]$geom, "GeomSegment")) {
        p$layers[[k]]$aes_params$linewidth <- NULL
        p$layers[[k]]$mapping$linewidth <- aes(linewidth = width_group)$linewidth
        p$layers[[k]]$show.legend <- TRUE
        
        p$layers[[k]]$geom_params$arrow <- grid::arrow(
          length = grid::unit(0.2, "cm"),
          type = "open",
          angle = 30,
          ends = "last"
        )
      }
    }
  }
  
  p <- p +
    scale_colour_manual(
      name = "Path direction",
      values = c(
        Positive = col_div[3],
        Negative = col_div[1]
      ),
      breaks = c("Positive", "Negative"),
      limits = c("Positive", "Negative"),
      drop = FALSE
    ) +
    scale_linewidth_manual(
      name = "Path coefficient",
      values = c("0-0.2" = 0.3, "0.2-0.5" = 0.8,
                 "0.5-1.0" = 1.5, "1.0+" = 2.5)
 #     guide = "none"
    ) + 
    figure_theme +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
  
  # Keep as a ggplot so patchwork can collect its legend
  sem_plots[[g]] <- p
}



# Heatmaps of SEM variance sources
response_names <- c(resilience = "Resilience", reactivity = "Reactivity")

partition_plot_df <- partition_summary %>%
  mutate(model = factor(model, levels = model_order),
         component = factor(component, levels = rev(names(component_labels))),
         label = sprintf("%.2f", R2))

fill_limits <- range(c(0, partition_plot_df$R2))
fill_max <- max(abs(fill_limits))

partition_plot_df <- partition_plot_df %>%
  mutate(text_colour = if_else(abs(R2) > 0.65 * fill_max, "white", "black"))

partition_plots <- list()

for (j in seq_along(stability_responses)) {
  resp <- stability_responses[j]
  df_j <- partition_plot_df %>% filter(response == resp)
  
  p <- ggplot(df_j, aes(model, component, fill = R2)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = label, colour = text_colour),
              size = 4.5, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = col_div[1], mid = col_div[2], high = col_div[3],
      midpoint = 0, limits = fill_limits,
      labels = scales::label_number(accuracy = 0.01),
      guide = guide_colourbar(direction = "horizontal")
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0.5)) +
    scale_y_discrete(labels = component_labels, expand = expansion(add = 0.5)) +
    labs(title = paste0("(", letters[length(model_names) + j], ") ",
                        response_names[[resp]]), x = NULL, y = NULL) +
    figure_theme +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.margin = margin(2, 2, 2, 2))
  
  partition_plots[[j]] <- wrap_elements(full = p)
}

# Extract the shared legend from an SEM plot
legend <- cowplot::get_legend(
  sem_plots[[1]] +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(0, 0, 0, 0)
    )
)

# Remove legends from individual SEM panels
sem_plots_clean <- lapply(sem_plots, function(p) {
  p + theme(legend.position = "none")
})

# Put the legend in its own row
legend_row <- patchwork::wrap_elements(full = legend)

combined_figure <- wrap_plots(c(sem_plots_clean, partition_plots),
                              ncol = 3,
                              widths = c(1, 1, 1)) /
  legend_row +
  plot_layout(heights = c(1, 0.06))

combined_figure

ggsave("../figures/SEM_and_variation_partitioning.png", combined_figure,
       width = 16, height = 14, dpi = 600)



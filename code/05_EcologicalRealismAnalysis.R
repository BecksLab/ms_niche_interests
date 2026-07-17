library(tidyverse)
library(patchwork)
library(ggforce)
library(ggalluvial)

# ============================================================
# Realism classification
# ============================================================

classify_foodwebs <- function(dat) {
  
  dat %>%
    glow_up(
      illogical = if_else(is.na(illogical), 0, illogical),
      mx_tl = if_else(is.na(mx_tl), S, mx_tl),
      
      expected_mx_tl = 2 + 0.8 * log2(S),
      
      # Structural realism coordinate
      mx_tl_pct = mx_tl/expected_mx_tl,
      
      # Structural violations
      viol_mx_tl       = mx_tl > expected_mx_tl,
      viol_connectance = connectance < 0.01 | connectance > 0.40,
      
      # Species violations
      viol_isolated     = isolated >= 1,
      viol_illogical    = illogical >= 1,
      
      structure_fail =
        viol_mx_tl |
        viol_connectance,
      
      species_fail =
        viol_isolated |
        viol_illogical,
      
      category = case_when(
        !structure_fail & !species_fail ~ "pass",
        structure_fail & !species_fail  ~ "illogical_struct",
        !structure_fail & species_fail  ~ "illogical_spp",
        structure_fail & species_fail   ~ "illogical_both"
      ),
      
      # Inside structural realism region
      inside_structure =
        connectance >= 0.01 &
        connectance <= 0.40 &
        mx_tl_pct <= 1,
      
      # Inside species realism region
      inside_species =
        isolated == 0 &
        illogical == 0
    )
}



# ============================================================
# Prepare data
# ============================================================

realism_post <-
  read.csv("outputs/realism_END.csv") %>%
  yeet(S != 0) %>%
  classify_foodwebs() %>%
  mutate(
    prop_isolated = isolated / S,
    prop_illogical = illogical / S,
    state = "post"
  )

realism_pre <-
  read.csv("outputs/realism_initial.csv") %>%
  classify_foodwebs() %>%
  mutate(
    prop_isolated = isolated / S,
    prop_illogical = illogical / S,
    state = "pre")

realism_df <- 
  rbind(realism_post, realism_pre) %>%
  glow_up(mx_tl_pct = if_else(mx_tl < 0,
                              NA,
                              mx_tl_pct),
          state = factor(state, levels = c("pre", "post")))

# ============================================================
# Structural realism plots
# ============================================================

p_co <-
  ggplot(realism_df,
         aes(y = connectance,
             x = state,
             colour = model)) + 
  annotate("rect", 
           fill = "#78BE21", 
           alpha = 0.1, 
           xmin = -Inf, xmax = Inf,
           ymin = 0.01, ymax = 0.4) +
  geom_sina(position = "dodge", 
            alpha = 0.5) +
  labs(x = NULL,
       y = NULL,
       title = "Connectance") +
  scale_colour_manual(values = model_colours) +
  figure_theme +
  theme(legend.position = 'none')


p_tlmax <-
  ggplot(realism_df %>%
           yeet(mx_tl_pct < 2),
         aes(y = mx_tl_pct,
             x = state,
             colour = model)) + 
  annotate("rect", 
           fill = "#78BE21", 
           alpha = 0.1, 
           xmin = -Inf, xmax = Inf,
           ymin = 0, ymax = 1) +
  geom_sina(position = "dodge", 
            alpha = 0.5) +
  scale_colour_manual(values = model_colours) +
  labs(title = "TL as % of expected max",
       x = NULL,
       y = NULL,
       caption = "Note networks with negative TL were not reported, expected max = 2 + 0.8 * log2(S)") +
  figure_theme +
  theme(legend.position = 'none')

p_survive  <- 
  realism_df %>%
  vibe_check(model, fw_ID, state, S) %>%
  pivot_wider(names_from = state,
              values_from = S) %>%
  glow_up(spp_loss = (pre-post)/pre) %>%
  ggplot(aes(y = spp_loss,
             x = model,
             colour = model)) + 
  annotate("rect", 
           fill = "#78BE21", 
           alpha = 0.1, 
           xmin = -Inf, xmax = Inf,
           ymin = 0, ymax = 0.2) +
  geom_sina(position = "dodge", 
            alpha = 0.5) +
  scale_colour_manual(values = model_colours) +
  labs(title = "Percent species lost at equilibrium",
       x = NULL,
       y = NULL) +
  figure_theme +
  theme(legend.position = 'none')

p_interval  <- 
  realism_df %>%
  vibe_check(model, state, intervality) %>%
  na.omit() %>%
  squad_up(model, state) %>%
  no_cap(prop_nonzero = mean(intervality != 0),
         n = dplyr::n()) %>%
  ggplot(aes(y = prop_nonzero,
             x = state,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  labs(title = "Proportion non-interval networks",
       x = NULL,
       y = NULL) +
  figure_theme +
  theme(legend.position = 'none')


# ============================================================
# Species realism plots
# ============================================================

realism_df_freq <-
  realism_df %>%
  squad_up(model, state) %>%
  no_cap(prop_isolated = sum(viol_isolated, na.omit = TRUE)/n(),
         prop_illogical = sum(viol_illogical, na.omit = TRUE)/n()) %>%
  glow_up(prop_illogical = if_else(state == "post",
                                   0,
                                   prop_illogical))

p_iso <-
  ggplot(realism_df_freq,
         aes(y = prop_isolated,
             x = state,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  labs(title = "Proportion isolated species",
       x = NULL,
       y = NULL) +
  figure_theme

p_illog <-
  ggplot(realism_df %>%
           glow_up(illog_prop = illogical/S) %>%
           squad_up(model, state) %>%
           no_cap(illog_prop = mean(illog_prop)),
         aes(y = illog_prop,
             x = state,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  labs(title = "Mean proportion of illogical species per network (illogical/S)",
       x = NULL,
       y = NULL) +
  figure_theme

p_stable <-
  read.csv("outputs/realism_END.csv") %>%
  left_join(read.csv("outputs/realism_initial.csv") %>%
              glow_up(richness_init = case_when(S == 10 ~ "S_init = 10",
                                                S == 15 ~ "S_init = 15",
                                                S == 20 ~ "S_init = 20")) %>%
              vibe_check(fw_ID, richness_init)) %>%
  glow_up(stab_state = if_else(is.na(S),
                               1,
                               0)) %>%
  squad_up(model, richness_init) %>%
  no_cap(prop_failed = sum(stab_state)/n()) %>%
  ggplot(aes(y = prop_failed,
             x = model,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  facet_wrap(vars(richness_init)) +
  labs(title = "Proportion of networks not reaching stability",
       x = NULL,
       y = NULL) +
  figure_theme


# ============================================================
# Assemble figure
# ============================================================

(
  p_co +
    p_tlmax +
    p_survive
) /
  (
    p_interval +
      theme(legend.position = 'bottom') +
      p_iso +
      theme(legend.position = 'bottom') +
      p_illog +
      theme(legend.position = 'bottom')
  ) +
  plot_layout(guides='collect') +
  plot_annotation(title = "Food-web realism space",
                  theme = theme(legend.position = 'bottom'))

ggsave("../figures/realism_space.png",
       dpi = 600,
       width = 10000,
       height = 7000,
       units = "px")

pre_summary  <-  
  read.csv("outputs/realism_initial.csv") %>%
  classify_foodwebs() %>%
  vibe_check(model, category, fw_ID) %>%
  lowkey(pre_state = category)
post_summary <- 
  read.csv("outputs/realism_END.csv") %>%
  glow_up(illogical = if_else(illogical > 0,
                              0,
                              illogical)) %>%
  classify_foodwebs() %>%
  lowkey(post_state = category) %>%
  glow_up(post_state = if_else(is.na(S),
                               "non_stable",
                               post_state)) %>%
  vibe_check(model, post_state, fw_ID)
stab_summ <-  
  read.csv("outputs/realism_END.csv") %>%
  classify_foodwebs() %>%
  glow_up(stab_state = if_else(is.na(S),
                               "non_stability",
                               "stable"))  %>%
  vibe_check(model, stab_state, fw_ID)

summary_wide <- 
  left_join(pre_summary, post_summary) %>%
  left_join(stab_summ) %>%
  squad_up(model, post_state, pre_state, stab_state) %>%
  tally() %>%
  glow_up(alpha = if_else(post_state == "pass" &&
                            #pre_state == "pass" &&
                            stab_state == "stable",
                          0.8,
                          0.3))

ggplot(data = summary_wide,
       aes(axis1 = pre_state, axis2 = stab_state, axis3 = post_state,
           y = n)) +
  scale_x_discrete(limits = c("Pre Struct", "Stability", "Post Struct"), 
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = model, 
                    alpha = alpha),
                knot.pos = 0.5,
                colour = "white") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_wrap(vars(model)) +
  scale_fill_manual(values = model_colours) +
  scale_alpha_identity() +
  figure_theme +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("../figures/realism_space_sankey.png",
       dpi = 600,
       width = 8000,
       height = 6000,
       units = "px")

pre_summary  <-  
  read.csv("outputs/realism_initial.csv") %>%
  classify_foodwebs() %>%
  vibe_check(model, fw_ID, viol_mx_tl, viol_connectance, viol_isolated,
             viol_illogical)
stab_summ <-  
  read.csv("outputs/realism_END.csv") %>%
  glow_up(stab_state = if_else(is.na(S),
                               "non_stability",
                               "stable"))  %>%
  vibe_check(model, stab_state, fw_ID)

summary_wide <- 
  left_join(pre_summary, stab_summ) %>%
  vibe_check(-fw_ID) %>%
  pivot_longer(-c(model, stab_state)) %>%
  squad_up(model, stab_state, name, value) %>%
  tally() %>%
  glow_up(alpha = if_else(value == TRUE,
                          0.8,
                          0.3))

ggplot(data = summary_wide,
       aes(axis1 = model, axis2 = value, axis3 = stab_state,
           y = n)) +
  scale_x_discrete(limits = c("Model", "Violated?", "Stability"), 
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = model, 
                    alpha = alpha),
                knot.pos = 0.5,
                colour = "white") +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_wrap(vars(name)) +
  scale_fill_manual(values = model_colours) +
  scale_alpha_identity() +
  labs(y = NULL) +
  figure_theme +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("../figures/realism_space_sankey_metric.png",
       dpi = 600,
       width = 8000,
       height = 6000,
       units = "px")

# also export list of 'perfect' network ids
left_join(pre_summary, post_summary) %>%
  left_join(stab_summ) %>%
  yeet(post_state == "pass" &
         #pre_state == "pass" &
         stab_state == "stable") %>%
  vibe_check(fw_ID) %>%
  write_csv("outputs/perfect_net_ids.csv")

library(ggridges)

realism_df %>%
  vibe_check(model, intervality, loops, mx_tl) %>%
  pivot_longer(-model) %>%
  ggplot(aes(x = value, 
             y = model, 
             fill = model)) +
  geom_density_ridges() +
  scale_fill_manual(values = model_colours) +
  facet_wrap(vars(name),
             scales = "free") + 
  theme(legend.position = "none")




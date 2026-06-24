library(tidyverse)
library(patchwork)

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
      viol_prop_basal  = prop_basal > 0.50,
      
      # Species violations
      viol_isolated     = isolated >= 1,
      viol_closed_loops = closed_loops >= 1,
      viol_illogical    = illogical >= 1,
      
      structure_fail =
        viol_mx_tl |
        viol_connectance |
        viol_prop_basal,
      
      species_fail =
        viol_isolated |
        viol_closed_loops |
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
        prop_basal <= 0.50 &   
        prop_basal >= 0.05 &
        mx_tl_pct <= 1,
      
      # Inside species realism region
      inside_species =
        isolated == 0 &
        closed_loops == 0 &
        illogical == 0
    )
}



# ============================================================
# Prepare data
# ============================================================

realism_post <-
  post_df %>%
  yeet(!is.na(S)) %>%
  classify_foodwebs() %>%
  mutate(
    prop_isolated = isolated / S,
    prop_loops    = closed_loops / S^2,
    prop_illogical = illogical / S,
    state = "post"
  )

realism_pre <-
  pre_df %>%
  classify_foodwebs() %>%
  mutate(
    prop_isolated = isolated / S,
    prop_loops    = closed_loops / S^2,
    prop_illogical = illogical / S,
    state = "pre")

realism_df <- rbind(realism_post, realism_pre) %>%
  glow_up(mx_tl_pct = if_else(closed_loops > 0,
                              NA,
                              mx_tl_pct))

# ============================================================
# Structural realism plots
# ============================================================

p_struct_1 <-
  ggplot(realism_df,
         aes(y = connectance,
             x = state,
             colour = model)) +
  geom_hline(yintercept = 0.4, 
             colour = shark_silver) +
  geom_hline(yintercept = 0.01, 
             colour = shark_silver) +
  geom_violin(position = "dodge", 
              alpha = 0.5) +
  scale_colour_manual(values = model_colours) +
  figure_theme


p_struct_2 <-
  ggplot(realism_df,
         aes(y = mx_tl_pct,
             x = state,
             colour = model)) +
  geom_hline(yintercept = 1, 
             colour = shark_silver) +
  geom_violin(position = "dodge", 
              alpha = 0.5) +
  scale_colour_manual(values = model_colours) +
  labs(y = "TL as % of expected max",
       caption = "Note for networks with loops TL was not reported, expected max = 2 + 0.8 * log2(S)") +
  figure_theme

p_struct_3 <-
  ggplot(realism_df,
         aes(y = prop_basal,
             x = state,
             colour = model)) +
  geom_hline(yintercept = 0.5, 
             colour = shark_silver) +
  geom_hline(yintercept = 0.05, 
             colour = shark_silver) +
  geom_violin(position = "dodge", 
              alpha = 0.5) +
  scale_colour_manual(values = model_colours) +
  labs(y = "Proportion basal") +
  figure_theme

# ============================================================
# Species realism plots
# ============================================================

realism_df_freq <-
  realism_df %>%
  squad_up(model, state) %>%
  no_cap(prop_isolated = sum(viol_isolated, na.omit = TRUE)/n(),
         prop_illogical = sum(viol_illogical, na.omit = TRUE)/n(),
         prop_loops = sum(viol_closed_loops, na.omit = TRUE)/n())

p_species_1 <-
  ggplot(realism_df_freq,
         aes(y = prop_isolated,
             x = state,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  labs(y = "Proportion isolated species") +
  figure_theme

p_species_2 <-
  ggplot(realism_df_freq,
         aes(y = prop_illogical,
             x = state,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  labs(y = "Proportion illogical species") +
  figure_theme

p_species_3 <-
  ggplot(realism_df_freq,
         aes(y = prop_loops,
             x = state,
             fill = model)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = model_colours) +
  labs(y = "Proportion of loops",
       caption = "calculated as n loops/S^1") +
  figure_theme

# ============================================================
# Assemble figure
# ============================================================

(
  p_struct_1 +
    p_struct_2 +
    p_struct_3
) /
  (
    p_species_1 +
      p_species_2 +
      p_species_3
  ) +
  plot_layout(guides='collect') +
  plot_annotation(
    title = "Food-web realism space",
    subtitle = "Shaded regions denote empirically realistic structural bounds"
  )

library(scatterplot3d)

scatterplot3d(
  x = realism_df$connectance,
  y = realism_df$prop_basal,
  z = realism_df$mx_tl_pct,
  color = model_colours[realism_df$model],
  pch = if_else(realism_df$state == "pre", 16, 1)
)


pre_summary  <- classify_foodwebs(pre_df) %>%
  vibe_check(model, category, fw_ID) %>%
  lowkey(pre_state = category)
post_summary <- classify_foodwebs(post_df) %>%
  lowkey(post_state = category) %>%
  glow_up(post_state = if_else(is.na(S),
                               "non_stable",
                               post_state)) %>%
  vibe_check(model, post_state, fw_ID)
stab_summ <- post_df %>%
  glow_up(stab_state = if_else(is.na(S),
                               "non_stability",
                               "stable"))  %>%
  vibe_check(model, stab_state, fw_ID)

summary_wide <- 
  left_join(pre_summary, post_summary) %>%
  left_join(stab_summ) %>%
  squad_up(model, post_state, pre_state, stab_state) %>%
  tally()

ggplot(data = summary_wide,
       aes(axis1 = pre_state, axis2 = stab_state, axis3 = post_state,
           y = n)) +
  scale_x_discrete(limits = c("Pre Struct", "Stability", "Post Struct"), 
                   expand = c(.2, .05)) +
  geom_alluvium(aes(fill = model)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_wrap(vars(model)) +
  scale_fill_manual(values = model_colours) +
  figure_theme 


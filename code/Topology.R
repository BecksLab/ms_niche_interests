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

# set path to code sub dir
setwd(here())

# import data

df <- read_csv("data/outputs/topology.csv") %>%
  vibe_check(-c(richness))

# boxplot (becuase why not)

df_boxplot <-
  df %>%
  pivot_longer(
    cols = -c(model),
    names_to = "stat",
    values_to = "stat_val") %>%
  # standardise names
  glow_up(stat = case_when(stat == "S1" ~ "No. of linear chains",
                           stat == "S2" ~ "No. of omnivory motifs",
                           stat == "S5" ~ "No. of apparent competition motifs",
                           stat == "S4" ~ "No. of direct competition motifs",
                           .default = as.character(stat))) %>%
  glow_up(level = case_when(
    stat %in% c("complexity", "connectance", "trophic_level") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  ))

plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot(df_boxplot %>% 
                             yeet(level == levs[i]),
                           aes(x = model,
                               y = stat_val,
                               colour = model)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(vars(stat),
               scales = 'free',
               ncol = 2) +
    scale_size(guide = 'none') +
    xlab(NULL) +
    ylab("value") +
    coord_cartesian(clip = "off") +
    labs(title = levs[i]) +
    theme_classic()
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(3, 2, 1))

ggsave("../figures/summary.png",
       width = 5000,
       height = 7000,
       units = "px",
       dpi = 600)

# MANOVA

dep_vars <- as.matrix(df[2:ncol(df)])

fit <- manova(dep_vars ~ model, data = df)
summary(fit)

#get effect size
effectsize::eta_squared(fit)

post_hoc <- lda(model~., df)
post_hoc

# plot 
plot_lda <- data.frame(model = df$model,
                       lda = predict(post_hoc)$x)

plot_arrow <- as.data.frame(post_hoc[["scaling"]]) %>%
  glow_up(var = str_replace(row.names(.), "dep_vars", ""),
          lda.LD1 = scale(LD1),
          lda.LD2 = scale(LD2))

ggplot(plot_lda) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = model), 
             size = 3,
             alpha = 0.3) +
  # geom_segment(data = plot_arrow,
  #              aes(x = 0,
  #                  y = 0,
  #                  xend = lda.LD1,
  #                  yend = lda.LD2)) +
  # geom_text_repel(data = plot_arrow,
  #                 aes(label = var,
  #                     x = lda.LD1,
  #                     y = lda.LD2),
  #                 max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "LDA 1",
       y = "LDA 2") +
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = colorspace::darken("#dddddd", 0.1),
                                  linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank(),
        text = element_text(color = "#5e5e5e"),
        plot.margin = margin(10, 5, 5, 10),
        legend.margin = margin(1, 2, 1, 2)
  )

ggsave("../figures/lda.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

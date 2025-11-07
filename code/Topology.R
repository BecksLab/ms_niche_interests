# libraries
library(ggpubr)
library(here)
library(effectsize)
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
  dplyr::select(-c(richness))

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
  mutate(var = str_replace(row.names(.), "dep_vars", ""),
         lda.LD1 = scale(LD1),
         lda.LD2 = scale(LD2))

ggplot(plot_lda) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = model), 
             size = 3,
             alpha = 0.3) +
  #=  geom_segment(data = plot_arrow,
  #=               aes(x = 0,
  #=                   y = 0,
  #=                   xend = lda.LD1,
  #=                   yend = lda.LD2)) +
  #=  geom_text_repel(data = plot_arrow,
  #=                  aes(label = var,
  #=                      x = lda.LD1,
  #=                      y = lda.LD2),
  #=                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
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

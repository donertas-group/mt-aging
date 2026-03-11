# main figure 7 panel a
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

margins = margin(5.5,5.5,5.5,5.5, "pt")  # t r b l 
significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

source('scripts/files.R')

pa =
  padat %>% 
  ggplot(aes(x=abv3q, y=rho)) +
  geom_violin(fill = "grey90") + 
  geom_boxplot(width = .2) +
  ggforce::geom_sina(method = "counts", alpha = .5) +
  geom_signif(test="wilcox.test", comparisons = list(c("Lowest 75%", "Highest 25%")), map_signif_level = significances) +
  theme_bw() +
  ylim(c(-0.34, 0.35)) +
  labs(x='', y='Age-expression correlation') +
  plot_annotation(title='(a)') +
  theme(plot.margin = margins, axis.text.x = element_text(size=10))

padat %>% 
  group_by(abv3q) %>% 
  count()

saveRDS(pa, file.path(figure_out_dir, 'figure7', 'figure7_a.rds'))

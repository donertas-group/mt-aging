# main figure 7 panel b
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

margins = margin(5.5,5.5,5.5,5.5, "pt")  # t r b l 
significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

source('scripts/files.R')

pb = pbdat %>% 
  ggplot(aes(x=age, y=expression)) +
  geom_violin(aes(fill=avgexp)) +
  geom_boxplot(width = .2, alpha=1, outliers=F) +
  scale_fill_gradient(name='Mean\nnormalized\nexpression', low='#FCFBFD', high='#6A51A3') +
  ggforce::geom_sina(method = "counts", alpha = .2) +
  theme_bw() +
  labs(x='Age (in years)', y='Normalized expression level')+
  plot_annotation(title='(b)') +
  theme(plot.margin = margins)

saveRDS(pb, file.path(figure_out_dir, 'figure7', 'figure7_b.rds'))

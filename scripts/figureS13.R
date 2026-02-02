library(tidyverse)
library(reshape2)
library(ggsignif)
library(patchwork)

source('scripts/files.R')

pa = psfig13a %>% 
  ggplot(aes(x=age, y=expression)) +
  geom_violin(aes(fill=avgexp)) +
  geom_boxplot(width = .2, alpha=1, outliers=F) +
  scale_fill_gradient(name='Mean\nnormalized\nexpression', low='#DADAEB', high='#6A51A3') +
  ggforce::geom_sina(method = "counts", alpha = .2) +
  theme_bw() +
  labs(x='Age (in years)', y='Normalized expression level') +
  plot_annotation(title='(a)', subtitle='\tAdipose - visceral')


pb = psfig13b %>%
  ggplot(aes(x=age, y=expression)) +
  geom_violin(aes(fill=avgexp)) +
  geom_boxplot(width = .2, alpha=1, outliers=F) +
  scale_fill_gradient(name='Mean\nnormalized\nexpression', low='#DADAEB', high='#6A51A3') +
  ggforce::geom_sina(method = "counts", alpha = .2) +
  theme_bw() +
  labs(x='Age (in years)', y='Normalized expression level') +
  plot_annotation(title='(b)', subtitle='\tOvary')

p = ggarrange(pa, pb, ncol=1)
ggsave(file.path(figure_out_dir, 'figureS13.pdf'), p, width=6, height=9)

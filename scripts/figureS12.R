library(tidyverse)
library(reshape2)
library(ggsignif)
library(patchwork)

source('scripts/files.R')

p =
  psfig12dat %>% 
  mutate(rho_range = cut(rho, 
                         breaks = c(-Inf, -0.3, -0.2, -0.1, 0, 0.1, 0.2, Inf),
                         labels = c('-0.3', '-0.2', '-0.1', '0', '0.1', '0.2', '0.3'),
                         right = FALSE)) %>% 
  ggplot(aes(y=reorder(minor_tissue, expression), x=expression)) +
  geom_violin(aes(fill=rho_range), trim = TRUE, scale = 'width', na.rm = FALSE, alpha=.8) +
  # geom_hline(yintercept ='Skin - Sun Exposed (Lower leg)', vjust=1.5) +
  geom_hline(yintercept = 34.5, linetype=2, alpha=.5, color='gray50') +
  annotate("text", x = 400, y=37, label = "Highest 25%", size=6) +
  scale_fill_manual(name='Age-related\nexpression change',
                    values=RColorBrewer::brewer.pal(8, 'RdBu')[1:7]) +
  geom_boxplot(outlier.colour = NULL, outlier.fill = NULL, outlier.shape = NA, outlier.stroke = 0.5, width = 0.2) +
  theme_bw() +
  scale_x_continuous(position = "top") +
  labs(x='TPM', y=' ')

ggsave(file.path(figure_out_dir, 'figureS12.pdf'), p, width=8, height=11)
  
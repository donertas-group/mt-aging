# combine all the panels of figure 7
library(tidyverse)
library(ggpubr)
library(patchwork)

source('scripts/files.R')

pa = readRDS(file.path(figure_out_dir, 'figure7', 'figure7_a.rds'))
pb = readRDS(file.path(figure_out_dir, 'figure7', 'figure7_b.rds'))
pc = readRDS(file.path(figure_out_dir, 'figure7', 'figure7_c.rds'))
pde = readRDS(file.path(figure_out_dir, 'figure7', 'figure7_de.rds'))
pfg = readRDS(file.path(figure_out_dir, 'figure7', 'figure7_fg.rds'))
ph = readRDS(file.path(figure_out_dir, 'figure7', 'figure7_h.rds'))

pab = ggarrange(pa, pb, nrow=1, widths = c(.8, 1))

pdefg = ggarrange(pde, pfg, ncol=1)
pcdefg = ggarrange(pc, pdefg, widths = c(.8, 1))
pabcdefg = ggarrange(pab, pcdefg, ncol=1, heights=c(1, 1.4))
pabcdefgh = ggarrange(pabcdefg, ph, ncol=1, heights = c(1.8, 1))

to_scale=.95
ggsave(file.path(figure_out_dir, 'figure7_final.pdf'), pabcdefgh, width = 210 * to_scale, height = 297 * to_scale, units = "mm")

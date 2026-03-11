# main figure 7 panels d-e
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

margins = margin(5.5,5.5,5.5,5.5, "pt")   # t r b l 
significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

source('scripts/files.R')

fields = c('sex', 'age_at_rec', 'PC', 'POFA_Total_Ratio', 'TFA', 'Lactate')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

mydf3=mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  rename(age=age_at_rec) %>% 
  mutate(lactate_quartile = ntile(Lactate, 4)) %>% 
  filter(lactate_quartile %in% c(1,4)) %>% 
  mutate(lactate_quartile = factor(lactate_quartile, levels=1:4)) %>% 
  transmute(lactate_quartile, PC_TFA_Ratio, PUFA_Total_Ratio=POFA_Total_Ratio) %>% 
  drop_na() 

pd1 =
  mydf3 %>% 
  ggplot(aes(x=lactate_quartile, y=PC_TFA_Ratio, fill=lactate_quartile)) +
  geom_boxplot(alpha=.9, outlier.size=.3) +
  geom_signif(test="wilcox.test", comparisons = list(c('4', '1')), map_signif_level = significances) +
  scale_y_continuous(expand=expand_scale(.25,0)) +
  scale_fill_manual(name='Lactate Quartile', 
                    breaks=c(1,4),
                    labels=c('Q1', 'Q4'),
                    values=RColorBrewer::brewer.pal(3, 'Set3')) +
  scale_x_discrete(breaks=c(1,4), 
                   labels=c('Q1', 'Q4')) +
  theme_bw() +
  theme(legend.position='none',
        plot.margin = margins) +
  xlab('Lactate Quartile') + ylab('PC/TFA Ratio')

pd2=
  mydf3 %>% 
  ggplot(aes(x=lactate_quartile, y=PUFA_Total_Ratio, fill=lactate_quartile)) +
  geom_boxplot(alpha=.9, outlier.size=.3) +
  geom_signif(test="wilcox.test", comparisons = list(c('4', '1')), map_signif_level = significances) +
  scale_y_continuous(expand=expand_scale(.25,0)) +
  scale_fill_manual(name='Lactate Quartile', 
                    breaks=c(1,4),
                    labels=c('Q1', 'Q4'),
                    values=RColorBrewer::brewer.pal(3, 'Set3')) +
  scale_x_discrete(breaks=c(1,4), 
                   labels=c('Q1', 'Q4')) +
  theme_bw() +
  theme(legend.position='none',
        plot.margin = margins) +
  xlab('Lactate Quartile') + ylab('PUFA/TFA Ratio')

pd = (pd1 | pd2) + plot_annotation(tag_levels = list(c('d', 'e')), tag_prefix='(', tag_suffix=')')
saveRDS(pd, file.path(figure_out_dir, 'figure7', 'figure7_de.rds'))

# calculate density estimates for boxplots
smydf3 = split(mydf3, as.character(mydf3$lactate_quartile))
mydf3_save = lapply(smydf3, function(df){
  data.frame('value.pc_tfa_ratio' = density(df$PC_TFA_Ratio)$x,
             'density_est.pc_tfa_ratio' = density(df$PC_TFA_Ratio)$y,
             
             'value.pufa_tfa_ratio' = density(df$PUFA_Total_Ratio)$x,
             'density_est.pufa_tfa_ratio' = density(df$PUFA_Total_Ratio)$y)
}) %>%
  melt(measure.vars=c()) %>% 
  rename(lactate_quartile=L1)
write.csv(mydf3_save, file=file.path(table_out_dir, 'f7d-e.csv'), row.names = F)

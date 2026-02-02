library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

margins = margin(5.5,0.5,5.5,5.5, "pt")  # t r b l 
significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

source('scripts/files.R')

fields = c('sex', 'age_at_rec', 'PC', 'TFA')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

mydf3 = mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  select(-TFA) %>% 
  mutate(age = as.numeric(age_at_rec)) %>%
  drop_na() %>% 
  select(-age_at_rec) %>% 
  rename(Sex=sex)

pc1 = mydf3 %>% 
  ggplot(aes(x=age, y=PC, color=Sex)) +
  geom_smooth(se=F) +
  geom_quantile(quantiles=c(.25, .75), formula=y ~ poly(x, 3), linetype=2) +
  # geom_quantile(quantiles=c(.5), formula=y ~ poly(x, 3), colour="red") +
  # facet_wrap(~variable, scales='free_y', ncol=5) + 
  ylab('PC level') + xlab('Age (in years)') +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  # ggtitle('Age-associated changes') +
  theme_bw() +
  theme(plot.margin = margins)

pc2 = mydf3 %>% 
  ggplot(aes(x=age, y=PC_TFA_Ratio, color=Sex)) +
  geom_smooth(se=F) +
  geom_quantile(quantiles=c(.25, .75), formula=y ~ poly(x, 3), linetype=2) +
  # geom_quantile(quantiles=c(.5), formula=y ~ poly(x, 3), colour="red") +
  # facet_wrap(~variable, scales='free_y', ncol=5) + 
  ylab('PC/TFA Ratio') + xlab('Age (in years)') +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  # ggtitle('Age-associated changes') +
  theme_bw() +
  theme(plot.margin = margins)

# pc = ggarrange(pc1, pc2, ncol=1, heights = c(1, 1), common.legend=T, legend = 'right')
pc = pc1 / pc2 + plot_annotation(title='(c)') + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

saveRDS(pc, file.path(figure_out_dir, 'figure7', 'figure7_c.rds'))


# get loess preds
smydf3 = split(mydf3, mydf3$Sex)
mydf3_save = lapply(smydf3, function(df){
  df$pc.loess_preds = stats::predict(stats::loess(PC ~ age, data=df))
  df$pc_tfa_ratio.loess_preds = stats::predict(stats::loess(PC_TFA_Ratio ~ age, data=df))
  df
}) %>% 
  bind_rows() %>% 
  select(-PC, -PC_TFA_Ratio)
write.csv(mydf3_save, file=file.path(table_out_dir, 'f7c.csv'), row.names = F)
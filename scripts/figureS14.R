library(tidyverse)
library(ggpubr)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

# panel a - age dist
cols = mymetadf %>% 
  filter(name %in% c('sex', 'age_at_rec')) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

pa = mydf2 %>% 
  mutate(ageints = cut(as.numeric(age_at_rec), breaks = seq(30,80,by=5))) %>% 
  group_by(sex, ageints) %>% 
  summarise(n=n()) %>% 
  ggplot(aes(x=ageints, y=n, fill=sex)) +
  geom_bar(stat='identity', alpha=.7) +
  geom_text(aes(label=n), vjust=-.4) +
  # scale_y_continuous(expand=expand_scale(.25,0)) +
  facet_grid(sex~.) +
  theme_bw() +
  ylim(c(0, 4500)) +
  theme(legend.position='none') +
  xlab('Age interval') +
  ylab('Num. of individuals')
write.csv(select(pa$data, c(sex, ageints, n)), file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/fs14a.csv', row.names = F)

# panel b - corrs
rm(list=setdiff(ls(), c('mydf', 'mymetadf', 'pa')))
fields = c('sex', 'age_at_rec', 'POFA', 'MOFA', 'SFA', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio', 'TFA')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) %>% drop_na() 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

mydf3 = mydf2 %>% 
  melt(id.vars=c('age_at_rec', 'sex')) %>% 
  mutate(age = as.numeric(age_at_rec)) %>%
  drop_na() %>% 
  mutate(variable = gsub('POFA', 'PUFA', variable),
         variable = gsub('MOFA', 'MUFA', variable),
         variable = gsub('Total', 'TFA', variable)) %>% 
  mutate(variable = setNames(c('TFA', 'PUFA', 'MUFA', 'SFA', 'PUFA/MUFA Ratio', 'PUFA/TFA Ratio', 'MUFA/TFA Ratio', 'SFA/TFA Ratio'),
                             c('TFA', 'PUFA', 'MUFA', 'SFA', 'PUFA_MUFA_Ratio', 'PUFA_TFA_Ratio', 'MUFA_TFA_Ratio', 'SFA_TFA_Ratio'))[variable]) 

pb = mydf3 %>% 
  mutate(variable = factor(variable, levels=c('TFA', 'PUFA', 'MUFA', 'SFA', 'PUFA/MUFA Ratio', 'PUFA/TFA Ratio', 'MUFA/TFA Ratio', 'SFA/TFA Ratio'))) %>% 
  ggplot(aes(x=age, y=value, color=sex)) +
  geom_smooth(se=F) +
  geom_quantile(quantiles=c(.25, .75), formula=y ~ poly(x, 3), linetype=2) +
  # geom_quantile(quantiles=c(.5), formula=y ~ poly(x, 3), colour="red") +
  facet_wrap(~variable, scales='free_y', ncol=4, strip.position = "left") + 
  labs(y='', x='Age (in years)', color='Sex') +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  theme_bw() +
  theme(legend.position = 'bottom',
        strip.background = element_blank(),
        strip.placement = "outside")

# get loess preds
smydf2 = split(mydf2, mydf2$sex)
mydf3_save = lapply(smydf2, function(df){
  df$tfa.loess_preds = stats::predict(stats::loess(TFA ~ age_at_rec, data=df))
  df$pufa.loess_preds = stats::predict(stats::loess(POFA ~ age_at_rec, data=df))
  df$mufa.loess_preds = stats::predict(stats::loess(MOFA ~ age_at_rec, data=df))
  df$sfa.loess_preds = stats::predict(stats::loess(SFA ~ age_at_rec, data=df))
  df$pufa_mufa_ratio.loess_preds = stats::predict(stats::loess(POFA_MOFA_Ratio ~ age_at_rec, data=df))
  df$pufa_tfa_ratio.loess_preds = stats::predict(stats::loess(POFA_Total_Ratio ~ age_at_rec, data=df))
  df$mufa_tfa_ratio.loess_preds = stats::predict(stats::loess(MOFA_Total_Ratio ~ age_at_rec, data=df))
  df$sfa_tfa_ratio.loess_preds = stats::predict(stats::loess(SFA_Total_Ratio ~ age_at_rec, data=df))
  df
}) %>% 
  bind_rows() %>% 
  select(-c(TFA, POFA, MOFA, SFA, POFA_MOFA_Ratio, POFA_Total_Ratio, MOFA_Total_Ratio, SFA_Total_Ratio)) %>% 
  rename(age=age_at_rec)

write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/fs14b.csv', row.names = F)

p = ggarrange(pa, pb, nrow=2, heights = c(1, 1.2))
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/figureS14.pdf', p, width=8, height=6)

library(tidyverse)
library(reshape2)
library(ggsignif)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

# panel a
fields = c('age_at_rec', 'POFA_MOFA_Ratio', 'SFA_Total_Ratio', 'MOFA_Total_Ratio', 'Lactate')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

mydf3 = mydf2 %>% 
  mutate(lactate_quartile = ntile(Lactate, 4)) %>% 
  filter(lactate_quartile %in% c(1,4)) %>% 
  mutate(lactate_quartile = factor(lactate_quartile, levels=1:4)) %>% 
  select(-Lactate, -age_at_rec) %>% 
  drop_na()
  
pa = 
  mydf3 %>% 
  melt(id.vars=c('lactate_quartile')) %>%
  drop_na() %>% 
  mutate(variable = setNames(c('PUFA/MUFA Ratio', 'SFA/TFA Ratio', 'MUFA/TFA Ratio'),
                             c('POFA_MOFA_Ratio', 'SFA_Total_Ratio', 'MOFA_Total_Ratio'))[variable]) %>% 
  mutate(variable = factor(variable, levels=c('PUFA/MUFA Ratio', 'MUFA/TFA Ratio', 'SFA/TFA Ratio'))) %>% 
  ggplot(aes(x=lactate_quartile, y=value, fill=lactate_quartile)) +
  geom_boxplot(alpha=.9, outlier.size=.3) +
  geom_signif(test="wilcox.test", comparisons = list(c('4', '1')), map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.25,0)) +
  scale_fill_manual(name='Lactate Quartile', 
                    breaks=c(1,4),
                    labels=c('Q1', 'Q4'),
                    values=RColorBrewer::brewer.pal(3, 'Set3')) +
  scale_x_discrete(breaks=c(1,4), 
                   labels=c('Q1', 'Q4')) +
  theme_bw() +
  facet_wrap(variable~., scales='free_y', nrow=1, strip.position = "left") +
  theme(legend.position='none', strip.background = element_blank(),
        strip.placement = "outside") +
  xlab('Lactate Quartile') + ylab(' ')

#calculate density estimates for boxplots
smydf3 = split(mydf3, as.character(mydf3$lactate_quartile))
mydf3_save = lapply(smydf3, function(df){
  data.frame('value.pufa_mufa_ratio' = density(df$POFA_MOFA_Ratio)$x,
             'density_est.pufa_mufa_ratio' = density(df$POFA_MOFA_Ratio)$y,
             
             'value.mufa_tfa_ratio' = density(df$MOFA_Total_Ratio)$x,
             'density_est.mufa_tfa_ratio' = density(df$MOFA_Total_Ratio)$y,
             
             'value.sfa_tfa_ratio' = density(df$SFA_Total_Ratio)$x,
             'density_est.sfa_tfa_ratio' = density(df$SFA_Total_Ratio)$y)
}) %>%
  melt(measure.vars=c()) %>% 
  rename(lactate_quartile=L1)
write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/fs15a.csv', row.names = F)

## panel b
rm(list=setdiff(ls(), c('mydf', 'mymetadf', 'pa')))
met_fields = c('sex', 'age_at_rec', 'PC', 'TFA', 'POFA_Total_Ratio', 'MOFA_Total_Ratio')
health_fields = c('weight_change')

cols = mymetadf %>% 
  filter(name %in% c(health_fields, met_fields)) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(c(cols$field_id), 'f.eid') %>% drop_na()
colnames(mydf2) = setNames(c(cols$name, 'f.eid'), c(cols$field_id, 'f.eid'))[colnames(mydf2)]
# 
mydf3 =
  mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  select(-PC, -TFA, -age_at_rec, -sex, -f.eid) %>% 
  mutate(weight_change = as.character(weight_change)) %>% 
  filter(weight_change %in% c('No - weigh about the same', 'Yes - gained weight', 'Yes - lost weight')) %>% 
  mutate(weight_change = setNames(c('No change', 'Gained weight', 'Lost weight'), c('No - weigh about the same', 'Yes - gained weight', 'Yes - lost weight'))[weight_change])
pb = mydf3 %>% 
  melt(measure.vars=c('PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio')) %>% 
  mutate(variable = setNames(c('PC/TFA Ratio', 'PUFA/TFA Ratio', 'MUFA/TFA Ratio'),
                             c('PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio'))[variable]) %>% 
  drop_na() %>% 
  mutate(weight_change = factor(weight_change, levels=rev(c('Gained weight', 'No change', 'Lost weight')))) %>%
  ggplot(aes(y=weight_change, x=value, fill=weight_change)) +
  geom_boxplot(alpha=.6, outliers = F) +
  geom_signif(test="wilcox.test", comparisons = combn(c('Gained weight', 'No change', 'Lost weight'), 2, simplify = F)[-4], step_increase = 0.5, map_signif_level = T) +
  scale_x_continuous(expand=expand_scale(.22,0)) +
  facet_wrap(variable~., scales='free_x', strip.position = "bottom", nrow=1) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(y='Weight change', x=' ')

#calculate density estimates for boxplots
smydf3 = split(mydf3, as.character(mydf3$weight_change))
mydf3_save = lapply(smydf3, function(df){
  data.frame('value.pc_tfa_ratio' = density(df$PC_TFA_Ratio)$x,
             'density_est.pc_tfa_ratio' = density(df$PC_TFA_Ratio)$y,
             
             'value.mufa_tfa_ratio' = density(df$MOFA_Total_Ratio)$x,
             'density_est.mufa_tfa_ratio' = density(df$MOFA_Total_Ratio)$y,
             
             'value.pufa_tfa_ratio' = density(df$POFA_Total_Ratio)$x,
             'density_est.pufa_tfa_ratio' = density(df$POFA_Total_Ratio)$y)
}) %>%
  melt(measure.vars=c()) %>% 
  rename(weight_change=L1)
write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/fs15b.csv', row.names = F)

## panel c-d
rm(list=setdiff(ls(), c('mydf', 'mymetadf', 'pa', 'pb')))
illness_codes = read.delim('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/illness_codings.tsv')
diabetes_codes = illness_codes %>% 
  filter(startsWith(meaning, 'E10') | startsWith(meaning, 'E11') | startsWith(meaning, 'E12') | startsWith(meaning, 'E13') | startsWith(meaning, 'E14')) %>% 
  filter(selectable == 'Y') %>% 
  pull(coding)

illness_ids = mymetadf %>% 
  filter(name == 'illnesses') %>% 
  pull(field_id)

diabetes_df = mydf %>% 
  select(f.eid, all_of(illness_ids)) %>% 
  melt(id.var='f.eid') %>% 
  drop_na() %>% 
  filter(value %in% diabetes_codes) %>% 
  select(f.eid) %>% 
  distinct() %>% 
  mutate(diabetes = 'patients')

met_fields = c('sex', 'age_at_rec', 'PC', 'POFA', 'MOFA', 'TFA', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')

cols = mymetadf %>% 
  filter(name %in% c(met_fields)) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(c(cols$field_id), 'f.eid') %>% drop_na()
colnames(mydf2) = setNames(c(cols$name, 'f.eid'), c(cols$field_id, 'f.eid'))[colnames(mydf2)]

mydf3 = mydf2 %>% 
  left_join(diabetes_df, by='f.eid') %>% 
  replace_na(list(diabetes='healthy')) %>% 
  select(-POFA_Total_Ratio, -TFA, -f.eid, -sex)

pc = mydf3 %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = setNames(c('Age', 'PC', 'PUFA', 'MUFA', 'MUFA/TFA Ratio', 'PUFA/MUFA Ratio'),
                             c('age_at_rec', 'PC', 'POFA', 'MOFA', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))[variable]) %>% 
  mutate(diabetes = gsub('patients', 'Diabetic', diabetes),
         diabetes = gsub('healthy', 'Non-diabetic', diabetes),
         diabetes = factor(diabetes, levels=c('Non-diabetic', 'Diabetic')),
         variable = factor(variable, levels=c('Age', 'PC', 'PUFA', 'MUFA', 'MUFA/TFA Ratio', 'PUFA/MUFA Ratio'))) %>% 
  drop_na() %>% 
  ggplot(aes(x=diabetes, y=value, fill=diabetes)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  facet_wrap(variable~., scales='free', nrow=2, strip.position = "left") +
  geom_signif(test="wilcox.test", comparisons = list(c("Diabetic", "Non-diabetic")), map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.15,0)) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.placement = "outside") +
  xlab("") + ylab(' ')


# calculate density estimates for boxplots
smydf3 = split(mydf3, as.character(mydf3$diabetes))
mydf3_save = lapply(smydf3, function(df){
  data.frame('value.age' = density(df$age_at_rec)$x,
             'density_est.age' = density(df$age_at_rec)$y,
             
             'value.pc' = density(df$PC)$x,
             'density_est.pc' = density(df$PC)$y,
             
             'value.pufa' = density(df$POFA)$x,
             'density_est.pufa' = density(df$POFA)$y,
             
             'value.mufa' = density(df$MOFA)$x,
             'density_est.mufa' = density(df$MOFA)$y,
             
             
             'value.mufa_total_ratio' = density(df$MOFA_Total_Ratio)$x,
             'density_est.mufa_total_ratio' = density(df$MOFA_Total_Ratio)$y,
             
             'value.pufa_mufa_ratio' = density(df$POFA_MOFA_Ratio)$x,
             'density_est.pufa_mufa_ratio' = density(df$POFA_MOFA_Ratio)$y
             )
}) %>%
  melt(measure.vars=c()) %>% 
  rename(diabetes=L1)
write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/fs15cd.csv', row.names = F)


#combine plots
p = ggarrange(pa,pb,pc,ncol=1, heights=c(1,1.4,1.8))
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/figureS15.pdf', p, width=8, height=8)

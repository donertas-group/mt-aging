# main figure 7 panels c-h and associated tables
library(tidyverse)
library(reshape2)
library(ggpubr)

# panel c
mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

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
  theme(legend.position = 'right')

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
  theme(legend.position = 'right')

pc = ggarrange(pc1, pc2, ncol=1, heights = c(1, 1), common.legend=T, legend = 'right')

# get loess preds
smydf3 = split(mydf3, mydf3$Sex)
mydf3_save = lapply(smydf3, function(df){
  df$pc.loess_preds = stats::predict(stats::loess(PC ~ age, data=df))
  df$pc_tfa_ratio.loess_preds = stats::predict(stats::loess(PC_TFA_Ratio ~ age, data=df))
  df
}) %>% 
  bind_rows() %>% 
  select(-PC, -PC_TFA_Ratio)

write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/f7c.csv', row.names = F)

# panel d
rm(list=setdiff(ls(), c('mydf', 'mymetadf', 'pc')))
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
  
pd1 = mydf3 %>% 
  ggplot(aes(x=lactate_quartile, y=PC_TFA_Ratio, fill=lactate_quartile)) +
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
  theme(legend.position='none') +
  xlab('Lactate Quartile') + ylab('PC/TFA Ratio')

pd2=mydf3 %>% 
  ggplot(aes(x=lactate_quartile, y=PUFA_Total_Ratio, fill=lactate_quartile)) +
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
  theme(legend.position='none') +
  xlab('Lactate Quartile') + ylab('PUFA/TFA Ratio')
pd = ggarrange(pd1, pd2, ncol=2)

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
write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/f7d-e.csv', row.names = F)

## panels f and g
# get illnesses
rm(list=setdiff(ls(), c('mydf', 'mymetadf', 'pd', 'pc')))
illness_codes = read.delim('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/illness_codings.tsv')
diabetes_codes = illness_codes %>% 
  filter(startsWith(meaning, 'E10') | startsWith(meaning, 'E11') | startsWith(meaning, 'E12') | startsWith(meaning, 'E13') | startsWith(meaning, 'E14')) %>% 
  filter(selectable == 'Y') %>% 
  pull(coding)

illness_ids = mymetadf %>% 
  filter(name == 'illnesses') %>% 
  pull(field_id)

diabetes_df =
  mydf %>% 
  select(f.eid, all_of(illness_ids)) %>% 
  melt(id.var='f.eid') %>% 
  drop_na() %>% 
  filter(value %in% diabetes_codes) %>% 
  select(f.eid) %>% 
  distinct() %>% 
  mutate(diabetes = 'patients')

# get other conditions and join diabetes info
met_fields = c('sex', 'age_at_rec', 'PC', 'TFA', 'POFA_Total_Ratio')

cols = mymetadf %>% 
  filter(name %in% c(met_fields)) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(c(cols$field_id), 'f.eid') 
colnames(mydf2) = setNames(c(cols$name, 'f.eid'), c(cols$field_id, 'f.eid'))[colnames(mydf2)]

mydf3 = mydf2 %>% 
  left_join(diabetes_df, by='f.eid') %>% 
  replace_na(list(diabetes='healthy')) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  select(-f.eid, -age_at_rec) %>% 
  transmute(diabetes = gsub('patients', 'Diabetic', diabetes),
            diabetes = gsub('healthy', 'Non-diabetic', diabetes),
            diabetes = factor(diabetes, levels=c('Non-diabetic', 'Diabetic')),
            PC_TFA_Ratio, PUFA_TFA_Ratio=POFA_Total_Ratio) %>% 
  drop_na() 

pf1=mydf3 %>% 
  ggplot(aes(x=diabetes, y=PC_TFA_Ratio, fill=diabetes)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  # facet_wrap(variable~., scales='free_y', nrow=2) +
  geom_signif(test="wilcox.test", comparisons = list(c("Diabetic", "Non-diabetic")), map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.15,0)) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  theme(legend.position = 'none') +
  xlab("Diabetes status") + ylab('PC/TFA Ratio')


pf2=mydf3 %>% 
  ggplot(aes(x=diabetes, y=PUFA_TFA_Ratio, fill=diabetes)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  # facet_wrap(variable~., scales='free_y', nrow=2) +
  geom_signif(test="wilcox.test", comparisons = list(c("Diabetic", "Non-diabetic")), map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.15,0)) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  theme(legend.position = 'none') +
  xlab("Diabetes status") + ylab('PUFA/TFA Ratio')

pf = ggarrange(pf1, pf2, ncol=2)
mydf3 %>% head
# calculate density estimates for boxplots
smydf3 = split(mydf3, as.character(mydf3$diabetes))
mydf3_save = lapply(smydf3, function(df){
  data.frame('value.pc_tfa_ratio' = density(df$PC_TFA_Ratio)$x,
             'density_est.pc_tfa_ratio' = density(df$PC_TFA_Ratio)$y,
             
             'value.pufa_tfa_ratio' = density(df$PUFA_TFA_Ratio)$x,
             'density_est.pufa_tfa_ratio' = density(df$PUFA_TFA_Ratio)$y)
}) %>%
  melt(measure.vars=c()) %>% 
  rename(diabetes=L1)

write.csv(mydf3_save, file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/f7f-g.csv', row.names = F)

# panel h 
rm(list=setdiff(ls(), c('mydf', 'mymetadf', 'pc', 'pd', 'pf')))
met_fields = c('age_at_rec', 'PC', 'POFA', 'MOFA', 'TFA', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')
health_fields = c('basal_metabolic_rate', 'cci', 'max_digits_remembered', 'walking_pace')

# change colnames
cols = mymetadf %>% 
  filter(name %in% c(health_fields, met_fields)) %>% 
  select(field_id, name) %>% 
  filter(!field_id %in% c('f.399.0.2', 'f.399.0.3'))

mydf2 = mydf %>% select(c(cols$field_id, 'cci')) 
colnames(mydf2) = setNames(c(cols$name, 'cci'), c(cols$field_id, 'cci'))[colnames(mydf2)]

mydf4 = mydf2 %>% 
  mutate(walking_pace = setNames(c(NA, NA, 1, 2, 3), 
                                 c('None of the above', 'Prefer not to answer', 'Slow pace', 'Steady average pace', 'Brisk pace'))[walking_pace]) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  select(-TFA) %>% 
  mutate_all(as.numeric)

combdf = expand.grid(c(health_fields, "cci"), setdiff(colnames(mydf4), c(health_fields, "cci", "age_at_rec", "age_int")))

corrsdf =
  apply(combdf, 1, function(x){
    v1 = x[1] %>% as.character()
    v2 = x[2] %>% as.character()
    cres = cor.test(mydf4[,v1], mydf4[,v2], m='s')
    data.frame(param=v1, met=v2, 'rho'=cres$estimate, 'p.value'=cres$p.value)
  }) %>% 
  melt(measure.vars=c()) %>% 
  select(-L1) %>% 
  mutate(rho = as.numeric(rho),
         p.value = as.numeric(p.value))

ph = corrsdf %>% 
  mutate(lab = as.character(round(rho,2))) %>% 
  mutate(sig = ifelse(p.value < 0.1, '*', ' '),
         sig = ifelse(p.value < 0.05, '**', ' '),
         sig = ifelse(p.value < 0.01, '***', ' ')) %>%
  mutate(met=gsub('POFA', 'PUFA', met),
         met=gsub('MOFA', 'MUFA', met),
         met=gsub('Total', 'TFA', met)) %>% 
  mutate(met = factor(met, levels=c('PUFA', 'MUFA', 'PC', 'PC_TFA_Ratio', 'PUFA_TFA_Ratio', 'MUFA_TFA_Ratio', 'PUFA_MUFA_Ratio'))) %>% 
  mutate(param= factor(param, levels=rev(c('basal_metabolic_rate', 'cci', 'max_digits_remembered', 'walking_pace')))) %>%
  mutate(lbl = paste0(lab, sig)) %>% 
  ggplot(aes(x=met, y=param, fill=rho)) +
  geom_tile() +
  geom_text(aes(label=lbl)) +
  labs(fill = 'Correlation\ncoeff.') +
  scale_fill_gradient2(high = 'firebrick2', mid = 'white', low = 'dodgerblue4', 
                       limits = c(-.42, .42), midpoint = 0)  +
  scale_x_discrete(name='Metabolites',
                   breaks=c('PUFA', 'MUFA', 'PC', 'PC_TFA_Ratio', 'PUFA_TFA_Ratio', 'MUFA_TFA_Ratio', 'PUFA_MUFA_Ratio'),
                   labels=c('PUFA', 'MUFA', 'PC', 'PC/TFA Ratio', 'PUFA/TFA Ratio', 'MUFA/TFA Ratio', 'PUFA/MUFA Ratio')) +
  scale_y_discrete(name='Health parameters',
                   breaks=c('basal_metabolic_rate', 'cci', 'max_digits_remembered', 'walking_pace'),
                   labels=c('Basal Metabolic Rate', 'CCI', 'Max Digits Remembered', 'Walking Pace')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=60, hjust = 1)) +
  xlab('Metabolites') + ylab('Health parameters')

write.csv(select(ph$data, -c(lab, sig)), file = '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/tables/f7h.csv', row.names = F)
# final figure
p = ggarrange(ggarrange(pc, ggarrange(pd, pf, ncol=1), widths = c(.8, 1)), ph, nrow=2, heights = c(1.2, 1))
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/figure7c-h.pdf', p, width=8, height=8)

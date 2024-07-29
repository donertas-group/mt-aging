# association between health parameters and metabolites
library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggridges)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

# get illnesses
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

# get other conditions and join diabetes info
met_fields = c('sex', 'age_at_rec', 'PC', 'POFA', 'MOFA', 'TFA', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')
health_fields = c('falls_in_last_year', 'fractures_in_5_years', 'fractures_from_falls', 'weight_change')

cols = mymetadf %>% 
  filter(name %in% c(health_fields, met_fields)) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(c(cols$field_id), 'f.eid') 
colnames(mydf2) = setNames(c(cols$name, 'f.eid'), c(cols$field_id, 'f.eid'))[colnames(mydf2)]

mydf3 = mydf2 %>% 
  left_join(diabetes_df, by='f.eid') %>% 
  replace_na(list(diabetes='healthy')) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  select(-f.eid) %>% 
  mutate(age_int = cut(age_at_rec, breaks=c(35,40,45,50,55,60,65,70,75)),
         age_int = as.character(age_int)) 

## TODO: calculate correlation on age-stratified inds



# stratify by medical condition
# diabetes
mydf3 %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         diabetes = gsub('patients', 'Diabetic', diabetes),
         diabetes = gsub('healthy', 'Non-diabetic', diabetes),
         diabetes = factor(diabetes, levels=c('Non-diabetic', 'Diabetic')),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  ggplot(aes(x=age_int, y=value, fill=diabetes)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  facet_wrap(variable~., scales='free_y', nrow=2) +
  geom_signif(test="wilcox.test", comparisons = list(c("Diabetic", "Non-diabetic")), map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.15,0)) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # theme(legend.position = 'none') +
  xlab("") + ylab('Value')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_diabetes.png', width=13, height=4)

p1=mydf3 %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         diabetes = gsub('patients', 'Diabetic', diabetes),
         diabetes = gsub('healthy', 'Non-diabetic', diabetes),
         diabetes = factor(diabetes, levels=c('Non-diabetic', 'Diabetic')),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  ggplot(aes(x=value, y='', fill=diabetes)) +
  geom_density_ridges(alpha=.4, scale=10, panel_scaling = F) +
  facet_wrap(~variable, scales='free', nrow=1) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  # theme(legend.position = 'none') +
  ylab("") + xlab('Value') +
  ggtitle('a.Diabetes')

# weight change
mydf3 %>% 
  mutate(weight_change = as.character(weight_change)) %>% 
  filter(weight_change %in% c('No - weigh about the same', 'Yes - gained weight', 'Yes - lost weight')) %>% 
  mutate(weight_change = setNames(c('No change', 'Gained weight', 'Lost weight'), c('No - weigh about the same', 'Yes - gained weight', 'Yes - lost weight'))[weight_change]) %>%
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  mutate(weight_change = factor(weight_change, levels=rev(c('Gained weight', 'No change', 'Lost weight')))) %>%
  ggplot(aes(x=age_int, y=value, fill=weight_change)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  geom_signif(test="wilcox.test", comparisons = combn(c('Gained weight', 'No change', 'Lost weight'), 2, simplify = F)[-4], step_increase = 0.5, map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.22,0)) +
  facet_wrap(variable~., scales='free_x', nrow=2) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # theme(legend.position = 'none') +
  xlab("Weight change") + ylab('Value') +
  coord_flip()
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_weight_change.png', width=14, height=4)

p2 = mydf3 %>% 
  mutate(weight_change = as.character(weight_change)) %>% 
  filter(weight_change %in% c('No - weigh about the same', 'Yes - gained weight', 'Yes - lost weight')) %>% 
  mutate(weight_change = setNames(c('No change', 'Gained weight', 'Lost weight'), c('No - weigh about the same', 'Yes - gained weight', 'Yes - lost weight'))[weight_change]) %>%
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  mutate(weight_change = factor(weight_change, levels=rev(c('Gained weight', 'No change', 'Lost weight')))) %>% 
  ggplot(aes(x=value, y='', fill=weight_change)) +
  geom_density_ridges(alpha=.4, scale=10, panel_scaling = F) +
  facet_wrap(~variable, scales='free', nrow=1) +
  theme_bw() +
  scale_fill_manual(name='Weight change', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  # theme(legend.position = 'none') +
  ylab("") + xlab('Value') +
  ggtitle('b.Weight change')

# falls
mydf3 %>% 
  filter(falls_in_last_year %in% c('No falls', 'Only one fall', 'More than one fall')) %>% 
  mutate(falls_in_last_year = gsub('Only one fall', 'One fall', falls_in_last_year),
         falls_in_last_year = gsub('More than one fall', '1+ falls', falls_in_last_year)) %>%
  mutate(falls_in_last_year = factor(falls_in_last_year, levels=c('No falls', 'One fall', '1+ falls'))) %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  ggplot(aes(x=age_int, y=value, fill=falls_in_last_year)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  geom_signif(test="wilcox.test", comparisons = combn(c('No falls', 'One fall', '1+ falls'), 2, simplify = F)[-4], step_increase = 0.5, map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.22,0)) +
  facet_wrap(variable~., scales='free_x', nrow=2) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # theme(legend.position = 'none') +
  xlab("Falls in last year") + ylab('Value') +
  coord_flip()
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_falls_in_last_year.png', width=12, height=4)

p3 = mydf3 %>% 
  filter(falls_in_last_year %in% c('No falls', 'Only one fall', 'More than one fall')) %>% 
  mutate(falls_in_last_year = gsub('Only one fall', 'One fall', falls_in_last_year),
         falls_in_last_year = gsub('More than one fall', '1+ falls', falls_in_last_year)) %>%
  mutate(falls_in_last_year = factor(falls_in_last_year, levels=c('No falls', 'One fall', '1+ falls'))) %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  # mutate(weight_change = factor(weight_change, levels=rev(c('Gained weight', 'No change', 'Lost weight')))) %>% 
  ggplot(aes(x=value, y='', fill=falls_in_last_year)) +
  geom_density_ridges(alpha=.4, scale=10, panel_scaling = F) +
  facet_wrap(~variable, scales='free', nrow=1) +
  theme_bw() +
  scale_fill_manual(name='Falls in last year', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  # theme(legend.position = 'none') +
  ylab("") + xlab('Value') +
  ggtitle('c.Falls in last year')

# fractures from falls
mydf3 %>% 
  filter(fractures_from_falls %in% c('Yes', 'No')) %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  ggplot(aes(x=age_int, y=value, fill=fractures_from_falls)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  geom_signif(test="wilcox.test", comparisons = list(c('No', 'Yes')), map_signif_level = T) +
  scale_y_continuous(expand=expand_scale(.22,0)) +
  facet_wrap(variable~., scales='free_y', nrow=2) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # theme(legend.position = 'none') +
  xlab("Fractures from falls") + ylab('Value')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_fractures_from_falls.png', width=12, height=4)


p4 = mydf3 %>% 
  filter(fractures_from_falls %in% c('Yes', 'No')) %>% 
  melt(measure.vars=c('age_at_rec', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable),
         variable = factor(variable, levels=c('Age', 'PC', 'POFA', 'MOFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>%
  # mutate(weight_change = factor(weight_change, levels=rev(c('Gained weight', 'No change', 'Lost weight')))) %>% 
  ggplot(aes(x=value, y='', fill=fractures_from_falls)) +
  geom_density_ridges(alpha=.4, scale=10, panel_scaling = F) +
  facet_wrap(~variable, scales='free', nrow=1) +
  theme_bw() +
  scale_fill_manual(name='Fractures\nfrom falls', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  # scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  # theme(legend.position = 'none') +
  ylab("") + xlab('Value') +
  ggtitle('d.Fractures from falls')


ggarrange(p1, p2, p3, p4, ncol=1)
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_densities.png', width=14, height=7)

# association between health parameters and metabolites
library(tidyverse)
library(reshape2)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

met_fields = c('sex', 'age_at_rec', 'PC', 'POFA', 'MOFA', 'TFA', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')
health_fields = c('basal_metabolic_rate', 'n_incorr_matches',  
                 'max_digits_remembered', 'fluid_int_score', 'grip_strength_left', 'grip_strength_right', 'walking_pace')

# three measurements for n_incorr_matches and n_incorr_matches_pilot
# get the average of all three measurements
df_399 = mydf %>% 
  select(f.eid, `f.399.0.1`, `f.399.0.2`, `f.399.0.3`) %>% 
  melt(id.var='f.eid') %>% 
  drop_na() %>% 
  group_by(f.eid) %>% 
  summarise(f.399.0.1 = mean(value))
mydf1 = mydf %>% 
  select(-c(`f.399.0.1`, `f.399.0.2`, `f.399.0.3`)) %>% 
  left_join(df_399, by='f.eid')

# change colnames
cols = mymetadf %>% 
  filter(name %in% c(health_fields, met_fields)) %>% 
  select(field_id, name) %>% 
  filter(!field_id %in% c('f.399.0.2', 'f.399.0.3'))

mydf2 = mydf1 %>% select(c(cols$field_id, 'cci')) 
colnames(mydf2) = setNames(c(cols$name, 'cci'), c(cols$field_id, 'cci'))[colnames(mydf2)]

# stratify by health outcome quartiles
# #check dist of quantiles
# mydf2 %>% 
#   select(-walking_pace) %>% 
#   melt(id.vars=c(met_fields)) %>% 
#   group_by(variable) %>% 
#   mutate(quantile = ntile(value,4)) %>% 
#   mutate(value = as.numeric(value)) %>% 
#   mutate(quantile = as.factor(quantile)) %>% 
#   ggplot(aes(x=variable, y=value, fill=quantile)) +
#   geom_boxplot(outlier.size = .4) +
#   facet_wrap(~variable, scales='free')

for (v in setdiff(c(health_fields), 'walking_pace')){
  mydf3 =
    mydf2 %>% 
    mutate(PC_TFA_Ratio = PC / TFA) %>% 
    mutate(walking_pace = setNames(c(NA, NA, 1, 2, 3), 
                                   c('None of the above', 'Prefer not to answer', 'Slow pace', 'Steady average pace', 'Brisk pace'))[walking_pace]) %>% 
    melt(id.vars=c(met_fields, 'PC_TFA_Ratio', 'cci')) %>% 
    group_by(variable) %>%
    mutate(quantile = ntile(value,4)) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(quantile = as.factor(quantile)) %>% 
    filter(variable == v) %>% 
    pivot_wider(names_from = variable, values_from = value) %>% 
    melt(measure.vars=c(setdiff(met_fields, 'sex'), 'PC_TFA_Ratio', 'cci')) %>% 
    filter(variable != 'TFA') %>% 
    filter(quantile %in% c(1,4)) %>% 
    mutate(variable = gsub('age_at_rec', 'Age', variable)) %>% 
    mutate(variable = factor(variable, levels=c('Age', 'POFA', 'MOFA', 'PC', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
    drop_na() 
  
  mydf3 %>% 
    ggplot(aes(x=quantile, y=value, fill=quantile)) +
    ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, normalize='panels', point_colour = NA) + 
    geom_boxplot(width = .1, outlier.shape = NA) +
    gghalves::geom_half_point(aes(color=quantile), side = "l", range_scale = .4, alpha = .1, size=.1) +
    facet_wrap(variable~., scales='free_y', nrow=2) +
    theme_bw() +
    scale_x_discrete(labels = c('Q1', 'Q4')) +
    xlab('Quantiles') + ylab('Value') +
    ggtitle(paste0('Distribution of metabolites for individuals grouped by ', v)) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
    scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
    theme(legend.position = 'none')
  ggsave(paste0('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_', v, '.png'), width=12, height=6)
}
# walking pace:
mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  mutate(walking_pace = setNames(c(NA, NA, 1, 2, 3), 
                                 c('None of the above', 'Prefer not to answer', 'Slow pace', 'Steady average pace', 'Brisk pace'))[walking_pace]) %>% 
  melt(measure.vars=setdiff(c(met_fields, 'PC_TFA_Ratio', 'cci'), c('sex'))) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable)) %>% 
  mutate(variable = factor(variable, levels=c('Age', 'POFA', 'MOFA', 'PC', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  mutate(walking_pace = as.factor(walking_pace)) %>% 
  ggplot(aes(x=walking_pace, y=value, fill=walking_pace)) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, normalize='panels', point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  gghalves::geom_half_point(aes(color=walking_pace), side = "l", range_scale = .4, alpha = .1, size=.1) +
  facet_wrap(variable~., scales='free_y', nrow=2) +
  theme_bw() +
  scale_x_discrete(breaks=1:3,
                   labels = c('Slow pace', 'Steady pace', 'Brisk pace')) +
  xlab('Usual walking pace') + ylab('Value') +
  ggtitle(paste0('Distribution of metabolites for individuals grouped by walking pace')) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  theme(legend.position = 'none',
        axis.text.x.bottom = element_text(angle=90, vjust = 0.5, hjust=1))
ggsave(paste0('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_walking_pace.png'), width=12, height=6)

# CCI
mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  mutate(quantile = ntile(cci,4)) %>%
  # mutate(value = as.numeric(value)) %>%
  mutate(quantile = as.factor(quantile)) %>% 
  melt(measure.vars=c(setdiff(met_fields, 'sex'), 'PC_TFA_Ratio')) %>% 
  filter(variable != 'TFA') %>% 
  filter(quantile %in% c(1,4)) %>% 
  mutate(variable = gsub('age_at_rec', 'Age', variable)) %>% 
  mutate(variable = factor(variable, levels=c('Age', 'POFA', 'MOFA', 'PC', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  drop_na() %>% 
  ggplot(aes(x=quantile, y=value, fill=quantile)) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, normalize='panels', point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  gghalves::geom_half_point(aes(color=quantile), side = "l", range_scale = .4, alpha = .1, size=.1) +
  facet_wrap(variable~., scales='free_y', nrow=2) +
  theme_bw() +
  scale_x_discrete(labels = c('Q1', 'Q4')) +
  xlab('Quantiles') + ylab('Value') +
  ggtitle(paste0('Distribution of metabolites for individuals grouped by CCI')) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  theme(legend.position = 'none')
ggsave(paste0('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/met_dists_for_cci.png'), width=12, height=6)


## Correlation test between health parameters and metabolites
mydf4 = mydf2 %>% 
  mutate(walking_pace = setNames(c(NA, NA, 1, 2, 3), 
                                 c('None of the above', 'Prefer not to answer', 'Slow pace', 'Steady average pace', 'Brisk pace'))[walking_pace]) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  select(-sex, -TFA) %>% 
  mutate_all(as.numeric)

combdf = expand.grid(c(health_fields, "cci"), setdiff(colnames(mydf4), c(health_fields, "cci", "age_at_rec")))
combdf

corrsdf = apply(combdf, 1, function(x){
  v1 = x[1] %>% as.character()
  v2 = x[2] %>% as.character()
  cres = cor.test(mydf4[,v1], mydf4[,v2], m='s')
  
  c(v1, v2, cres$estimate, cres$p.value)
}) %>% 
  t() %>% 
  as.data.frame() %>% 
  magrittr::set_colnames(c('param', 'met', 'rho', 'p.value')) %>% 
  mutate(rho = as.numeric(rho),
         p.value = as.numeric(p.value))


corrsdf %>% 
  mutate(lab = as.character(round(rho,2))) %>% 
  mutate(sig = ifelse(p.value < 0.1, '*', ' '),
         sig = ifelse(p.value < 0.05, '**', ' '),
         sig = ifelse(p.value < 0.01, '***', ' ')) %>%
  mutate(met = factor(met, levels=c('POFA', 'MOFA', 'PC', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio'))) %>% 
  mutate(param = gsub('cci', 'CCI', param)) %>% 
  mutate(param= factor(param, levels=rev(c('basal_metabolic_rate', 'walking_pace', 'grip_strength_left', 'grip_strength_right', 'CCI', 'n_incorr_matches', 'max_digits_remembered', 'fluid_int_score')))) %>% 
  mutate(label = paste0(lab, sig)) %>% 
  ggplot(aes(x=met, y=param, fill=rho)) +
  geom_tile() +
  geom_text(aes(label=label)) +
  labs(fill = 'Correlation\ncoeff.') +
  scale_fill_gradient2(high = 'firebrick2', mid = 'white', low = 'dodgerblue4', 
                       limits = c(-.42, .42), midpoint = 0)  +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) +
  xlab('Metabolites') + ylab('Health parameters')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/health_params_nmrs_correlation_heatmaps.png', width=8, height=7)



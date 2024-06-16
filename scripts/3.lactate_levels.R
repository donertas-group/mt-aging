# association of lactate with age and other variables
# lactate -> as a proxy for mt-health
library(tidyverse)
library(reshape2)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

fields = c('sex', 'age_at_rec', 'PC', 'POFA', 'MOFA', 'SFA', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio', 'TFA', 'Lactate')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

# Correlation between lactate and other metabolites
bind_rows(
  mydf2 %>%
    mutate(PC_TFA_Ratio = PC / TFA) %>% 
    rename(age = age_at_rec) %>% 
    filter(sex == 'Female') %>% 
    select(-sex) %>% 
    mutate_all(as.numeric) %>% 
    drop_na() %>% 
    cor(method = 'spearman') %>% 
    as.data.frame() %>% 
    rownames_to_column('var1') %>% 
    melt(id.var='var1') %>% 
    filter(var1 == 'Lactate') %>% 
    filter(variable != 'Lactate') %>% 
    mutate(gr = 'Female'),
  mydf2 %>%
    mutate(PC_TFA_Ratio = PC / TFA) %>% 
    rename(age = age_at_rec) %>% 
    filter(sex == 'Male') %>% 
    select(-sex) %>% 
    mutate_all(as.numeric) %>% 
    drop_na() %>% 
    cor(method = 'spearman') %>% 
    as.data.frame() %>% 
    rownames_to_column('var1') %>% 
    melt(id.var='var1') %>% 
    filter(var1 == 'Lactate') %>% 
    filter(variable != 'Lactate') %>% 
    mutate(gr = 'Male'),
  mydf2 %>%
    mutate(PC_TFA_Ratio = PC / TFA) %>% 
    rename(age = age_at_rec) %>% 
    # filter(sex == 'Male') %>% 
    select(-sex) %>% 
    mutate_all(as.numeric) %>% 
    drop_na() %>% 
    cor(method = 'spearman') %>% 
    as.data.frame() %>% 
    rownames_to_column('var1') %>% 
    melt(id.var='var1') %>% 
    filter(var1 == 'Lactate') %>% 
    filter(variable != 'Lactate') %>% 
    mutate(gr = 'Combined')
) %>% 
  ggplot(aes(x=value, y=reorder(variable, value))) +
  geom_bar(stat='identity', alpha=.8) +
  geom_label(aes(label=round(value,3)), size=3) +
  facet_grid(~gr) +
  theme_bw() +
  xlim(c(-0.18, 0.18)) +
  ylab('Variables') +
  xlab("Spearman's correlation coefficient")
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/lactate_corrs.png', width=10, height=5)

# #check p-values
# cordf = mydf2 %>%
#   mutate(PC_TFA_Ratio = PC / TFA) %>% 
#   rename(age = age_at_rec) %>% 
#   select(-sex) %>% 
#   mutate_all(as.numeric) %>% 
#   drop_na()
# 
# apply(cordf, 2, function(x){
#   corres = cor.test(x, cordf[,'Lactate'], m='s')
#   corres$p.value
# })

# Lactate vs metabolites scatter-plot
# pp = 
  mydf2 %>%
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  rename(age = age_at_rec) %>% 
  drop_na() %>% 
  melt(id.vars=c('sex', 'Lactate')) %>% 
  ggplot(aes(x=Lactate, y=value, color=sex)) +
  geom_point(size=.2, alpha=.1) +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  geom_smooth(se=FALSE) +
  # geom_quantile(quantiles=c(.25, .75), formula=y ~ x, linetype=2) +
  # geom_quantile(quantiles=c(5), formula=y ~ x) +
  facet_wrap(variable~., scales = 'free_y', ncol=4) +
  theme_bw()
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/lactate_vs_nmrs.png', width=10, height=5)

# per quantile for sexes seperately
mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  rename(age=age_at_rec) %>% 
  mutate(lactate_quartile = ntile(Lactate, 4)) %>% 
  filter(lactate_quartile %in% c(1,4)) %>% 
  mutate(lactate_quartile = factor(lactate_quartile, levels=1:4)) %>% 
  select(-Lactate) %>% 
  melt(id.vars=c('sex', 'lactate_quartile')) %>% 
  drop_na() %>% 
  mutate(variable = factor(variable, levels=c('age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  ggplot(aes(x=lactate_quartile, y=value, fill=sex)) +
  geom_boxplot(alpha=.6, outlier.size=.3) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  facet_wrap(variable~., scales='free_y') +
  theme_bw() +
  xlab('Lactate Quartile') + ylab('Value')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/lactate_quartiles_boxplots.png', width=10, height=5)


# per quantile for sexes seperately
mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  rename(age=age_at_rec) %>% 
  mutate(lactate_quartile = ntile(Lactate, 4)) %>% 
  filter(lactate_quartile %in% c(1,4)) %>% 
  mutate(lactate_quartile = factor(lactate_quartile, levels=1:4)) %>% 
  select(-Lactate, -sex) %>% 
  melt(id.vars=c('lactate_quartile')) %>% 
  drop_na() %>% 
  mutate(variable = factor(variable, levels=c('age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
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
  theme(legend.position='none') +
  facet_wrap(variable~., scales='free_y', nrow=2) +
  xlab('Lactate Quartile') + ylab('Value')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/lactate_quartiles_boxplots2.png', width=10, height=5)


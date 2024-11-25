# Plot age_vs_nmr metabolomics' correlation
library(tidyverse)
library(reshape2)
library(ggpubr)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

fields = c('sex', 'age_at_rec', 'PC', 'POFA', 'MOFA', 'SFA', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio', 'TFA')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

# #trajectories
# mydf2 %>% 
#   mutate(PC_TFA_Ratio = PC / TFA) %>% 
#   melt(id.vars=c('age_at_rec', 'sex')) %>% 
#   mutate(age = as.numeric(age_at_rec)) %>%
#   drop_na() %>% 
#   mutate(variable = factor(variable, levels=c('PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
#   ggplot(aes(x=age, y=value, color=sex)) +
#   geom_smooth() +
#   facet_wrap(~variable, scales='free_y', ncol=5) + 
#   ylab('Value') + xlab('Age') +
#   ggtitle('Age-associated changes') +
#   theme_bw() +
#   theme(legend.position = 'bottom')
# ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/age_vs_nmrs.png', width=10, height=4)

#trajectories - with quantiles
mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  melt(id.vars=c('age_at_rec', 'sex')) %>% 
  mutate(age = as.numeric(age_at_rec)) %>%
  drop_na() %>% 
  mutate(variable = factor(variable, levels=c('PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  ggplot(aes(x=age, y=value, color=sex)) +
  geom_smooth(se=F) +
  geom_quantile(quantiles=c(.25, .75), formula=y ~ poly(x, 3), linetype=2) +
  # geom_quantile(quantiles=c(.5), formula=y ~ poly(x, 3), colour="red") +
  facet_wrap(~variable, scales='free_y', ncol=5) + 
  ylab('Value') + xlab('Age') +
  scale_color_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  ggtitle('Age-associated changes') +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/age_vs_nmrs.pdf', width=10, height=4)

# boxplots
mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  melt(id.vars=c('age_at_rec', 'sex')) %>% 
  mutate(age = factor(age_at_rec, levels=37:73)) %>% 
  mutate(variable = factor(variable, levels=c('PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  ggplot(aes(x=age, y=value, fill=sex)) +
  geom_boxplot(outlier.size = 0.001, size=.2, alpha=.9, outlier.shape=NA) +
  facet_wrap(~variable, scales = 'free_y', ncol=5) +
  ylab('Value') + xlab('Age') +
  scale_fill_manual(values=RColorBrewer::brewer.pal(3, 'Set2')) +
  scale_x_discrete(breaks=c(37, 40, 45, 50, 55, 60, 65, 70, 73)) +
  theme_bw() +
  theme(legend.position = 'none')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/age_vs_nmrs_boxplots.png', width=11, height=5)


# correlations between metabolites
ph1 =
  mydf2 %>% 
  filter(sex == 'Female') %>% 
  select(-sex) %>% 
  rename(Age = age_at_rec) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  mutate_all(as.numeric) %>%
  drop_na() %>% 
  # group_by(sex) %>%
  cor(method = 'spearman') %>%
  as.data.frame() %>% 
  rownames_to_column('v1') %>% 
  melt(id.vars = 'v1') %>% 
  mutate(variable = factor(variable, levels=c('Age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  mutate(v1 = factor(v1, levels=c('Age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  ggplot(aes(x=v1, y=variable, fill=value)) +
  geom_tile(color = "white") + 
  geom_text(aes(label=round(value, 2))) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       name = 'Spearmans\ncorr coeff.') +
  coord_fixed() +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle('Correlation heatmap', subtitle = 'Females') +
  xlab('') + ylab('')

ph2 = mydf2 %>% 
  filter(sex == 'Male') %>% 
  select(-sex) %>% 
  mutate_all(as.numeric) %>%
  drop_na() %>% 
  rename(Age = age_at_rec) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  # group_by(sex) %>%
  cor(method = 'spearman') %>%
  as.data.frame() %>% 
  rownames_to_column('v1') %>%
  melt(id.vars = 'v1') %>% 
  mutate(variable = factor(variable, levels=c('Age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  mutate(v1 = factor(v1, levels=c('Age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  ggplot(aes(x=v1, y=variable, fill=value)) +
  geom_tile(color = "white") + 
  geom_text(aes(label=round(value, 2))) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       name = 'Spearmans\ncorr coeff.') +
  coord_fixed() +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle('Correlation heatmap', subtitle = 'Males') +
  xlab('') + ylab('')
ggarrange(ph1, ph2, ncol=2, common.legend = T)
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/nmrs_correlation_heatmaps_persex.png', width=14, height=7)


mydf2 %>% 
  select(-sex) %>% 
  rename(Age = age_at_rec) %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>%
  mutate_all(as.numeric) %>%
  drop_na() %>% 
  # group_by(sex) %>%
  cor(method = 'spearman') %>%
  as.data.frame() %>% 
  rownames_to_column('v1') %>%
  melt(id.vars = 'v1') %>% 
  mutate(variable = factor(variable, levels=c('Age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  mutate(v1 = factor(v1, levels=c('Age', 'PC', 'TFA', 'POFA', 'MOFA', 'SFA', 'PC_TFA_Ratio', 'POFA_MOFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'SFA_Total_Ratio'))) %>% 
  ggplot(aes(x=v1, y=variable, fill=value)) +
  geom_tile(color = "white") + 
  geom_text(aes(label=round(value, 2))) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       name = 'Spearmans\ncorr coeff.') +
  coord_fixed() +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle('Correlation heatmap', subtitle = 'Females & Males Combined') +
  xlab('') + ylab('')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/nmrs_correlation_heatmaps_combined.png', width=8, height=7)

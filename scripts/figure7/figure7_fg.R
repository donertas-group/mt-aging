# main figure 7 panels f-g
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

margins = margin(5.5,5.5,5.5,5.5, "pt")   # t r b l 
significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

source('scripts/files.R')

illness_codes = read.delim(file.path(project_dir, 'data', 'illness_codings.tsv'))
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
  geom_signif(test="wilcox.test", comparisons = list(c("Diabetic", "Non-diabetic")), map_signif_level = significances) +
  scale_y_continuous(expand=expand_scale(.15,0)) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  theme(legend.position = 'none',
        plot.margin = margins) +
  xlab("Diabetes status") + ylab('PC/TFA Ratio')


pf2=mydf3 %>% 
  ggplot(aes(x=diabetes, y=PUFA_TFA_Ratio, fill=diabetes)) +
  geom_boxplot(outlier.shape = NA, alpha=.6) +
  # facet_wrap(variable~., scales='free_y', nrow=2) +
  geom_signif(test="wilcox.test", comparisons = list(c("Diabetic", "Non-diabetic")), map_signif_level = significances) +
  scale_y_continuous(expand=expand_scale(.15,0)) +
  theme_bw() +
  scale_fill_manual(name='Medical condition', 
                    values=RColorBrewer::brewer.pal(3, 'Set1')) +
  theme(legend.position = 'none',
        plot.margin = margins) +
  xlab("Diabetes status") + ylab('PUFA/TFA Ratio')

pfg = (pf1 | pf2) + plot_annotation(tag_levels = list(c('f', 'g')), tag_prefix='(', tag_suffix=')')
saveRDS(pfg, file.path(figure_out_dir, 'figure7', 'figure7_fg.rds'))

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
write.csv(mydf3_save, file=file.path(table_out_dir, 'f7f-g.csv'), row.names = F)

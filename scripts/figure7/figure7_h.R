# main figure 7 panel h 
library(tidyverse)
library(reshape2)
library(ggpubr)
library(patchwork)

margins = margin(5.5,5.5,5.5,5.5, "pt")  # t r b l 
significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

source('scripts/files.R')

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

ph =
  corrsdf %>% 
  mutate(lab = as.character(round(rho,2))) %>% 
  mutate(sig = ' ',
         sig = ifelse(p.value < 0.05, '*', sig),
         sig = ifelse(p.value < 0.01, '**', sig),
         sig = ifelse(p.value < 0.001, '***', sig),
         sig = ifelse(p.value < 0.0001, '****', sig),
  ) %>%
  mutate(met=gsub('POFA', 'PUFA', met),
         met=gsub('MOFA', 'MUFA', met),
         met=gsub('Total', 'TFA', met)) %>% 
  mutate(met = factor(met, levels=c('PUFA', 'MUFA', 'PC', 'PC_TFA_Ratio', 'PUFA_TFA_Ratio', 'MUFA_TFA_Ratio', 'PUFA_MUFA_Ratio'))) %>% 
  mutate(param= factor(param, levels=rev(c('basal_metabolic_rate', 'cci', 'max_digits_remembered', 'walking_pace')))) %>%
  mutate(lbl = paste0(lab, sig)) %>% 
  ggplot(aes(x=met, y=param, fill=rho)) +
  geom_tile() +
  geom_text(aes(label=lbl), size=3) +
  labs(fill = 'Correlation\ncoeff.') +
  scale_fill_gradient2(high = 'firebrick2', mid = 'white', low = 'dodgerblue4', 
                       limits = c(-.42, .42), midpoint = 0)  +
  scale_x_discrete(name='Metabolites',
                   breaks=c('PUFA', 'MUFA', 'PC', 'PC_TFA_Ratio', 'PUFA_TFA_Ratio', 'MUFA_TFA_Ratio', 'PUFA_MUFA_Ratio'),
                   labels=c('PUFA', 'MUFA', 'PC', 'PC/TFA Ratio', 'PUFA/TFA Ratio', 'MUFA/TFA Ratio', 'PUFA/MUFA Ratio')) +
  scale_y_discrete(name='Health parameters',
                   breaks=c('basal_metabolic_rate', 'cci', 'max_digits_remembered', 'walking_pace'),
                   labels=c('Basal Metabolic Rate (n = 29,594)', 'CCI (n = 30,269)', 'Max Digits Remembered (n = 2,655)', 'Walking Pace (n = 29,941)')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        plot.margin = margins) +
  xlab('Metabolites') + ylab('Health parameters') + 
  plot_annotation(title='(h)')

saveRDS(ph, file.path(figure_out_dir, 'figure7', 'figure7_h.rds'))

## stat table
ph_stat_table = apply(combdf, 1, function(x){
  v1 = x[1] %>% as.character()
  v2 = x[2] %>% as.character()
  cres = cor.test(mydf4[,v1], mydf4[,v2], m='s')
  n = mydf4 %>% select(v1,v2) %>% drop_na %>% nrow()
  data.frame(param=v1, met=v2, 'rho'=cres$estimate, 'p.value'=cres$p.value, n=n)
}) %>%
  melt(measure.vars=c()) %>%
  select(-L1) %>%
  mutate(rho = as.numeric(rho),
         p.value = as.numeric(p.value)) %>%
  mutate(lab = as.character(round(rho,2))) %>%
  mutate(sig = ' ',
         sig = ifelse(p.value < 0.05, '*', sig),
         sig = ifelse(p.value < 0.01, '**', sig),
         sig = ifelse(p.value < 0.001, '***', sig),
         sig = ifelse(p.value < 0.0001, '****', sig),
  ) %>%
  mutate(met=gsub('POFA', 'PUFA', met),
         met=gsub('MOFA', 'MUFA', met),
         met=gsub('Total', 'TFA', met)) %>%
  mutate(met = factor(met, levels=c('PUFA', 'MUFA', 'PC', 'PC_TFA_Ratio', 'PUFA_TFA_Ratio', 'MUFA_TFA_Ratio', 'PUFA_MUFA_Ratio'))) %>%
  mutate(param = setNames(c('Basal Metabolic Rate', 'CCI', 'Max Digits Remembered', 'Walking Pace'),
                          c('basal_metabolic_rate', 'cci', 'max_digits_remembered', 'walking_pace'))[param]) %>%
  transmute('Health parameters' = param,
            'Metabolites' = met,
            'P-value' = p.value,
            'Correlation coefficient' = round(rho, 5),
            'Significance' = sig,
            'Sample size' = paste0('n=', n)
  ) %>%
  distinct()
xlsx::write.xlsx(ph_stat_table, file.path(table_out_dir, 'stat_table_fig7h.xlsx'))

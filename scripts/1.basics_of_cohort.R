# plat basics of the data (i.e., age dist, sex, vs.)
library(tidyverse)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

# age dist
cols = mymetadf %>% 
  filter(name %in% c('sex', 'age_at_rec')) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

mydf2 %>% 
  mutate(ageints = cut(as.numeric(age_at_rec), breaks = seq(30,80,by=5))) %>% 
  group_by(sex, ageints) %>% 
  summarise(n=n()) %>% 
  ggplot(aes(x=ageints, y=n, fill=sex)) +
  geom_bar(stat='identity', alpha=.7) +
  geom_text(aes(label=n), vjust=-.4) +
  facet_grid(sex~.) +
  theme_bw() +
  ylim(c(0, 4300)) +
  theme(legend.position='none') +
  xlab('Age interval') +
  ylab('Num. of individuals') +
  ggtitle('Age distribution')
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/age_dist.png', width=8, height=4)


mydf2 %>% group_by(sex) %>% summarise(n())
mydf2$age_at_rec %>% summary


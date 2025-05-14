library(tidyverse)

source('scripts/files.R')

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

ggsave(files.path(figure_out_dir, 'age_dist.png'), width=8, height=4)


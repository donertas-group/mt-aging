library(tidyverse)
library(reshape2)
library(ggsignif)

source('scripts/files.R')

significances = c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05)

# distribution of lactate levels
fields = c('age_at_rec', 'Lactate')
cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(cols$field_id) 
colnames(mydf2) = setNames(cols$name, cols$field_id)[colnames(mydf2)]

mydf3 = mydf2 %>% 
  select(Lactate, age_at_rec) %>% 
  drop_na()

mydf3 %>% 
  ggplot(aes(x=Lactate)) +
  geom_histogram(color='gray80', fill='lightblue4') +
  geom_vline(xintercept = quantile(mydf3$Lactate, probs = c(.25)), linetype=2, color='gray25') +
  annotate("text", x=2.9, y=4800, label = "Q1") +
  annotate("text", x=4.8, y=4500, label = "Q4") +
  geom_vline(xintercept = quantile(mydf3$Lactate, probs = c(.75)), linetype=2, color='gray25') +
  theme_bw() +
  labs(x='Lactate levels (mmol/l)', y='Count')
ggsave(files.path(figure_out_dir, 'figureSXX.pdf'), width=6, height=4)


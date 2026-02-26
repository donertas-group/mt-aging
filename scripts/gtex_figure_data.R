# this script is to create figure data files for each figure in .csv format. 
# 
# here we create a minimal table that can reproduce all figures, all other information for further analysis can be obtained using .rds files and $data object. 
# 
library(tidyverse)

p = readRDS('./results/alltissues/ExpressionLevel_AgeChangeRho.rds')

p$data %>%
  select(1,2,3,7) %>%
  write_excel_csv('./results/alltissues/ExpressionLevel_AgeChangeRho.csv')

p = readRDS('./results/change_w_age/adipose_subcutaneous.rds')

p$data %>%
  select(2,3,7) %>%
  write_excel_csv('./results/change_w_age/adipose_subcutaneous.csv')

p = readRDS('./results/change_w_age/adipose_vis.rds')

p$data %>%
  select(2,3,7) %>%
  write_excel_csv('./results/change_w_age/adipose_vis.csv')

p = readRDS('./results/change_w_age/breast.rds')

p$data %>%
  select(2,3,7) %>%
  write_excel_csv('./results/change_w_age/breast.csv')

p = readRDS('./results/change_w_age/expression_by_3q.rds')

p$data %>%
  select(2,3,7) %>%
  write_excel_csv('./results/change_w_age/expression_by_3q.csv')


p = readRDS('./results/change_w_age/ovary.rds')

p$data %>%
  select(2,3,7) %>%
  write_excel_csv('./results/change_w_age/ovary.csv')

p = readRDS('./results/meanexpression/above3Q_highlighted.rds')

p$data %>%
  select(1,2,3,7) %>%
  write_excel_csv('./results/meanexpression/above3Q_highlighted.csv')


# calculate charlson comorbidity index (cci) and append to main data frame.
library(reshape2)
library(tidyverse)
library(comorbidity)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all.rds')
metadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

# get field_ids for icd-10 illnesses
illnesses = metadf %>% 
  filter(name == 'illnesses') %>% 
  pull(field_id)

# create tidy df for comorbidity calculation
mynewdf = mydf %>% 
  select(f.eid, all_of(illnesses)) %>% 
  melt(id.var='f.eid') %>% 
  drop_na() %>% 
  transmute(f.eid, code=value)

# calculate cci
ccis = comorbidity(x=mynewdf, id = 'f.eid', code = 'code', map = 'charlson_icd10_quan', assign0 = TRUE)  
ccis_scores = score(ccis, weights = 'charlson', assign0 = TRUE)

# append to the main df
mydf2 = data.frame(f.eid = ccis$f.eid, cci = ccis_scores) %>% 
  right_join(mydf, by='f.eid') %>% 
  replace_na(list(cci = 0))
saveRDS(mydf2, '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')

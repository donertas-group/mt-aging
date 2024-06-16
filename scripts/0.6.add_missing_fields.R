# This script adds missing fields (basal metabolic rate and cognitive functions) to main dataframe
# This is to not wait for the new basket to be approved - after getting a new basket with all fields this will be obselete.
library(tidyverse)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_columns.rds')

myadddf = read_csv('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/scripts/get_missing_fields/additional_fields.csv')

# field names
cols_adddf = bind_rows(data.frame(f.eid = '399', name = 'n_incorr_matches'),
                      data.frame(f.eid = '10137', name = 'n_incorr_matches_pilot'),
                      data.frame(f.eid = '4282', name = 'max_digits_remembered'),
                      data.frame(f.eid = '20016', name = 'fluid_int_score'),
                      data.frame(f.eid = '23105', name = 'basal_metabolic_rate'))

#convert field names to match with main df
myaddmetadf = data.frame(var1=colnames(myadddf)) %>% 
  mutate(var = gsub('participant.', 'f.', var1)) %>% 
  mutate(var = gsub('_i', '.', var)) %>% 
  mutate(var = gsub('_a', '.', var)) %>% 
  mutate(var = gsub('p', '', var)) %>% 
  transmute(var1, f.eid = str_split_i(var, fixed('.'), 2),
            instance_id = str_split_i(var, fixed('.'), 3),
            array_id = str_split_i(var, fixed('.'), 4)) %>% 
  replace_na(list(array_id = '0')) %>% 
  left_join(cols_adddf, by='f.eid') %>% 
  transmute(field_id = paste0('f.', f.eid, '.', instance_id, '.', array_id),
            name, instance_id, array_id, f.eid, var1) %>% 
  filter(f.eid != 'eid') %>% 
  select(-f.eid)

# fix colnames of the data
colnames(myadddf) = setNames(c('f.eid', myaddmetadf$field_id), c('participant.eid', myaddmetadf$var1))[colnames(myadddf)]

#combine data and metadata and save
mydf2 = mydf %>% 
  left_join(myadddf, by='f.eid')
saveRDS(mydf2, '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all.rds')

mymetadf2 = bind_rows(select(myaddmetadf, -var1),
                      mymetadf) 

saveRDS(mymetadf2, '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')

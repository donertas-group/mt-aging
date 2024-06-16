# import required fields for mt-aging project from the whole table
library(tidyverse)
bd = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/data/raw/ukb678748_rdf.rds')

# get age & sex of inds
age_at_rec = colnames(bd)[which(startsWith(colnames(bd), 'f.21022'))]
sex = colnames(bd)[which(startsWith(colnames(bd), 'f.31.'))]

# get nmr metabolomics fields
# nmrdf = read.delim('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/Nightingale_biomarker_groups.txt') # to get all nmrs
# nmrs = sapply(mutate(nmrdf, field_id = paste0('f.', field_id))$field_id, function(x) colnames(bd)[which(startsWith(colnames(bd), as.character(x)))]) %>%
#   as.vector()
nmrs = c(
  'f.23437.0.0' = 'PC',
  'f.23446.0.0' = 'POFA',
  'f.23447.0.0' = 'MOFA',
  'f.23448.0.0' = 'SFA',
  'f.23458.0.0' = 'POFA_MOFA_Ratio',
  'f.23453.0.0' = 'POFA_Total_Ratio',
  'f.23454.0.0' = 'MOFA_Total_Ratio',
  'f.23455.0.0' = 'SFA_Total_Ratio',
  'f.23442.0.0' = 'TFA',
  'f.23471.0.0' = 'Lactate'
)


# get other related fields (grip strength etc.)
falls = colnames(bd)[which(startsWith(colnames(bd), 'f.2296.'))] # Falls in the last year
fractures_5_years = colnames(bd)[which(startsWith(colnames(bd), 'f.2463.'))] # Fractured/broken bones in last 5 years
fractures_falls = colnames(bd)[which(startsWith(colnames(bd), 'f.3005.'))] # Fracture resulting from simple fall
weight_change = colnames(bd)[which(startsWith(colnames(bd), 'f.2306.'))] # Weight change compared with 1 year ago

## health parameters
# cognitive measures
n_incorr_matches = colnames(bd)[which(startsWith(colnames(bd), 'f.399.'))] #Number of incorrect matches in round
n_incorr_matches_pilot = colnames(bd)[which(startsWith(colnames(bd), 'f.10137.'))] #Number of incorrect matches in round (pilot)
max_digits_remembered = colnames(bd)[which(startsWith(colnames(bd), 'f.4282.'))]  # Maximum digits remembered correctly
fluid_int_score = colnames(bd)[which(startsWith(colnames(bd), 'f.20016.'))] # fluid inteligence score

# other measures
basal_metabolic_rate = colnames(bd)[which(startsWith(colnames(bd), 'f.23105.'))] # basal metabolic rate
grip_strength_left = colnames(bd)[which(startsWith(colnames(bd), 'f.46.'))] # grip strength left hand
grip_strength_right = colnames(bd)[which(startsWith(colnames(bd), 'f.47.'))] # grip strength right hand
walking_pace = colnames(bd)[which(startsWith(colnames(bd), 'f.924.'))] # usual walking pace
# illnesses = colnames(bd)[which(startsWith(colnames(bd), 'f.20002.'))] # non-cancer illness
illnesses = colnames(bd)[which(startsWith(colnames(bd), 'f.41270.'))] # ICD10 disease codes
health_params = c(falls, fractures_5_years, fractures_falls, weight_change, 
                  n_incorr_matches,  n_incorr_matches_pilot, max_digits_remembered, fluid_int_score, 
                  grip_strength_left,  grip_strength_right,  walking_pace, basal_metabolic_rate, illnesses) 

# get eids with proteomics data
eids_proteomics = read.csv('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/olink_expressions.csv')$olink_instance_0.eid

## get metadf with column name codings
mymetadf = bind_rows(
  data.frame(field_id = c("f.eid", "f.21022.0.0", "f.31.0.0"),
             name = c('id', 'age_at_rec', 'sex')) %>% 
    mutate(instance_id = str_split_i(field_id, fixed('.'), 3),
           array_id = str_split_i(field_id, fixed('.'), 4)),
  data.frame(field_id = names(nmrs),
             name = unname(nmrs)) %>% 
    mutate(instance_id = str_split_i(field_id, fixed('.'), 3),
           array_id = str_split_i(field_id, fixed('.'), 4)),
  data.frame(c(
    paste0(falls, '-falls_in_last_year'),
    paste0(fractures_5_years, '-fractures_in_5_years'),
    paste0(fractures_falls, '-fractures_from_falls'),
    paste0(weight_change, '-weight_change'),
    paste0(n_incorr_matches, '-n_incorr_matches'),
    paste0(n_incorr_matches_pilot, '-n_incorr_matches_pilot'),
    paste0(max_digits_remembered, '-max_digits_remembered'),
    paste0(fluid_int_score, '-fluid_int_score'),
    paste0(basal_metabolic_rate, '-basal_metabolic_rate'),
    paste0(grip_strength_left, '-grip_strength_left'),
    paste0(grip_strength_right, '-grip_strength_right'),
    paste0(walking_pace, '-walking_pace'),
    paste0(illnesses, '-illnesses')
  )) %>% 
    magrittr::set_colnames('var') %>% 
    transmute(field_id = str_split_i(var, '-', 1),
              instance_id = str_split_i(field_id, fixed('.'), 3),
              array_id = str_split_i(field_id, fixed('.'), 4),
              name = str_split_i(var, '-', 2)) %>% 
    filter(instance_id == 0)
) 

# create df and only select inds with proteomics & nmr metabolomics data, as well as only instances 0
mydf =
  bd %>%
  select('f.eid', age_at_rec, sex, all_of(names(nmrs)), all_of(health_params)) %>% 
  filter(if_any(names(nmrs), ~ !is.na(.))) %>% 
  filter(`f.eid` %in% eids_proteomics) %>% 
  select(mymetadf$field_id)

saveRDS(mydf, '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df.rds')
saveRDS(mymetadf, '/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_columns.rds')

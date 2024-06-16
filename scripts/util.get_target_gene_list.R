#extract gene names from melikes file to be used in ukbb while extracting protemics expression
library(tidyverse)

basedir = "/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data"

target = read.delim(file.path(basedir, 'ukbb_blood_to_liver_correlated_genes.tsv'))
available = read.delim(file.path(basedir, 'samples_fields_list.txt'), sep = ',', header = FALSE)

mydf = available %>% 
  t() %>% 
  magrittr::set_colnames('field_name') %>% 
  as.data.frame() %>% 
  mutate(external_gene_name = gsub('olink_instance_0.', '', field_name)) %>% 
  filter(external_gene_name != 'eid') %>% 
  mutate(external_gene_name = str_to_upper(external_gene_name)) %>% 
  right_join(target, by='external_gene_name') 
  
write_csv(mydf, file.path(basedir, 'ukbb_available_targets.csv'))

# use these names to extract proteonmics data from ukbb rap
write_lines(paste0(c('olink_instance_0.eid', mydf$field_name), collapse=','), file=file.path(basedir, 'ukbb_available_targets_query.txt'))


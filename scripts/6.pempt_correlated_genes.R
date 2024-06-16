# Analysis of genes whose expression in blood is correlated with pempt expression in liver
library(tidyverse)
library(reshape2)

genecorrs = read.table('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/proteomics_meta/ukbb_blood_to_liver_correlated_genes.tsv', sep = '\t', header = T)
protein_exprr = read.csv('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/olink_expressions.csv')

colnames(protein_exprr) = str_to_upper(gsub('olink_instance_0.', '', colnames(protein_exprr)))
exprdf = protein_exprr %>% 
  melt(id.vars='EID') %>% 
  transmute(f.eid = EID,
            external_gene_name = variable,
            expression = value) %>% 
  full_join(genecorrs, by='external_gene_name') %>% 
  transmute(f.eid, gene = external_gene_name, expression, rho)

mydf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_final.rds')
mymetadf = readRDS('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed/ukb678748_subset_df_all_columns.rds')
fields = c('f.eid', 'sex', 'age_at_rec', 'PC', 'POFA', 'MOFA', 'TFA', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')

cols = mymetadf %>% 
  filter(name %in% fields) %>% 
  select(field_id, name)

mydf2 = mydf %>% select(c(cols$field_id), 'f.eid') 
colnames(mydf2) = setNames(c(cols$name, 'f.eid'), c(cols$field_id, 'f.eid'))[colnames(mydf2)]

mydf3 = mydf2 %>% 
  mutate(PC_TFA_Ratio = PC / TFA) %>% 
  rename(Age = age_at_rec)


# correlations
gnames = unique(exprdf$gene)
vnames = c('PC', 'POFA', 'MOFA', 'TFA', 'PC_TFA_Ratio', 'POFA_Total_Ratio', 'MOFA_Total_Ratio', 'POFA_MOFA_Ratio')

cordf = apply(expand.grid(gnames, vnames), 1, function(x){
  g = x[1] %>% as.character()
  v = x[2] %>% as.character()
  
  pexp = 
    exprdf %>% 
    filter(gene == g) %>% 
    select(f.eid, expression) %>% 
    distinct() %>% 
    full_join(select(mydf3, 'f.eid', all_of(v)), by='f.eid') %>% 
    drop_na() %>% 
    rename(var = eval(v))
  
  corres = cor.test(pexp$expression, pexp$var, m='s')
  c(g, v, corres$estimate, corres$p.value)
}) %>% 
  t() %>% 
  magrittr::set_colnames(c('gene', 'met', 'rho', 'p.value')) %>% 
  as.data.frame() %>% 
  mutate_at(.vars = vars(rho, p.value),
            .funs = list(as.numeric))

library(ComplexHeatmap)
cormat = cordf %>% 
  select(-p.value) %>% 
  pivot_wider(names_from=met, values_from=rho) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  t()

annodf = exprdf %>% 
  select(gene, rho) %>% 
  distinct() %>% 
  arrange(desc(rho))

annot_col = circlize::colorRamp2(c(-.43, 0, .43), c("#9970AB", "#F5F5F5", "#01665E"))
row_annots = HeatmapAnnotation('PEMT Corr.' = round(annodf$rho, 2), col = list('PEMT Corr.' = annot_col))

pdf('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/pemt_correlated_genes_corrs.pdf', width = 11, height = 4)
Heatmap(cormat[,annodf$gene], name='Correlation\ncoeff.', 
        col = circlize::colorRamp2(c(-.43, 0, .43), c("blue", "white", "red")), 
        top_annotation = row_annots, cluster_columns = FALSE)
dev.off()


## correlation across pemt correlated genes
protein_exprr %>% 
  drop_na() %>% 
  select(-EID) %>% 
  mutate_all(as.numeric) %>% 
  cor(m='s') %>% 
  as.data.frame() %>% 
  rownames_to_column('var1') %>% 
  melt(id.vars = 'var1') %>% 
  ggplot(aes(x=var1, y=variable, fill=value)) +
  geom_tile() +
  # geom_text(aes(label=value)) +
  labs(fill = 'Correlation\ncoeff.') +
  scale_fill_gradient2(high = 'firebrick2', mid = 'white', low = 'dodgerblue4', 
                       limits = c(-1, 1), midpoint = 0) +
  xlab('') + ylab('') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1)) +
  coord_fixed()
ggsave('/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/results/figures/corrs_among_pemt_correlated_genes.png', width = 8, height = 8)
# they all are positively correlated


# Description

# Setup ----

## Folder Organisation ----
# sapply(c('docs', 'data/processed', 'data/raw', 'scripts', 'results'),
#        function(folder){
#          system(paste('mkdir -p',folder))
#        }) # run only at the beginning of the project


## Libraries ----
easypackages::libraries(
  'tidyverse'
)

## Visualization helpers ----
easypackages::libraries('ggpubr')
theme_set(theme_pubr(base_size = 6))
pntnorm <- (1/0.352777778)
geom_point2 <- function(...)geom_point(shape=21,...)

# Path variables ----

gtex_tpm = '../../processed_data/external/bulkRNAseq/GTEx_v8/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
gtex_processed_folder = '../../processed_data/custom/bulkRNAseq/GTEx_v8/expression/tpm_lm_l2_qn/'
gtex_attr = '../../processed_data/custom/bulkRNAseq/GTEx_v8/attr.rds'

# Functions ----

plotsave <- function(ggobj, prefix, width=8, height=8, ...){
  path = strsplit(prefix,'/')[[1]]
  path = paste(path[-length(path)], collapse = '/')
  if(!file.exists(path)){
    system(paste('mkdir -p',path))
  }
  saveRDS(object = ggobj, file = paste(prefix,'.rds',sep=''))
  ggsave(file = paste(prefix, '.pdf', sep = ''), plot = ggobj, units = 'cm', 
         width = width, height = height, useDingbats = F, limitsize = F)
  ggsave(file = paste(prefix, '.png', sep = ''), plot = ggobj, units = 'cm', 
         width = width, height = height, limitsize = F)
}

# Analysis ----

## Variables ----

pemt = 'ENSG00000133027'

## Data ----
# GTEx v8 dataset was previously preprocessed to analyse age-related changes across tissues in Donertas 2021. We import this preprocessed data (TPM counts, corrected for batch effects using linear model, log2 transformed and quantile normalized). 

tissues = list.files(gtex_processed_folder)
attributes = readRDS(gtex_attr)
tissue_exp = lapply(tissues, function(mytis){
  tisname = gsub('.rds','',mytis)
  expr = readRDS(paste(gtex_processed_folder,mytis,sep='/'))
  pemt_gene_num = sum(rownames(expr) %in% pemt)
  if(length(pemt_gene_num)==1){
    data.frame(GeneID = pemt, sample_id = colnames(expr), expression = expr[pemt,], Tissue = tisname)
  } else if(length(pemt_gene_num)==0){
    NULL
  } else if(length(pemt_gene_num)>1){
    data.frame(GeneID = pemt, sample_id = colnames(expr), expression = colMeans(expr[pemt,]), Tissue = tisname)
  }
})

tissue_exp = reshape2::melt(tissue_exp, id.var=colnames(tissue_exp[[1]])) %>%
  select(-L1) %>%
  mutate(Gene = 'PEMT') %>%
  left_join(attributes)

## Age-related expression change ----

tissue_exp2 = tissue_exp %>%
  ungroup() %>%
  group_by(minor_tissue, age) %>%
  summarise(n = length(unique(id))) %>%
  filter(n>5) %>%
  ungroup() %>%
  left_join(tissue_exp) 

corvals = tissue_exp2 %>%
  group_by(Gene, minor_tissue) %>%
  summarise(rho = cor(expression, as.numeric(age), method = 's'),
            p = cor.test(expression, as.numeric(age), method = 's')$p.val) %>%
  mutate(padj = p.adjust(p, method='fdr'))

## Avg expression in tissues ----
# Since all tissues were processed separately and normalized within themselves, the mean expression levels are not informative. Instead we use the exact calculation GTEx provides in their website but limiting the sample set to only healthy individuals and removing outliers -  the same cohort we use for age-related expression change analysis -

pemtexp = data.table::fread(gtex_tpm, select = c('Name', unique(tissue_exp$sample_id)))

pemtexp = as.data.frame(pemtexp)
geneindex = which(sapply(strsplit(pemtexp$Name, '[.]'),function(x)x[1]%in%pemt))
rownames(pemtexp) = pemtexp$Name
pemtexp$Name = NULL
pemtexp = data.frame(expression = c(unlist(pemtexp[geneindex,])), sample_id = colnames(pemtexp))
pemtexp = pemtexp %>%
  left_join(attributes)

pemtexp_across_tissues = pemtexp %>% 
  group_by(id, sex, age, minor_tissue) %>%
  summarise(expression = median(expression)) %>%
  ungroup() 

### Mean expression plot ----

meanexpression_vals = pemtexp_across_tissues %>%
  group_by(minor_tissue) %>%
  summarise(avg_exp = mean(expression)) %>%
  mutate(abv3q = c('Lowest 75%','Highest 25%')[1+ (avg_exp>quantile(avg_exp,0.75))]) 

meanexpression_plot_abv3q = meanexpression_vals %>%
  left_join(pemtexp_across_tissues) %>%
  mutate(minor_tissue = fct_reorder(minor_tissue,-avg_exp)) %>%
  ggplot(aes(x = minor_tissue,  y = expression)) +
  geom_violin(scale = 'width', aes(fill = abv3q, color = abv3q)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = 'gray45', linewidth = 0.1) +
  scale_fill_brewer(palette = 2,direction = -1) +
  scale_color_brewer(palette = 2, direction = -1) +
  guides(fill = guide_legend('Expression', override.aes = list(color = NA)), color = 'none') +
  xlab(NULL) +
  ylab('TPM') +
  theme_pubr(base_size = 6) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.key.size = unit(0.25, 'cm'))
plotsave(meanexpression_plot_abv3q, './results/meanexpression/above3Q_highlighted',width = 16, height = 10)  

## Age-related change across all tissues ----

alltissues_mean_rho_violin = meanexpression_vals %>%
  left_join(pemtexp_across_tissues) %>%
  left_join(corvals) %>%
  mutate(minor_tissue = fct_reorder(minor_tissue,-avg_exp)) %>%
  ggplot(aes(y = expression)) +
  geom_violin(scale = 'width', aes(fill = rho, x = minor_tissue), linewidth = 0.1) +
  geom_boxplot(aes(x = minor_tissue), width = 0.2, outlier.shape = NA, color = 'gray45', linewidth = 0.1) +
  scale_fill_gradient2(midpoint = 0) +
  scale_color_gradient2(midpoint = 0) +
  guides(fill = guide_legend('Age-related\nexpression change', override.aes = list(color = NA)), color = 'none') +
  xlab(NULL) +
  ylab('TPM') +
  theme_pubr(base_size = 6) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = c(1,1), legend.justification = c(1,1),
        legend.key.size = unit(0.25, 'cm')) +
  geom_vline(xintercept = 12.5, linetype = 'dashed', size = 0.2) +
  annotate('text', x = 12, y = 500, hjust = 1, label = 'Hightest 25%', vjust = 1, size = 8/pntnorm) 

plotsave(alltissues_mean_rho_violin, './results/alltissues/ExpressionLevel_AgeChangeRho',width = 16, height = 12)  

corvals %>%
  left_join(meanexpression_vals) %>%
  filter(abv3q=='Highest 25%') %>%
  arrange(-avg_exp)

# Joining with `by = join_by(minor_tissue)`
# # A tibble: 12 × 7
# # Groups:   Gene [1]
# Gene  minor_tissue                       rho       p  padj avg_exp abv3q      
# <chr> <chr>                            <dbl>   <dbl> <dbl>   <dbl> <chr>      
# 1 PEMT  Liver                          -0.0471 0.638   0.821   126.  Highest 25%
# 2 PEMT  Breast - Mammary Tissue        -0.177  0.0480  0.245    66.7 Highest 25%
# 3 PEMT  Pituitary                       0.0298 0.693   0.827    65.5 Highest 25%
# 4 PEMT  Adipose - Subcutaneous         -0.197  0.00617 0.175    62.9 Highest 25%
# 5 PEMT  Testis                          0.0130 0.888   0.894    54.5 Highest 25%
# 6 PEMT  Thyroid                        -0.0424 0.560   0.791    50.3 Highest 25%
# 7 PEMT  Adipose - Visceral (Omentum)   -0.219  0.00760 0.175    44.6 Highest 25%
# 8 PEMT  Uterus                         -0.146  0.575   0.791    44.5 Highest 25%
# 9 PEMT  Prostate                        0.135  0.300   0.575    42.9 Highest 25%
# 10 PEMT  Ovary                          -0.333  0.0896  0.343    39.8 Highest 25%
# 11 PEMT  Nerve - Tibial                 -0.132  0.0742  0.310    36.0 Highest 25%
# 12 PEMT  Skin - Sun Exposed (Lower leg) -0.0169 0.804   0.881    35.4 Highest 25%

### Visualisation of specific tissues ----

adipose_sub_age = tissue_exp2 %>%
  filter(minor_tissue == 'Adipose - Subcutaneous') %>%
  group_by(minor_tissue, age) %>% summarise(avgexp = mean(expression)) %>%
  ungroup() %>%
  left_join(tissue_exp2) %>%
  ggplot(aes(x = age,  y= expression)) +
  geom_violin(scale = 'width', aes(fill = avgexp), color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.2) + 
  geom_jitter(width = 0.1, size = 0.05, color = 'gray45', alpha = 0.3) + 
  scale_fill_gradient(low = "#EFEDF5", high = '#6A51A3') +xlab(NULL) +
  ylab('Scaled Expression Level') +
  ggtitle('Adipose - Subcutaneous') +
  guides(fill = guide_colorbar('Mean\nScaled\nExpression')) +
  theme_pubr(base_size = 4.5) + 
  theme(legend.key.width = unit(0.25,'cm'),
        legend.key.height = unit(0.25,'cm'),
        # legend.position = c(0.025,1), legend.justification = c(0,1),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

plotsave(adipose_sub_age, './results/change_w_age/adipose_subcutaneous',width = 4, height = 4)  

breast_age = tissue_exp2 %>%
  filter(minor_tissue == 'Breast - Mammary Tissue') %>%
  group_by(minor_tissue, age) %>% summarise(avgexp = mean(expression)) %>%
  ungroup() %>%
  left_join(tissue_exp2) %>%
  ggplot(aes(x = age,  y= expression)) +
  geom_violin(scale = 'width', aes(fill = avgexp), color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.2) + 
  geom_jitter(width = 0.1, size = 0.05, color = 'gray45', alpha = 0.3) + 
  scale_fill_gradient(low = "#EFEDF5", high = '#6A51A3') +xlab(NULL) +
  ylab('Scaled Expression Level') +
  ggtitle('Breast - Mammary Tissue') +
  guides(fill = guide_colorbar('Mean\nScaled\nExpression')) +
  theme_pubr(base_size = 4.5) + 
  theme(legend.key.width = unit(0.25,'cm'),
        legend.key.height = unit(0.25,'cm'),
        # legend.position = c(0.025,1), legend.justification = c(0,1),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

plotsave(breast_age, './results/change_w_age/breast',width = 5, height = 4)  

ovary_age = tissue_exp2 %>%
  filter(minor_tissue == 'Ovary') %>%
  group_by(minor_tissue, age) %>% summarise(avgexp = mean(expression)) %>%
  ungroup() %>%
  left_join(tissue_exp2) %>%
  ggplot(aes(x = age,  y= expression)) +
  geom_violin(scale = 'width', aes(fill = avgexp), color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.2) + 
  geom_jitter(width = 0.1, size = 0.05, color = 'gray45', alpha = 0.3) + 
  scale_fill_gradient(low = "#EFEDF5", high = '#6A51A3') +xlab(NULL) +
  ylab('Scaled Expression Level') +
  ggtitle('Ovary') +
  guides(fill = guide_colorbar('Mean\nScaled\nExpression')) +
  theme_pubr(base_size = 4.5) + 
  theme(legend.key.width = unit(0.25,'cm'),
        legend.key.height = unit(0.25,'cm'),
        # legend.position = c(0.025,1), legend.justification = c(0,1),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

plotsave(ovary_age, './results/change_w_age/ovary',width = 3.5, height = 4)  

adipose_vis_age = tissue_exp2 %>%
  filter(minor_tissue == 'Adipose - Visceral (Omentum)') %>%
  group_by(minor_tissue, age) %>% summarise(avgexp = mean(expression)) %>%
  ungroup() %>%
  left_join(tissue_exp2) %>%
  ggplot(aes(x = age,  y= expression)) +
  geom_violin(scale = 'width', aes(fill = avgexp), color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.2) + 
  geom_jitter(width = 0.1, size = 0.05, color = 'gray45', alpha = 0.3) + 
  scale_fill_gradient(low = "#EFEDF5", high = '#6A51A3') +
  xlab(NULL) +
  ylab('Scaled Expression Level') +
  ggtitle('Adipose - Visceral') +
  guides(fill = guide_colorbar('Mean\nScaled\nExpression')) +
  theme_pubr(base_size = 4.5) + 
  theme(legend.key.width = unit(0.25,'cm'),
        legend.key.height = unit(0.25,'cm'),
        # legend.position = c(0.025,1), legend.justification = c(0,1),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

plotsave(adipose_vis_age, './results/change_w_age/adipose_vis',width = 4, height = 4)  

### Compare highly expressing and lowly expressing tissues ----

rho_3q = corvals %>%
  left_join(meanexpression_vals) %>%
  na.omit() %>%
  mutate(abv3q = factor(abv3q, levels = c('Lowest 75%','Highest 25%'))) %>%
  ggplot(aes(x = abv3q, y = rho)) +
  geom_violin(scale = 'width', color = NA, fill = 'gray85') +
  geom_boxplot(width = 0.2, linewidth = 0.2) +
  geom_jitter(width = 0.1, size = 0.2, color = 'gray45') + 
  stat_compare_means(size = 5/pntnorm) +
  xlab(NULL) +
  ylab('Age-Expression Correlation') + 
  theme_pubr(base_size = 4.5) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

plotsave(rho_3q, './results/change_w_age/expression_by_3q',width = 4, height = 4)  

## Combine plots for the main ----

allinset_allhoriz = ggarrange(rho_3q , breast_age, adipose_sub_age, 
                     ncol = 3, nrow = 1, widths = c(1.5,2,2), align = 'h')

plotsave(allinset_allhoriz, './results/inset', width = 16, height = 4)


allinset_allhoriz2 = ggarrange(rho_3q , breast_age + theme(legend.position = 'bottom',
                                                           legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0.3)), adipose_sub_age+ theme(legend.position = 'bottom',
                                                                                                                         legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0.3)), 
          ncol = 3, nrow = 1, widths = c(1.5,2,2), align = 'h')

alltis_w_inset = alltissues_mean_rho_violin + theme(legend.direction = 'horizontal',
                                                    legend.position = c(0,1), legend.justification = c(0,1),
                                                    legend.background = element_blank()) + 
  annotation_custom(ggplotGrob(allinset_allhoriz2 + theme(panel.border = element_rect(fill = NA, colour = 'black'))), xmin = 14, xmax = 48, 
                    ymin = 100, ymax = 550)

plotsave(alltis_w_inset, './results/alltis_w_inset',width = 16, height = 10)  

# Save All ----

save(list = ls(), file = './data/processed/all.RData')

# Session Info ----

sessionInfo()

# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggthemes_4.2.4  ggpubr_0.6.0    lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.1     purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] pkgload_1.3.2      bit64_4.0.5        vroom_1.6.1        carData_3.0-5      shiny_1.7.4        assertthat_0.2.1   cellranger_1.1.0   remotes_2.4.2      sessioninfo_1.2.2  pillar_1.9.0       backports_1.4.1   
# [12] glue_1.6.2         digest_0.6.31      RColorBrewer_1.1-3 promises_1.2.0.1   ggsignif_0.6.4     colorspace_2.1-0   cowplot_1.1.1      htmltools_0.5.5    httpuv_1.6.9       plyr_1.8.8         pkgconfig_2.0.3   
# [23] devtools_2.4.5     pheatmap_1.0.12    broom_1.0.4        xtable_1.8-4       scales_1.2.1       processx_3.8.0     later_1.3.0        tzdb_0.3.0         timechange_0.2.0   generics_0.1.3     farver_2.1.1      
# [34] car_3.1-2          usethis_2.1.6      ellipsis_0.3.2     cachem_1.0.7       withr_2.5.0        cli_3.6.1          readxl_1.4.2       magrittr_2.0.3     crayon_1.5.2       mime_0.12          memoise_2.0.1     
# [45] ps_1.7.4           fs_1.6.1           fansi_1.0.4        rstatix_0.7.2      pkgbuild_1.4.0     textshaping_0.3.6  profvis_0.3.7      tools_4.2.2        data.table_1.14.8  prettyunits_1.1.1  hms_1.1.3         
# [56] lifecycle_1.0.3    munsell_0.5.0      callr_3.7.3        compiler_4.2.2     systemfonts_1.0.4  rlang_1.1.0        grid_4.2.2         rstudioapi_0.14    htmlwidgets_1.6.2  miniUI_0.1.1.1     labeling_0.4.2    
# [67] gtable_0.3.3       abind_1.4-5        reshape2_1.4.4     R6_2.5.1           bit_4.0.5          fastmap_1.1.1      utf8_1.2.3         ragg_1.2.5         stringi_1.7.12     parallel_4.2.2     easypackages_0.1.0
# [78] Rcpp_1.0.10        vctrs_0.6.1        tidyselect_1.2.0   urlchecker_1.0.1  

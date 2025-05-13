# mt-aging

Collaborative project with Ermolaeva group. See the preprint: https://www.biorxiv.org/content/10.1101/2024.04.25.591184v1

## Data
UK Biobank dataset is used including phenotype data, a subset of proteomics and nmr metabolomics data. Processed data can be found under `/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/data/processed`, while the raw dataframe with phenotype/metadata can be found under `/scratch/shire/data/biobank/ukbb_immunosenescence/data/raw`.

As not all fields were available in the raw dataframe (i.e., cognitive measures), missing fields were seperately retreieved via DNAnexus, where the associated scripts can be found under `/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/scripts/get_missing_fields`. The proteomics expression data was also seperately downloaded from dnanexus using scripts under `/scratch/shire/data/biobank/ukbb_immunosenescence/mt-aging/scripts/get_proteomics_dnanexus`. Scripts under both folder need to be run on DNAnexus or on HPC using `dx toolkit`.

## Scripts

`./scripts` folder contains scripts that are ready to be published. `./myscripts` folder contains scripts that are used in data analysis and might contain fullpaths or sample names.

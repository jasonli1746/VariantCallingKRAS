# VariantCallingKRAS
Code to annotate VCF data with allele frequency, base calls

Written by Jason J. Li
code last changed 12/15/21

## Objective
There is evidence that KRAS allele frequency affects MAP kinase dependence in cancer and has a relationship with oncogenic inhibitor resistance. We seek to delve deeper into the details of this relationship. In order to do so, we need to calculate Variant Allele Frequency (VAF) for the KRAS mutation. The data we are using is from patient derived xenograft mouse models from 4 different vendors. We have run them through a variant calling pipeline to obtain VCF files. 

This repository represents code that will create a data frame where each row represents a PDX model. It will be annotated with VAF of KRAS and the mutation in HGVSp_short format (e.g., p.G12C). In cases where there are multiple mutations for a model, they will all be listed on the same row(model) and the VAFs and mutations will be separated by semicolons (i.e. "0.4;0.67" "p.G12C;p.G13D"). The dataframes for each vendor are generated in different sections and then combined into one giant dataframe at the end.

## How to use
The main script is "KRAS-annotate-muts-RNAseq-AF-clean.R". The rest of the R scripts are helper functions that are called in the body of the main script.

Data for now is in my "~/Projects/VariantCallRNAseq/data/" directory and as such this code cannot be run from cloning this repo. Apologies.

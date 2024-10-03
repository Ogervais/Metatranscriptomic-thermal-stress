---
title: "Test normalisation data metatranscriptomic"
author: "Ophelie Gervais"
date: "29/01/24"

---
  
  
#### Library --------------------------------------------------------------------------------

library('tidyverse')


############IMPORT###########

setwd("D:/Project/3_Temperature_effect/gene_expression/")

countdata <- read.table("featurecounts/featurecounts_eggnog_exp2.txt", header = TRUE, sep = "\t", row.names = 1)

#remove part of file name
colnames(countdata) <- gsub("\\_sorted.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("X.media.ecolo.home.ogervais.project.thermal_stress.metatranscripto.gene_expression.exp2.", "", colnames(countdata))
countdata <- countdata %>% 
  select(6:13) 
countdata <- countdata[,-5] #remove sample 24.2

#calculate library size for each sample
library_sizes <- colSums(countdata)

#calculate cpm for each sample
cpm_matrix <- t(t(countdata) / library_sizes) * 1e6
#or
#cpm_matrix <- cpm(countdata) #edgeR package

#get modify annotation file
annotation <- read.delim("eggnog_mmseqs.emapper.annotations", row.names=1) %>%
  separate(max_annot_lvl, into = c("number", "taxonomy"), sep = "\\|") %>%
  select(taxonomy) %>%
  mutate_all(~ ifelse(. == '-', NA, .)) 

#associate countdata with annot file
cpm <- merge(cpm_matrix, annotation, all.x = TRUE, by = "row.names") %>%
  filter(!is.na(taxonomy))

#within taxon scaling
within_taxon_scaled <- cpm %>%
  group_by(taxonomy) %>%
  mutate(across(starts_with("RNA"), ~ . / sum(.) * 1e6)) %>%
  column_to_rownames(var = 'Row.names') %>%
  select(1:7) %>%
  replace(is.na(.), 0)

#handling zero values 
df_filtered <- within_taxon_scaled[rowSums(within_taxon_scaled != 0, na.rm = TRUE) > 0, ]

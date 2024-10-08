#title: "data preparation kegg completion"
#author: "Ogervais"
#date: "2024-02-19"
  
#### Library --------------------------------------------------------------------------------

library('tidyverse')
library('ggbreak')

############IMPORT###########

setwd("D:/Project/3_Temperature_effect/gene_expression/exp2")

countdata <- read.table("featurecounts/featurecounts_eggnog_exp2.txt", header = TRUE, sep = "\t", row.names = 1)

#remove part of file name
colnames(countdata) <- gsub("\\_sorted.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("X.media.ecolo.home.ogervais.project.thermal_stress.metatranscripto.gene_expression.exp2.", "", colnames(countdata))
countdata <- countdata %>% 
  select(6:13) 
count15 <- countdata[,-5:-8] #remove sample 24
count24 <- countdata[,-1:-5] #remove sample 15 and 24.2

annotation <- read.delim("eggnog_mmseqs.emapper.annotations", row.names=1) %>%
  mutate_all(~ ifelse(. == '-', NA, .)) %>%
  separate(max_annot_lvl, into = c("number", "taxonomy"), sep = "\\|") %>%
  mutate_all(~ ifelse(. == '-', NA, .)) %>%
  select(taxonomy, KEGG_ko)
annot <- annotation 

#script that extract KEGG ko from a selected taxon
calculate_row_sums <- function(df, cols, species) {
  df_output <- df %>%
    mutate(RowSum = rowSums(select(., all_of(cols)))) %>%
    filter(RowSum != 0) %>%
    select(-all_of(cols)) %>%
    merge(annotation, all.x = TRUE, by = "row.names") %>%
    filter(!is.na(taxonomy) & !is.na(KEGG_ko)) %>% 
    separate_rows(KEGG_ko, sep = ",") %>%
    mutate(KEGG_ko = sub('ko:', '', KEGG_ko)) %>% 
    filter(taxonomy == species) %>%
    select(KEGG_ko) 
  return(df_output)
}

count15_spiro <- calculate_row_sums(count15, c('RNA.15.4', 'RNA.15.6', 'RNA.15.7', 'RNA.15.8'),  'Spirochaetes')
count15_bactero <- calculate_row_sums(count15, c('RNA.15.4', 'RNA.15.6', 'RNA.15.7', 'RNA.15.8'),  'Bacteroidetes')

count24_spiro <- calculate_row_sums(count24, c('RNA.24.3', 'RNA.24.3', 'RNA.24.8'),  'Spirochaetes')
count24_bactero <- calculate_row_sums(count24, c('RNA.24.3', 'RNA.24.3', 'RNA.24.8'),  'Bacteroidetes')

write.table(count15_spiro, file ='kegg_completion/count15_spiro' , sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(count15_bactero, file ='kegg_completion/count15_bactero' , sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(count24_spiro, file ='kegg_completion/count24_spiro' , sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(count24_bactero, file ='kegg_completion/count24_bactero' , sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#data were analysed using microbeAnnotator v.2.0.5 on an HPC
# python3 /media/ecolo/home/ogervais/anaconda2/envs/microbeannotator/lib/python3.7/site-packages/microbeannotator/pipeline/ko_mapper.py \
# -i ${INPUT}/count_spiro ${INPUT}/count_bactero\
# -p ${OUTPUT}/spiro_bactero





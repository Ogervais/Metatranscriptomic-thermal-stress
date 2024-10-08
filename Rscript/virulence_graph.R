#title: "virulence graph"
#author: "Ophelie Gervais"
#date: "9/07/24"

#### Library --------------------------------------------------------------------------------
library('tidyverse')
library('ggpubr')

##FUNCTION--------------------------------------------------------------

# Function to calculate mean excluding zeros
mean_excluding_zeros <- function(x) {
  mean(x[x != 0])
}


####normalised data
#used df_filtered from the normalization.R script
source("D:/Project/3_Temperature_effect/R_script/Github_script/Normalization.R")

##VIRULENCE FACTOR--------------------------------------------------------------
#separate temperature in two table and calculate the mean separately
countdata15 <- as.data.frame(df_filtered) %>% 
  select(1:4)
countdata15$Geneid <- rownames(countdata15)
countdata24 <- as.data.frame(df_filtered) %>% 
  select(5:7)
countdata24$Geneid <- rownames(countdata24)

##test with full dataset
#all the data
setwd('D:/Project/3_Temperature_effect/annotation/diamond/exp2/virulence')

datafull <- read.delim("eggnog_vir_bact_setB_noseq_eval5.tsv") %>%
  filter(length > 50) %>%
  filter(pident > 70) %>% #identity > 70%
  separate(sseqid, into = c("VFC", "other"), sep = "\\(")

#at 15C
data15_full <- datafull %>%
  left_join(countdata15, by = c('qseqid'= 'Geneid')) %>%
  select(1, 14:17) %>%
  column_to_rownames(var = 'qseqid') %>%
  filter(rowSums(.) != 0) #%>%
#merge(annotation, by = 'row.names', all.x = TRUE) %>%
#column_to_rownames(var = 'Row.names')
#write.table(data15_full, 'virulence_15.txt', sep = '\t', row.names = FALSE)

#at 24C
data24_full <- datafull %>%
  left_join(countdata24, by = c('qseqid'= 'Geneid')) %>%
  select(1, 14:16) %>%
  column_to_rownames(var = 'qseqid') %>%
  filter(rowSums(.) != 0) #%>%
#merge(annotation, by = 'row.names', all.x = TRUE) %>%
#column_to_rownames(var = 'Row.names')
#write.table(data24_full, 'virulence_24.txt', sep = '\t', row.names = FALSE)

#column_means_15_full <- colMeans(data15_full)
column_means_15_full <- data15_full %>%
  summarise(across(everything(), mean_excluding_zeros)) %>%
  unlist()
column_count_15_full <- apply(data15_full, 2, function(x) sum(x != 0))
# Convert the result to a data frame
means_data_15_full <- data.frame(
  Column = names(column_means_15_full),
  Mean = as.vector(column_means_15_full),
  Count = as.vector(column_count_15_full)
) 

#column_means_24_full <- colMeans(data24_full)
column_means_24_full <- data24_full %>%
  summarise(across(everything(), mean_excluding_zeros)) %>%
  unlist()
column_count_24_full <- apply(data24_full, 2, function(x) sum(x != 0))
# Convert the result to a data frame
means_data_24_full <- data.frame(
  Column = names(column_means_24_full),
  Mean = as.vector(column_means_24_full),
  Count = as.vector(column_count_24_full)
) 

means_data_full <- bind_rows(means_data_15_full, means_data_24_full)
temp <- c('15', '15', '15', '15', '24', '24', '24')
vir <- c('VFG', 'VFG', 'VFG', 'VFG', 'VFG', 'VFG', 'VFG')
means_data_full$Temperature <- temp
means_data_full$Condition <- vir



##ANTIMICROBIAL RESISTANCE--------------------------------------------------------------
setwd('D:/Project/3_Temperature_effect/annotation/diamond/exp2/amr')

annot_amr <- read.delim("card-data/aro_index.tsv") 

#data result diamond
#core dataset
data_amr <- read.delim("eggnog_amr_bact_noseq_eval5.tsv") %>%
  filter(length > 50) %>% #length >50
  filter(pident > 70) %>% #identity > 70%
  separate(sseqid, into = paste0("col", 1:4), sep = "\\|", fill = "right", extra = "drop") %>%
  left_join(annot_amr, by = c("col3" = "ARO.Accession"))

#calculate mean and count of AMR
#at 15C
data15_amr <- data_amr %>%
  left_join(countdata15, by = c('qseqid'= 'Geneid')) %>%
  select(1, 27:30) %>%
  column_to_rownames(var = 'qseqid') %>%
  filter(rowSums(.) != 0) #%>%
#merge(annotation, by = 'row.names', all.x = TRUE) %>%
#column_to_rownames(var = 'Row.names')
#write.table(data15_amr, 'amr_15.txt', sep = '\t', row.names = FALSE)

#at 24C
data24_amr <- data_amr %>%
  left_join(countdata24, by = c('qseqid'= 'Geneid')) %>%
  select(1, 27:29) %>%
  column_to_rownames(var = 'qseqid') %>%
  filter(rowSums(.) != 0) #%>%
#merge(annotation, by = 'row.names', all.x = TRUE) %>%
#column_to_rownames(var = 'Row.names')
#write.table(data24_amr, 'amr_24.txt', sep = '\t', row.names = FALSE)


#column_means_15_amr <- colMeans(data15_amr)
column_means_15_amr <- data15_amr %>%
  summarise(across(everything(), mean_excluding_zeros)) %>%
  unlist()
column_count_15_amr <- apply(data15_amr, 2, function(x) sum(x != 0))
# Convert the result to a data frame
means_data_15_amr <- data.frame(
  Column = names(column_means_15_amr),
  Mean = as.vector(column_means_15_amr),
  Count = as.vector(column_count_15_amr)
) 

#column_means_24_amr <- colMeans(data24_amr)
column_means_24_amr <- data24_amr %>%
  summarise(across(everything(), mean_excluding_zeros)) %>%
  unlist()
column_count_24_amr <- apply(data24_amr, 2, function(x) sum(x != 0))
# Convert the result to a data frame
means_data_24_amr <- data.frame(
  Column = names(column_means_24_amr),
  Mean = as.vector(column_means_24_amr),
  Count = as.vector(column_count_24_amr)
) 

means_data_amr <- bind_rows(means_data_15_amr, means_data_24_amr)
temp <- c('15', '15', '15', '15', '24', '24', '24')
vir <- c('ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG')
means_data_amr$Temperature <- temp
means_data_amr$Condition <- vir
means_data <- bind_rows(means_data_full, means_data_amr)

setwd('D:/Project/3_Temperature_effect/annotation/diamond/exp2/graph')

graph_count <- ggplot(means_data, aes(Temperature, Count, fill = Temperature)) +
  geom_boxplot() +
  geom_point() +
  theme_minimal() +
  scale_fill_manual(values = c("#009dfe", "#fb3455")) +
  facet_wrap(~Condition) +
  ylab('Number of different genes')
graph_count
#ggsave('vir_amr_count_normalized.png', dpi = 500, height =2000, width = 2500, units = 'px')


graph_tmp <- ggplot(means_data, aes(Temperature, Mean, fill = Temperature)) +
  geom_boxplot() +
  geom_point() +
  theme_minimal() +
  scale_fill_manual(values = c("#009dfe", "#fb3455")) +
  facet_wrap(~Condition) +
  ylab('TPM')
graph_tmp
#ggsave('vir_amr_tpm_normalized.png', dpi = 500, height =2000, width = 2500, units = 'px')


ggarrange(graph_count, graph_tmp,
          nrow = 1, ncol = 2,
          labels = c('A', 'B'))

#ggsave('vir_amr_count_tpm_normalized.png', dpi = 500, height =2000, width = 4000, units = 'px')

## statistical test
#pairwise
# Mann-Whitney U test
wilcox.test(means_data_15_amr$Mean, means_data_24_amr$Mean)
wilcox.test(means_data_15_amr$Count, means_data_24_amr$Count)
wilcox.test(means_data_15_full$Mean, means_data_24_full$Mean)
wilcox.test(means_data_15_full$Count, means_data_24_full$Count)

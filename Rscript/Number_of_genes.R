#### Library --------------------------------------------------------------------------------

library('tidyverse')
library('ggbreak')

############IMPORT###########

setwd("D:/Project/3_Temperature_effect/gene_expression/exp2")

#transform name
annotation <- read.delim("eggnog_mmseqs.emapper.annotations", row.names=1) %>%
  mutate_all(~ ifelse(. == '-', NA, .)) %>%
  separate(max_annot_lvl, into = c("number", "taxonomy"), sep = "\\|") %>%
  mutate_all(~ ifelse(. == '-', NA, .)) %>%
  select(taxonomy, COG_category)
# annotation_ <- annotation %>%
#   separate_rows(COG_category, sep = '(?<=.)(?=.)')
# annot <- annotation_
annot <- annotation

annot$taxonomy <- fct_recode(annot$taxonomy, 'Other phylum' = 'Bacteria',
                             'Acidobacteriota' = 'Acidobacteriia',
                             'Acidobacteriota' = 'Acidobacteria',
                             'Actinomycetota' = 'Actinobacteria',
                             'Actinomycetota' = 'Acidimicrobiia',
                             'Actinomycetota' = 'Rubrobacteria',
                             'Actinomycetota' = 'Coriobacteriia',
                             'Bacteroidota' = 'Bacteroidetes',
                             'Chloroflexota' = 'Chloroflexi',
                             'Chloroflexota' = 'Chloroflexia',
                             'Chloroflexota' = 'Dehalococcoidia',
                             'Bacillota' = 'Firmicutes',
                             'Bacillota' = 'Bacilli',
                             'Bacillota' = 'Clostridia',
                             'Bacillota' = 'Negativicutes',
                             'Bacillota' = 'Erysipelotrichia',
                             'Deferribacterota' = 'Deferribacteres',
                             'Spirochaetota' = 'Spirochaetes',
                             'Synergistota' = 'Synergistetes',
                             'Verrucomicrobiota' = 'Verrucomicrobiae',
                             'Verrucomicrobiota' = 'Verrucomicrobia',
                             'Verrucomicrobiota' = 'Opitutae',
                             'Pseudomonadota' = 'Proteobacteria',
                             'Pseudomonadota' = 'Acidithiobacillales',
                             'Pseudomonadota' = 'Alphaproteobacteria',
                             'Pseudomonadota' = 'Caulobacterales',
                             'Pseudomonadota' = 'Rhodospirillales',
                             'Pseudomonadota' = 'Sphingomonadales',
                             'Pseudomonadota' = 'Betaproteobacteria',
                             'Pseudomonadota' = 'Neisseriales',
                             'Pseudomonadota' = 'Nitrosomonadales',
                             'Pseudomonadota' = 'Rhodocyclales',
                             'Pseudomonadota' = 'Deltaproteobacteria',
                             'Pseudomonadota' = 'Epsilonproteobacteria',
                             'Pseudomonadota' = 'Gammaproteobacteria',
                             'Pseudomonadota' = 'Aeromonadales',
                             'Pseudomonadota' = 'Chromatiales',
                             'Pseudomonadota' = 'Legionellales',
                             'Pseudomonadota' = 'Methylococcales',
                             'Pseudomonadota' = 'Oceanospirillales',
                             'Pseudomonadota' = 'Pasteurellales',
                             'Pseudomonadota' = 'Thiotrichales',
                             'Pseudomonadota' = 'Vibrionales',
                             'Pseudomonadota' = 'Xanthomonadales',
                             'Pseudomonadota' = 'Hydrogenophilales',
                             'Pseudomonadota' = 'Rickettsiales',
                             'Thermodesulfobacteriota' = 'Thermodesulfobacteria',
                             'Thermomicrobiota' = 'Thermomicrobia',
                             'Thermotogota' = 'Thermotogae',
                             'Mycoplasmota' = 'Tenericutes',
                             'Nitrospirota' = 'Nitrospirae',
                             'Deinococcota' = 'Deinococcus-Thermus',
                             'Bdellovibrionota' = 'Bdellovibrionales',
                             'Aquificota' = 'Aquificae') 

##Number of genes at 15C and 24--------------------------------------------------------

countdata <- read.table('featurecounts/featurecounts_eggnog_exp2.txt', row.names = 1, header = TRUE)
#remove part of file name
colnames(countdata) <- gsub("\\_sorted.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("X.media.ecolo.home.ogervais.project.thermal_stress.metatranscripto.gene_expression.exp2.", "", colnames(countdata))
countdata <- countdata %>% 
  select(6:13) 
countdata <- countdata[,-5] #remove sample 24.2

countdata_phylum <- merge(countdata, annot, all.x = TRUE, by = "row.names") %>%
  filter(!is.na(taxonomy))

#extract the count of each phylum at 15°C
count_15_phylum <- countdata_phylum %>%
  select(Row.names, RNA.15.4, RNA.15.6, RNA.15.7, RNA.15.8, taxonomy) %>%
  filter(!(RNA.15.4 == 0 & RNA.15.6 == 0 & RNA.15.7== 0 &  RNA.15.8 == 0))

phylum <- c('Pseudomonadota', 'Spirochaetota', 'Bacteroidota', 'Bacillota', 'Other phylum', 'Actinomycetota', 'Cyanobacteria', 'Planctomycetes', 
            'Chloroflexota', 'Thermotogota', 'Verrucomicrobiota', 'Fusobacteria', 'Acidobacteriota', 'Deinococcota', 'Aquificota', 'Chlorobi', 'Chlamydiae', 
            'Deferribacterota', 'Bdellovibrionota', 'Nitrospirota',
            'Synergistota', 'Thermodesulfobacteriota', 'Gemmatimonadetes', 'Thermomicrobiota', 'Mycoplasmota')

mutate(count_15_phylum, taxonomy = fct_relevel(taxonomy, phylum)) %>%
  ggplot(aes(x = taxonomy)) +
  geom_bar() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.8)) +
  ylab('Number of genes') +
  scale_y_break(breaks = c(1000, 5000), scales = 0.4)

#extract the count of each phylum at 24°C
count_24_phylum <- countdata_phylum %>%
  select(Row.names, RNA.24.3, RNA.24.4, RNA.24.8, taxonomy) %>%
  filter(!(RNA.24.3 == 0 & RNA.24.4 == 0 & RNA.24.8 == 0)) 

mutate(count_24_phylum, taxonomy = fct_relevel(taxonomy, phylum)) %>%
  ggplot(aes(x = taxonomy)) +
  geom_bar() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.8)) +
  ylab('Number of genes') +
  scale_y_break(breaks = c(1000, 6000), scales = 0.4)

#Get the sum of count at 15°C
sum_15 <- count_15_phylum %>%
  group_by(taxonomy) %>%
  summarise(count15 = n()) %>%
  ungroup()

#Get the sum of count at 24°C
sum_24 <- count_24_phylum %>%
  group_by(taxonomy) %>%
  summarise(count24 = n()) %>%
  ungroup()

sum <- merge(sum_15, sum_24, all.x = TRUE, by = 'taxonomy') %>%
  pivot_longer(cols = count15:count24, names_to = 'condition', values_to = 'value') 

alpha_vals <- c(0.4, 1)

mutate(sum, taxonomy = fct_relevel(taxonomy, phylum)) %>%
  ggplot(aes(x = taxonomy, y = value, alpha = condition)) +
  geom_bar(stat = 'identity', position=position_dodge()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.8),
    legend.position = 'none') +
  ylab('Number of different genes') +
  scale_alpha_manual('Temperature', values = alpha_vals, labels = c('15°C', '24°C')) +
  scale_y_break(breaks = c(850, 6000), scales = 0.3)
ggsave('graph/number_genes_phylum_temp.png', dpi = 500, height =2000, width = 3800, units = 'px')


##extract count for the pseudomonadota phylum-----------------------------------
countdata_pseudo <- merge(countdata, annotation, all.x = TRUE, by = "row.names") %>%
  filter(!is.na(taxonomy)) %>%
  filter(grepl('Proteobacteria|Acidithiobacillales|Alphaproteobacteria|Caulobacterales|Rhodospirillales|Sphingomonadales|Betaproteobacteria|Neisseriales|Nitrosomonadales|Rhodocyclales|Deltaproteobacteria|Epsilonproteobacteria|Gammaproteobacteria|Aeromonadales|Chromatiales|Legionellales|Methylococcales|Oceanospirillales|Pasteurellales|Thiotrichales|Vibrionales|Xanthomonadales|Hydrogenophilales|Rickettsiales', taxonomy))
countdata_pseudo$taxonomy <- fct_recode(countdata_pseudo$taxonomy,
                                        'Other Pseudomonadota' = 'Proteobacteria',
                                        'Other Alphaproteobacteria' = 'Alphaproteobacteria',
                                        'Other Betaproteobacteria' = 'Betaproteobacteria',
                                        'Other Gammaproteobacteria' = 'Gammaproteobacteria')

#extract the count of each pseudomononadota class at 15°C
countdata_pseudo_15 <- countdata_pseudo %>%
  select(Row.names, RNA.15.4, RNA.15.6, RNA.15.7, RNA.15.8, taxonomy) %>%
  filter(!(RNA.15.4 == 0 & RNA.15.6 == 0 & RNA.15.7== 0 &  RNA.15.8 == 0))

mutate(countdata_pseudo_15, taxonomy = fct_relevel(taxonomy, pseudo)) %>%
  ggplot(aes(x = taxonomy, fill = taxonomy)) +
  geom_bar() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, colour = pal),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.8)) +
  ylab('Number of different genes') +
  scale_fill_manual(values = pal) + 
  theme(legend.position='none') +
  scale_y_break(breaks = c(1500, 2000), scales = 0.3)

#extract the count of each pseudomononadota class at 24°C
countdata_pseudo_24 <- countdata_pseudo %>%
  select(Row.names, RNA.24.3, RNA.24.4, RNA.24.8, taxonomy) %>%
  filter(!(RNA.24.3 == 0 & RNA.24.4 == 0 & RNA.24.8 == 0))

pal <- c("#A034F0", "#A034F0", "#A034F0", "#A034F0", "#A034F0", "#A034F0", "#A034F0", "#A034F0", "#A034F0", "#A034F0", '#0bb4ff', '#0bb4ff', '#0bb4ff', '#0bb4ff', '#0bb4ff', "#FF8C00", "#FF8C00", "#FF8C00", "#FF8C00",  '#555555', '#555555', '#555555', '#555555', '#555555') 

pseudo<- c('Vibrionales', 'Other Gammaproteobacteria', 'Xanthomonadales', 'Oceanospirillales', 'Chromatiales', 'Thiotrichales', 'Methylococcales', 'Legionellales', 'Aeromonadales', 'Pasteurellales', 
           'Other Alphaproteobacteria', 'Rhodospirillales', 'Sphingomonadales', 'Rickettsiales', 'Caulobacterales',
           'Other Betaproteobacteria', 'Neisseriales', 'Rhodocyclales', 'Nitrosomonadales',
           'Deltaproteobacteria', 'Other Pseudomonadota', 'Epsilonproteobacteria', 'Acidithiobacillales', 'Hydrogenophilales')

mutate(countdata_pseudo_24, taxonomy = fct_relevel(taxonomy, pseudo)) %>%
  ggplot(aes(x = taxonomy, fill = taxonomy)) +
  geom_bar() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, colour = pal),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.8)) +
  ylab('Number of different genes') +
  scale_fill_manual(values = pal) + 
  theme(legend.position='none') +
  scale_y_break(breaks = c(2000, 3000), scales = 0.3)

#Get the sum of count at 15°C
sum_15_pseudo <- countdata_pseudo_15 %>%
  group_by(taxonomy) %>%
  summarise(count15 = n()) %>%
  ungroup()

#Get the sum of count at 24°C
sum_24_pseudo <- countdata_pseudo_24 %>%
  group_by(taxonomy) %>%
  summarise(count24 = n()) %>%
  ungroup()
sum_pseudo <- merge(sum_15_pseudo, sum_24_pseudo, all.x = TRUE, by = 'taxonomy') %>%
  pivot_longer(cols = count15:count24, names_to = 'condition', values_to = 'value') 

alpha_vals <- c(0.4, 1)

mutate(sum_pseudo, taxonomy = fct_relevel(taxonomy, pseudo)) %>%
  ggplot(aes(x = taxonomy, y = value, alpha = condition, fill = taxonomy)) +
  geom_bar(stat = 'identity', position=position_dodge()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, colour = pal),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.8)) +
  ylab('Number of different genes') +
  scale_alpha_manual('Temperature', values = alpha_vals, labels = c('15°C', '24°C')) +
  scale_fill_manual(values = pal) +
  guides(fill = FALSE) +
  scale_y_break(breaks = c(2500, 3000), scales = 0.3)
#ggsave('graph/number_genes_pseudomonadota_temp.png', dpi = 500, height = 2000, width = 4000, units = 'px')









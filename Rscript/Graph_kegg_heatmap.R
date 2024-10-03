---
title: "KEGG heatmap"
author: "Ophelie Gervais"
date: "31/01/24"
---
  
#### Library --------------------------------------------------------------------------------

library('tidyverse')
library('ggpubr')
library('RColorBrewer')
library('viridis')
library('scales')

############IMPORT###########

setwd("D:/Project/3_Temperature_effect/gene_expression/exp2")

annot_dram <- read_tsv('annotations_dram.tsv') %>%
  separate(name, into = c("merge", "name"), sep = "^(?:[^_]*_){2}") %>%
  select(2:24) %>%
  column_to_rownames(var = 'name')

countdata <- read.table("result/deseq2/15_vs_24_normalized_DE.txt", header = TRUE, sep = "\t", row.names = 1)
countdata <- merge(countdata, annot_dram, all.x = TRUE, by = "row.names")

#modify file and kept only interesting taxa
data <- countdata %>%
  select(log2FoldChange, taxonomy, Description, Preferred_name, EC, KEGG_ko, KEGG_Pathway, KEGG_Module, ko_id, kegg_hit) %>%
  mutate_all(~ ifelse(. == '-', NA, .)) %>%
  filter(grepl('Spirochaetes|Alphaproteobacteria|Caulobacterales|Rhodospirillales|Sphingomonadales|Gammaproteobacteria|Aeromonadales|Chromatiales|Legionellales|Methylococcales|Oceanospirillales|Pasteurellales|Thiotrichales|Vibrionales|Xanthomonadales|Bacteroidetes', taxonomy)) 
data$taxonomy <- fct_recode(data$taxonomy, 
                            'Spirochaetota' = 'Spirochaetes',
                            'Others Alphaproteobacteria' = 'Alphaproteobacteria',
                            'Others Gammaproteobacteria' = 'Gammaproteobacteria',
                            'Bacteroidota' = 'Bacteroidetes')

# Split the "KEGG_Module" column by ','
data_new <- data %>%
  separate_rows(KEGG_Module, sep = ",")

complete_table_brite_sum <- data_new %>% 
  group_by(KEGG_Module, taxonomy) %>% 
  summarize(mean(log2FoldChange))

colnames(complete_table_brite_sum) <- c(
  "Module",
  "taxonomy",
  "log2FoldChange")

#read names of modules
curated_names <- read.table("D:/Project/7_Database/KEGG_module_bacteria_.txt", header = FALSE, sep = "\t")

sel.gmm <- subset(curated_names, V1 %in% unique(complete_table_brite_sum$Module))

global.table <- complete_table_brite_sum %>% 
  left_join(sel.gmm, by = c("Module" = "V1"))
colnames(global.table) <- c("Module", "taxonomy", "log2FoldChange", "Names", "Category")
global.table <- na.omit(global.table)
global.table <- as.data.frame(global.table)

#write.table(global.table, file="result/Supplementary_KEGG.txt", sep="\t", row.names = FALSE)

order_cat = c('Metabolic capacity', 'Drug resistance', 'Aromatics degradation', 'Biosynthesis of other bacterial compounds', 
              'Polyketide sugar unit biosynthesis', 'Terpenoid backbone biosynthesis', 'Cofactor and vitamin metabolism', 'Lipopolysaccharide metabolism',
              'Glycosaminoglycan metabolism', 'Other amino acid metabolism', 'Aromatic amino acid metabolism', 'Histidine metabolism', 
              'Polyamine biosynthesis', 'Arginine and proline metabolism', 'Lysine metabolism', 'Branched-chain amino acid metabolism', 
              'Cysteine and methionine metabolism', 'Serine and threonine metabolism', 'Pyrimidine metabolism', 'Purine metabolism', 
              'Lipid metabolism', 'Fatty acid metabolism','ATP synthesis', 'Photosynthesis', 'Sulfur metabolism', 'Nitrogen metabolism', 
              'Methane metabolism', 'Carbon fixation', "Other carbohydrate metabolism", "Central carbohydrate metabolism") 
order_tax = c('Bacteroidota', 'Others Alphaproteobacteria', 'Sphingomonadales', 
              'Others Gammaproteobacteria', 'Chromatiales', 'Legionellales', 'Methylococcales','Oceanospirillales', 'Thiotrichales', 'Vibrionales', 'Xanthomonadales', 
              'Spirochaetota')
order_mod = c('Nitrate assimilation', 
              'Vancomycin resistance, D-Ala-D-Lac type', 'Multidrug resistance, efflux pump MexAB-OprM', 'Cationic antimicrobial peptide (CAMP) resistance, protease PgtE', 
              'Trans-cinnamate degradation', 'Toluene degradation', 'Catechol ortho-cleavage', 'Catechol meta-cleavage', 'Benzoate degradation, benzoate => catechol / methylbenzoate => methylcatechol', 'Anthranilate degradation',
              'Rebeccamycin biosynthesis', 'Pyrrolnitrin biosynthesis',
              'dTDP-L-rhamnose biosynthesis',
              'C5 isoprenoid biosynthesis', 'C10-C20 isoprenoid biosynthesis', 
              'Ubiquinone biosynthesis', 'Thiamine biosynthesis, prokaryotes, AIR (+ DXP/tyrosine) => TMP/TPP', 'Tetrahydrofolate biosynthesis, mediated by ribA and trpF', 
              'Tetrahydrofolate biosynthesis, mediated by PTPS', 'Tetrahydrofolate biosynthesis, GTP => THF', 'Pyridoxal-P biosynthesis, erythrose-4P => pyridoxal-P', 'Pimeloyl-ACP biosynthesis', 'Riboflavin biosynthesis', 'Pantothenate biosynthesis, valine/L-aspartate => pantothenate',
              'NAD biosynthesis, aspartate', 'Menaquinone biosynthesis, chorismate (+ polyprenyl-PP) => menaquinol', 'Heme biosynthesis, plants and bacteria, glutamate => heme',
              'Coenzyme A biosynthesis', 'Cobalamin biosynthesis, cobyrinate a', 'C1-unit interconversion', 'Biotin biosynthesis, pimeloyl-ACP/CoA => biotin', 'Biotin biosynthesis, BioW pathway', 'Biotin biosynthesis, BioI pathway',
              'KDO2-lipid A biosynthesis, Raetz pathway, LpxL-LpxM type', 'CMP-KDO biosynthesis', 'ADP-L-glycero-D-manno-heptose biosynthesis',
              'Keratan sulfate degradation',
              'Glutathione biosynthesis',
              'Tyrosine degradation', 'Tyrosine biosynthesis, chorismate => HPP => tyrosine','Tryptophan biosynthesis', 'Shikimate pathway',
              'Histidine biosynthesis', 
              'Polyamine biosynthesis, arginine => agmatine => putrescine => spermidine', 'Polyamine biosynthesis, arginine => ornithine => putrescine', 'GABA biosynthesis, prokaryotes',
              'Urea cycle', 'Proline biosynthesis', 'Ornithine biosynthesis, glutamate', 'Arginine biosynthesis, glutamate', 'Arginine biosynthesis, ornithine',
              'Lysine degradation', 'Lysine biosynthesis, succinyl-DAP pathway', 'Lysine biosynthesis, DAP dehydrogenase pathway','Lysine biosynthesis, DAP aminotransferase pathway', 'Lysine biosynthesis, acetyl-DAP pathway', 
              'Valine/isoleucine biosynthesis', 'Leucine degradation', 'Leucine biosynthesis','Isoleucine biosynthesis, threonine', 'Isoleucine biosynthesis, pyruvate',
              'Methionine salvage pathway', 'Methionine degradation', 'Methionine biosynthesis', 'Ethylene biosynthesis', 'Cysteine biosynthesis, serine', 'Cysteine biosynthesis, methionine',
              'Threonine biosynthesis', 'Serine biosynthesis', 'Ectoine biosynthesis',
              'Pyrimidine ribonucleotide biosynthesis', 'De novo pyrimidine biosynthesis',
              'Guanine ribonucleotide biosynthesis', 'Deoxyribonucleotide biosynthesis', 'De novo purine biosynthesis', 'Adenine ribonucleotide biosynthesis', 
              'Triacylglycerol biosynthesis', 'Ketone body biosynthesis', 'Jasmonic acid biosynthesis',
              'Fatty acid biosynthesis', 'beta-Oxidation, acyl-CoA synthesis', 'beta-Oxidation', 
              'Succinate dehydrogenase', 'NADH:quinone oxidoreductase', 'F-type ATPase', 'Cytochrome o ubiquinol oxidase', 'Cytochrome c oxidase, prokaryotes', 'Cytochrome c oxidase, cbb3-type', 'Cytochrome c oxidase',
              'Cytochrome bc1 complex respiratory unit', 'Cytochrome bc1 complex',
              'Anoxygenic photosystem II', 
              'Thiosulfate oxidation by SOX complex', 'Dissimilatory sulfate reduction', 'Denitrification', 'Assimilatory sulfate reduction', 
              'Dissimilatory nitrate reduction', 'Complete nitrification', 'Assimilatory nitrate reduction', 
              'Methanogenesis, methylamine/dimethylamine/trimethylamine => methane','Methanogenesis, acetate => methane', 'Formaldehyde assimilation, serine pathway',
              'Reductive pentose phosphate cycle, ribulose-5P => glyceraldehyde-3P', 'Reductive pentose phosphate cycle, glyceraldehyde-3P => ribulose-5P', 'Reductive pentose phosphate cycle (Calvin cycle)',
              'Reductive citrate cycle (Arnon-Buchanan cycle)', 'Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)', 'Phosphate acetyltransferase-acetate kinase pathway', 
              'Incomplete reductive citrate cycle', 'Hydroxypropionate-hydroxybutylate cycle', 'Dicarboxylate-hydroxybutyrate cycle', 'CAM (Crassulacean acid metabolism), light', 
              'CAM (Crassulacean acid metabolism), dark', 'C4-dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type', 'C4-dicarboxylic acid cycle, NAD - malic enzyme type', 
              'C4-dicarboxylic acid cycle, NADP - malic enzyme type', '3-Hydroxypropionate bi-cycle',
              "Propanoyl-CoA metabolism", 'Photorespiration', 'Nucleotide sugar biosynthesis', 'Methylaspartate cycle', 'Malonate semialdehyde pathway', 'Glyoxylate cycle',
              'Glucuronate pathway (uronate pathway)', 'Galactose degradation', 'Ethylmalonyl pathway', 'D-Glucuronate degradation', 'D-galactonate degradation', 'Ascorbate degradation', 'Ascorbate biosynthesis',
              'Trehalose biosynthesis', "Semi-phosphorylative Entner-Doudoroff pathway", 'Pyruvate oxidation', 'PRPP biosynthesis','Pentose phosphate pathway (Pentose phosphate cycle)', 'Pentose phosphate pathway',
              'Glycolysis (Embden-Meyerhof pathway)', 'Glycolysis', 'Gluconeogenesis', 'Entner-Doudoroff pathway', 'Citrate cycle (TCA cycle)', 'Citrate cycle')

cat.heat <- mutate(global.table, taxonomy = fct_relevel(taxonomy, order_tax),
                   Category = fct_relevel(Category, order_cat)) %>%
  ggplot(aes(taxonomy, Category)) +
  geom_tile(aes(fill = log2FoldChange), color = "white") + 
  theme_minimal() + 
  rotate_x_text() + 
  scale_fill_viridis(option="A") +
  scale_x_discrete(labels = label_wrap(15)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ylab("Metabolic Category") 
cat.heat
#ggsave("D:/Project/3_Temperature_effect/gene_expression/exp2/graph/category_heatmap_normalized.png", dpi = 500, width = 3400, height = 3500, units = 'px')

mod.heat <- mutate(global.table, taxonomy = fct_relevel(taxonomy, order_tax),
                   Names = fct_relevel(Names, order_mod)) %>%
  ggplot(aes(taxonomy, Names)) +
  geom_tile(aes(fill = log2FoldChange), color = "white") + 
  theme_minimal() + 
  rotate_x_text() + 
  scale_fill_viridis(option="A") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ylab("Metabolic Module") 
mod.heat
#ggsave("D:/Project/3_Temperature_effect/gene_expression/exp2/graph/module_heatmap_normalized.pdf", dpi = 500, width = 4000, height = 9000, units = 'px')

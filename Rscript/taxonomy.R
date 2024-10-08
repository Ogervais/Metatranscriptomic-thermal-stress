#title: "taxonomy change during thermal stress"
#author: "Ophelie Gervais"
#date: "24/08/23"
  
#### Library --------------------------------------------------------------------
   
library('tidyverse')
library('phyloseq')
library('colorRamp')
library('vegan')
library('scales')
library('ggpubr')
library('Tweedieverse')

#set directory

setwd("D:/Project/3_Temperature_effect/taxonomy")

#create matrix with all the sample----------------------------------------------

RNA15_4 <- read.csv('RNA-15-4_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)
RNA15_6 <- read.csv('RNA-15-6_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)
RNA15_7 <- read.csv('RNA-15-7_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)
RNA15_8 <- read.csv('RNA-15-8_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)
RNA24_3 <- read.csv('RNA-24-3_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)
RNA24_4 <- read.csv('RNA-24-4_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)
RNA24_8 <- read.csv('RNA-24-8_kaiju_summary.tsv', sep = '\t') %>%
  na.omit() %>%
  select(taxon_name, reads)

otu <- full_join(RNA15_4, RNA15_6, by = 'taxon_name') %>%
  full_join(RNA15_7, by = 'taxon_name') %>%
  full_join(RNA15_8, by = 'taxon_name') %>%
  full_join(RNA24_3, by = 'taxon_name') %>%
  full_join(RNA24_4, by = 'taxon_name') %>%
  full_join(RNA24_8, by = 'taxon_name') %>%
  rename(SampleID = taxon_name, 'RNA15-4' = reads.x, 'RNA15-6' = reads.y, 'RNA15-7' = reads.x.x, 'RNA15-8' = reads.y.y, 
         'RNA24-3' = reads.x.x.x, 'RNA24-4' = reads.y.y.y, 'RNA24-8' = reads) %>%
  replace(is.na(.), 0) %>%
  filter(SampleID != 'Viruses') %>%
  filter(!str_detect(SampleID, "Archaea"))

#create OTU table---------------------------------------------------------------

str(otu)

#OTU ID's have to be the rownames
row.names(otu) <- otu$SampleID

#check if the number of reads per sample is correct 
count <- colSums(otu[, c(2:ncol(otu))])

count

#create taxonomy table----------------------------------------------------------

taxa<- otu %>% 
  select(SampleID) %>% 
  separate(SampleID, c("Organism", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           ";")  
taxa<- taxa %>% mutate(Class = ifelse(str_detect(Phylum, "unclassified"), Phylum, Class)) %>%
  mutate(Order = ifelse(str_detect(Class, "unclassified"), Class, Order)) %>%
  mutate(Family = ifelse(str_detect(Order, "unclassified"), Order, Family)) %>%
  mutate(Genus = ifelse(str_detect(Family, "unclassified"), Family, Genus)) %>%
  mutate(Species = ifelse(str_detect(Genus, "unclassified"), Genus, Species)) 

taxa <- select(taxa, -c("Organism")) %>%
  replace(is.na(.), "unclassified Bacteria")

taxa[taxa == ""] <- "unclassified Bacteria"

str(taxa)

#the output is a data frame of characters, and we need taxa to be recognized as factors
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)

str(taxa)

#add the first column from the OTU table so that OTU and taxa tables match
taxa<- cbind(otu$SampleID, taxa)

#rename the first column
colnames(taxa)[1]<- "SampleID"
str(taxa)

#OTU ID's have to be the rownames, same as in the OTU table
rownames(taxa)<- taxa$SampleID

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
otu<- otu %>% 
  select(-SampleID)

taxa<- taxa %>% 
  select(-SampleID)

#Metadata table----------------------------------------------------------------
#contain description of sample names

meta <- read.table('sampleID.txt', sep = '\t', header = TRUE)

str(meta)

#make sample names row names 
rownames(meta)<- meta$SampleID

#and delete the first column because it is now redundant
meta<- meta %>% 
  select(-SampleID)

#make proper levels of the factor so that they are in the wanted order in figures
meta$Temp<- factor(meta$Temp, levels = c("15", "24"))

#Phyloseq-----------------------------------------------------------------------

otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

#transform data to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(meta)

#and put them in one object
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)

#check if everything looks good
sample_sums(phylo_object)       #should sum up the number of all reads per sample
## RNA15-4 RNA15-6 RNA15-7 RNA15-8 RNA24-3 RNA24-4 RNA24-8 
## 208219  166391  410494  208354  111420   71419   52761

sample_names(phylo_object)      #sample names

rank_names(phylo_object)        #taxa levels  

sample_variables(phylo_object)  #factors
## [1] "Temp"

otu_table(phylo_object)[1:3, 1:5]
## OTU Table:          [3 taxa and 5 samples]
##                      taxa are rows

taxa_names(phylo_object)[1:5]

#relativisation
phylo_rel<- transform_sample_counts(phylo_object, function(x) x*100 / sum(x))

#check if everything looks good
sample_sums(phylo_rel)          #should sum up to 100% per sample
## RNA15-4 RNA15-6 RNA15-7 RNA15-8 RNA24-3 RNA24-4 RNA24-8 
##100        100        100        100           100        100        100      

sample_names(phylo_rel)         #sample names
## [1] RNA15-4 RNA15-6 RNA15-7 RNA15-8 RNA24-3 RNA24-4 RNA24-8 

rank_names(phylo_rel)           #taxa levels  
## [1] "Organism" "Domain"   "Phylum"   "Class"    "Order"    "Family"   "Genus"    "Species" 

sample_variables(phylo_rel)     #factors
## [1] "Temp"

otu_table(phylo_rel)[1:3, 1:4]
## OTU Table:          [3 taxa and 4 samples]
##                      taxa are rows

taxa_names(phylo_rel)[1:16]

#to start getting the feel of the data, let's check what had the highest abundance

max(phylo_rel@otu_table)        #the highest relative abundance

rownames(phylo_rel@otu_table)[which.max(apply(phylo_rel@otu_table,MARGIN=1,max))] #row name

#the NAs appear when reads were not assigned down to the respective taxonomy level.  

phylo_rel@otu_table["cellular organisms;Bacteria;FCB group;Bacteroidota/Chlorobiota group;Bacteroidota;Sphingobacteriia;Sphingobacteriales;Sphingobacteriaceae;unclassified Sphingobacteriaceae;Sphingobacteriaceae bacterium;"] #extract that whole row

#Data visualization

#NMDS plot
phylo_rel_nmds<- ordinate(phylo_rel, method = "NMDS", distance = "bray")
vegan::stressplot(phylo_rel_nmds)

#NMDS plot
plot_ordination(phylo_rel, phylo_rel_nmds,
                color = "Temp") +
  theme_bw()

#PCoA
phylo_rel_pcoa<- ordinate(phylo_rel, method = "PCoA", distance = "bray")

#PCoA plot
plot_ordination(phylo_rel, phylo_rel_pcoa,
                color = "Temp") +
  theme_bw()


##means replicate---------------------------------------------------------------
phylo_rel_mean <- merge_samples(phylo_rel, "Temp")

sample_sums(phylo_rel_mean) #sums up to 400!
##  15  24 
## 400 300 

#now divide all relative OTU abundances by the number of replicates
phylo_rel_mean<- transform_sample_counts(phylo_rel_mean, function(x) x*100 / sum(x))

sample_sums(phylo_rel_mean)
##  15  24 
## 100 100

#check taxonomy level names
rank_names(phylo_rel_mean)


## agglomeration on the phylum level
phylo_rel_rel_mean_phylum<- tax_glom(phylo_rel_mean, taxrank = "Phylum")

#check the number of taxa in the whole phyloseq object
phylo_rel_mean 


#check to how many phyla were these assigned to
phylo_rel_rel_mean_phylum 

#check sample sums to make sure nothing was deleted due to NAs (should sum up 100%)
sample_sums(phylo_rel_rel_mean_phylum)


#if we try to color by Phylum to see their names
plot_bar(phylo_rel_rel_mean_phylum, fill = "Phylum")

#transform phyloseq object to a data frame
phylo_rel_rel_mean_phylumDF<- psmelt(phylo_rel_rel_mean_phylum)

str(phylo_rel_rel_mean_phylumDF)

#make the phyla characters, not factors
phylo_rel_rel_mean_phylumDF$Phylum<- as.character(phylo_rel_rel_mean_phylumDF$Phylum)

#add new column with renamed low abundant taxa
phylo_rel_rel_mean_phylumDF<- phylo_rel_rel_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum, Abundance < 1, "< 1%"))

#check all phyla names
unique(phylo_rel_rel_mean_phylumDF$Phylum2)

#if NA on the phylum level, so we will rename them
phylo_rel_rel_mean_phylumDF<- phylo_rel_rel_mean_phylumDF %>% 
  mutate(Phylum2 = replace(Phylum2, Phylum2 == "NA", "unassigned Bacteria"))

#reorder the phyla so that they are stacked according to abundance
phylo_rel_rel_mean_phylumDF$Phylum2<- reorder(phylo_rel_rel_mean_phylumDF$Phylum2,
                                              phylo_rel_rel_mean_phylumDF$Abundance)

#check how many unique phyla are there to find discrete colors for them
unique(phylo_rel_rel_mean_phylumDF$Phylum2)

graphe <- phylo_rel_rel_mean_phylumDF %>%
  ggplot(aes(Sample, Abundance, fill = Phylum2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#333333", 
                               "#e6d800",
                               "#50e991",
                               "#e60049", 
                               "#9b19f5", 
                               "#ffa300", 
                               "#dc0ab4")) +
  labs(y = "Relative abundance [%]",
       fill = "Phylum") +
  theme_bw() +
  theme(
    axis.text = element_text(family = 'Roboto Mono'
    ),
    axis.title.y = element_text(margin = margin(t = 10),
                                size = 15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = 'bold', size = 14
    ),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 13),
    plot.title.position = 'plot',
  ) 

#png("D:/Project/3_Temperature_effect/graph/exp2/taxonomy_phylum_bact_mean.png", res = 500, width = 3800, height = 2400)
graphe
#dev.off()


## agglomeration on the class level
phylo_rel_rel_mean_class<- tax_glom(phylo_rel_mean, taxrank = "Class")

#check to how many class were these assigned to
phylo_rel_rel_mean_class 


#check sample sums to make sure nothing was deleted due to NAs (should sum up 100%)
sample_sums(phylo_rel_rel_mean_class)


#if we try to color by Phylum to see their names, 
# we can see that 135 are too many to see in a stacked barplot
plot_bar(phylo_rel_rel_mean_class, fill = "Class")

#transform phyloseq object to a data frame
phylo_rel_rel_mean_classDF<- psmelt(phylo_rel_rel_mean_class)

str(phylo_rel_rel_mean_classDF)

#make the phyla characters, not factors
phylo_rel_rel_mean_classDF$Class<- as.character(phylo_rel_rel_mean_classDF$Class)

#add new column with renamed low abundant taxa
phylo_rel_rel_mean_classDF<- phylo_rel_rel_mean_classDF %>% 
  mutate(Class2 = replace(Class, Abundance < 2, "< 2%"))

#check all phyla names
unique(phylo_rel_rel_mean_classDF$Class2)


# if NA on the phylum level, so we will rename them
phylo_rel_rel_mean_classDF<- phylo_rel_rel_mean_classDF %>% 
  mutate(Class2 = replace(Class2, Class2 == "NA", "unassigned Bacteria"))

#reorder the phyla so that they are stacked according to abundance
phylo_rel_rel_mean_classDF$Class2<- reorder(phylo_rel_rel_mean_classDF$Class2,
                                            phylo_rel_rel_mean_classDF$Abundance)


#check how many unique phyla are there to find discrete colors for them
unique(phylo_rel_rel_mean_classDF$Class2)

pal <- c("#333333", "#e60049", "#0bb4ff", "#50e991","#9b19f5", "#ffa300", "#dc0ab4", "#00bfa0")

graphe <- phylo_rel_rel_mean_classDF %>%
  ggplot(aes(Sample, Abundance, fill = Class2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal) +
  labs(y = "Relative abundance [%]",
       fill = "Phylum/Class") +
  theme_bw() +
  theme(
    axis.text = element_text(family = 'Roboto Mono'
    ),
    axis.title.y = element_text(margin = margin(t = 10),
                                size = 15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = 'bold', size = 14
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.title.position = 'plot',
  ) 

#png("D:/Project/3_Temperature_effect/graph/exp2/taxonomy_class_bact_mean.png", res = 500, width = 3400, height = 2400)
graphe
#dev.off()


## agglomeration on the order level
phylo_rel_rel_mean_order<- tax_glom(phylo_rel_mean, taxrank = "Order")

#check the number of taxa in the whole phyloseq object
phylo_rel_mean 


#check to how many phyla were these assigned to
phylo_rel_rel_mean_order


#check sample sums to make sure nothing was deleted due to NAs (should sum up 100%)
sample_sums(phylo_rel_rel_mean_order)

#transform phyloseq object to a data frame
phylo_rel_rel_mean_orderDF<- psmelt(phylo_rel_rel_mean_order)

str(phylo_rel_rel_mean_orderDF)

#make the phyla characters, not factors
phylo_rel_rel_mean_orderDF$Order<- as.character(phylo_rel_rel_mean_orderDF$Order)

#add new column with renamed low abundant taxa
phylo_rel_rel_mean_orderDF<- phylo_rel_rel_mean_orderDF %>% 
  mutate(Order2 = replace(Order, Abundance < 2, "< 2%"))

#check all phyla names
unique(phylo_rel_rel_mean_orderDF$Order2)

# if NA on the phylum level, so we will rename them
phylo_rel_rel_mean_orderDF<- phylo_rel_rel_mean_orderDF %>% 
  mutate(Order2 = replace(Order2, Order2 == "NA", "unassigned Bacteria"))

#reorder the phyla so that they are stacked according to abundance
phylo_rel_rel_mean_orderDF$Order2<- reorder(phylo_rel_rel_mean_orderDF$Order2,
                                            phylo_rel_rel_mean_orderDF$Abundance)

#check how many unique phyla are there to find discrete colors for them
unique(phylo_rel_rel_mean_orderDF$Order2)

pal1 <- c('#333333', '#8A17DA', '#0B4EFF', '#3DAB6C', '#fec20c', '#FF5733', '#C40027', '#E649BE', '#B273E3', '#0bb4ff', '#50e991', '#F5E16E')

graphe <- phylo_rel_rel_mean_orderDF %>%
  ggplot(aes(Sample, Abundance, fill = Order2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal1)  +
  labs(y = "Relative abundance [%]",
       fill = "Class/Order") +
  theme_bw() +
  theme(
    axis.text = element_text(family = 'Roboto Mono'
    ),
    axis.title.y = element_text(margin = margin(t = 10),
                                size = 15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = 'bold', size = 14
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.title.position = 'plot',
  ) 

#png("D:/Project/3_Temperature_effect/graph/exp2/taxonomy_order_bact_mean.png", res = 500, width = 3000, height = 2400)
graphe
#dev.off()



### statistical test

#relativisation
phylo_scale<- transform_sample_counts(phylo_object, function(x) x / sum(x))

set.seed(1)

# Calculate bray curtis distance matrix
phylo_bray <- distance(phylo_scale, method = 'bray')

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phylo_object))

## Adonis test
adonis2(phylo_bray ~ Temp, data = sampledf)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 5039
# 
# adonis2(formula = phylo_bray ~ Temp, data = sampledf)
# Df SumOfSqs      R2      F Pr(>F)  
# Temp      1  0.18116 0.54167 5.9091   0.03 *
#   Residual  5  0.15329 0.45833                
# Total     6  0.33446 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##tweedieverse 

#change table to have taxo for genus and not species
otu_tweedie <- full_join(RNA15_4, RNA15_6, by = 'taxon_name') %>%
  full_join(RNA15_7, by = 'taxon_name') %>%
  full_join(RNA15_8, by = 'taxon_name') %>%
  full_join(RNA24_3, by = 'taxon_name') %>%
  full_join(RNA24_4, by = 'taxon_name') %>%
  full_join(RNA24_8, by = 'taxon_name') %>%
  rename(SampleID = taxon_name, 'RNA15-4' = reads.x, 'RNA15-6' = reads.y, 'RNA15-7' = reads.x.x, 'RNA15-8' = reads.y.y, 
         'RNA24-3' = reads.x.x.x, 'RNA24-4' = reads.y.y.y, 'RNA24-8' = reads) %>%
  replace(is.na(.), 0) %>%
  filter(SampleID != 'Viruses') %>%
  filter(!str_detect(SampleID, "Archaea"))

#create OTU table-
#check if the number of reads per sample is correct 
count <- colSums(otu_tweedie[, c(2:ncol(otu_tweedie))])
otu_tweedie_gen <- otu_tweedie %>%
  separate(SampleID, c("Organism", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),";") %>%
  select(!Species) %>%
  unite(SampleID, c("Organism", "Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  group_by(SampleID) %>%
  summarise(across(starts_with('RNA'), sum)) %>%
  column_to_rownames(var = 'SampleID')



tweedi_test <- Tweedieverse(otu_tweedie_gen,
                            meta,
                            output = './tweedie',
                            fixed_effects = 'Temp',
                            base_model = 'CPLM'
)

tweedi_test <- tweedi_test %>%
  separate(feature, c("Organism", "Domain", "Phylum", "Class", "Order", "Family", "Genus"),";")
tweedi_test <- tweedi_test %>%
  filter(Genus != 'environmental samples')

palette_loli <- c("#C40027", "#9b19f5", "#F5E16E", "#0bb4ff", "#50e991")
ggdotchart(tweedi_test, 'Genus', 'coef',
           color = 'Order',
           palette = palette_loli,
           sorting = 'descending',
           add = 'segments',
           rotate = TRUE,
           dot.size = 4,
           ylab = FALSE) +
  guides(color = guide_legend(nrow = 3, 
                              title.position = 'top')) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.position = 'top',  
        legend.justification = 'right') +
  xlab('coefficient') 

#ggsave("./tweedie/lolipopchart.png", dpi = 500, width = 2600, height = 3200, units = 'px')

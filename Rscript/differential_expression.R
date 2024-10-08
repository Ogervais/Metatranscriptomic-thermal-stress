#title: "Differential gene expression"
#author: "Ophelie Gervais"
#date: "31/01/24"

  
#Test Deseq2 after previous normalisation taxa

#### Library --------------------------------------------------------------------------------

library('tidyverse')
library('DESeq2')
library("RColorBrewer")
library("pheatmap")
library("ggplotify")
library('gridExtra')

#used df_filtered from the normalization.R script
source("D:/Project/3_Temperature_effect/R_script/Github_script/Normalization.R")

# Convert to matrix
countdata <- as.matrix(df_filtered)
countdata <- matrix(as.integer(round(countdata)), nrow = nrow(countdata),
                    dimnames = dimnames(countdata))
#or storage.mode(countdata) <- 'integer'
head(df_filtered)

# Assign condition 
metadata <- read.table("SraRunTable.csv", header = TRUE, sep = ";", row.names = 1)
metadata <- metadata[-5,] #remove sample 24.2
metadata$Temperature <- factor(metadata$Temperature)
metadata$Condition <- factor(metadata$Condition)

all(rownames(metadata) %in% colnames(countdata))

#########DESEQ2#############
data.col = ncol(countdata)
data.row = nrow(countdata)
data_all <- DESeqDataSetFromMatrix(countData = countdata,
                                   colData = metadata,
                                   design = ~Temperature)

#PCA
vsd <- varianceStabilizingTransformation(data_all)

PCA <- plotPCA(vsd, intgroup=c("Temperature"), returnData = TRUE)
percentVar <- round(100 * attr(PCA, "percentVar"))

gplot_PCA <- PCA %>%
  mutate(Temperature = fct_relevel(Temperature, "15", '24')) %>%
  ggplot(aes(PC1, PC2, color = Temperature)) +
  geom_point(size = 4, alpha = 1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold"), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size =14)) +
  scale_color_manual(values = c("#8d99ae", "#8c0d07")) 
gplot_PCA
#ggsave("graph/pca_vsd.png", dpi = 500, width=2800, height=2800, units = 'px')

##Heatmaps distance
##Clustering
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
color = colorRampPalette(brewer.pal(11,"Spectral"))(100)
df <- as.data.frame(colData(vsd)[,-3])
subset_df <- df[, "Temperature", drop = FALSE]
ann_colors = list(
  Temperature = c("15" = "#8d99ae", "24" = "#8c0d07"))

pheatmap_ <- as.grob(pheatmap(sampleDistMatrix, 
                              clustering_distance_rows = sampleDists, 
                              clustering_distance_cols = sampleDists, 
                              col = color,
                              fontsize = 11,
                              fontsize_col = 0.1,
                              annotation_col = subset_df, 
                              annotation_colors = ann_colors))
my_plot <- grid.arrange(pheatmap_, nrow=1, ncol=1)
#ggsave("graph/clustering_vsd.png", my_plot, dpi = 500, units = 'px', width = 3000, height = 2000)

#Differential tests
data_all$Temperature <- relevel(data_all$Temperature, ref = "15")

#Stop DESeq2 from performing additional normalization
normFactors <- matrix(1, ncol = data.col, nrow = data.row)
normalizationFactors(data_all) <- normFactors
dds <- DESeq(data_all, quiet = TRUE)

axe <- read.delim("eggnog_mmseqs.emapper.annotations") 
axe <- axe[,1]
dds2 <- dds[rownames(dds) %in% axe,] #kept only gene who are annotated

res <- results(dds2)
res

#Result extraction
res_CTL_vs_24 <- results(dds2, contrast=c("Temperature","24","15") , test="Wald") #ref is 15

annotation <- read.delim("eggnog_mmseqs.emapper.annotations", row.names=1) %>%
  separate(max_annot_lvl, into = c("number", "taxonomy"), sep = "\\|") %>%
  mutate_all(~ ifelse(. == '-', NA, .)) 


res_CTL_vs_24_annot <- merge(DataFrame(res_CTL_vs_24), annotation, all.x = TRUE, by = "row.names")

#write.table(as.data.frame(res_CTL_vs_24_annot),file="result/deseq2/15_vs_24_normalized_allDE.txt", sep="\t", row.names = FALSE)

#Filter statisticaly significant genes
res_CTL_vs_24_annot_pval <- subset(res_CTL_vs_24_annot[order(res_CTL_vs_24_annot$padj),], padj < 0.05)

#write.table(as.data.frame(res_CTL_vs_24_annot_pval),file="result/deseq2/15_vs_24_normalized_DE.txt", sep="\t", row.names = FALSE)

########################################################
# Cardiac Differentiation RNA-seq Analysis
# Dataset: GSE116574
# Author: Shameem Fathima
########################################################

############################
# Load Libraries
############################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)

############################
# Create Folders
############################

dir.create("Results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

############################
# Load Data
############################

counts <- read.csv(
  "Raw_Counts/GSE116574_count_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "Metadata/metadata.csv",
  row.names = 1
)

############################
# Select only 3D samples
############################

metadata <- metadata[
  metadata$condition %in% c("3D4","3D5","3D8"),
]

counts <- counts[
  ,
  rownames(metadata)
]

metadata$condition <- factor(
  metadata$condition,
  levels = c("3D4","3D5","3D8")
)

############################
# DESeq2 Analysis
############################

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

############################
# Differential Expression
############################

res_3D5_vs_3D4 <- results(
  dds,
  contrast = c("condition","3D5","3D4")
)

res_3D8_vs_3D5 <- results(
  dds,
  contrast = c("condition","3D8","3D5")
)

res_3D8_vs_3D4 <- results(
  dds,
  contrast = c("condition","3D8","3D4")
)

write.csv(
  as.data.frame(res_3D5_vs_3D4),
  "Results/DEGs_3D5_vs_3D4.csv"
)

write.csv(
  as.data.frame(res_3D8_vs_3D5),
  "Results/DEGs_3D8_vs_3D5.csv"
)

write.csv(
  as.data.frame(res_3D8_vs_3D4),
  "Results/DEGs_3D8_vs_3D4.csv"
)

############################
# Variance Stabilization
############################

vsd <- vst(dds)

############################
# PCA Plot
############################

pcaData <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

p <- ggplot(
  pcaData,
  aes(PC1, PC2, color = condition)
) +
  geom_point(size = 4) +
  xlab(
    paste0("PC1: ", percentVar[1], "%")
  ) +
  ylab(
    paste0("PC2: ", percentVar[2], "%")
  ) +
  theme_bw()

ggsave(
  "Figures/PCA_plot.png",
  p,
  width = 8,
  height = 6
)

############################
# Sample Correlation Heatmap
############################

sampleDists <- dist(
  t(assay(vsd))
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

png(
  "Figures/Sample_correlation_heatmap.png",
  width = 1200,
  height = 1000
)

pheatmap(sampleDistMatrix)

dev.off()

############################
# Significant DEGs
############################

res <- as.data.frame(
  res_3D8_vs_3D4
)

sig <- subset(
  res,
  padj < 0.05
)

write.csv(
  sig,
  "Results/Significant_DEGs.csv"
)

############################
# Volcano Plot
############################

png(
  "Figures/volcano_plot.png",
  width = 1400,
  height = 1000
)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1
)

dev.off()

############################
# Top 20 DEG Heatmap
############################

top20 <- head(
  rownames(sig),
  20
)

png(
  "Figures/Heatmap_Top20_DEGs.png",
  width = 1200,
  height = 1000
)

pheatmap(
  assay(vsd)[top20, ],
  scale = "row"
)

dev.off()

############################
# Top 50 DEG Heatmap
############################

top50 <- head(
  rownames(sig),
  50
)

png(
  "Figures/Heatmap_Top50_DEGs.png",
  width = 1200,
  height = 1000
)

pheatmap(
  assay(vsd)[top50, ],
  scale = "row"
)

dev.off()

############################
# GO Enrichment
############################

sig_genes <- rownames(sig)

gene_df <- bitr(
  sig_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

write.csv(
  as.data.frame(ego),
  "Results/GO_Results.csv"
)

png(
  "Figures/GO_enrichment_barplot.png",
  width = 1200,
  height = 900
)

barplot(
  ego,
  showCategory = 15
)

dev.off()

png(
  "Figures/GO_enrichment_dotplot.png",
  width = 1200,
  height = 900
)

dotplot(
  ego,
  showCategory = 10
)

dev.off()

############################
# KEGG Enrichment
############################

ekegg <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa"
)

write.csv(
  as.data.frame(ekegg),
  "Results/KEGG_Results.csv"
)

png(
  "Figures/KEGG_barplot.png",
  width = 1200,
  height = 900
)

barplot(
  ekegg,
  showCategory = 15
)

dev.off()

png(
  "Figures/KEGG_dotplot.png",
  width = 1200,
  height = 900
)

dotplot(
  ekegg,
  showCategory = 15
)

dev.off()

png(
  "Figures/KEGG_Ridgeplot.png",
  width = 1200,
  height = 900
)

ridgeplot(
  ego
)

dev.off()

############################
# Save Session Info
############################

capture.output(
  sessionInfo(),
  file = "Results/sessionInfo.txt"
)

########################################################
# END OF ANALYSIS
########################################################

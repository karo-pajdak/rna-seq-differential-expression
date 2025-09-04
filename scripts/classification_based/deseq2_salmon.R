# =======================================================================================
# DESeq2 Script (Salmon Version)
# =======================================================================================
# This script ...
#
# =======================================================================================

# Install necessary packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")
BiocManager::install("tximport")
BiocManager::install("txdbmaker")
BiocManager::install("apeglm")
BiocManager::install("ggrepel")

library(tximport)
library(txdbmaker)
library(DESeq2)
library(ggrepel)

## Set Up ##
# Define working directory paths
input_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/results/salmon"
meta_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/config/"
annotation_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/annotation"
output_dir <- file.path(Sys.getenv("HOME"), "BioinformaticsPortfolio/rna-seq-differential-expression/results/deseq")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


## Import and Gene-level Summary ##
setwd(meta_dir)
# List Salmon output files
metadata <- read.csv("metadata.csv", row.names = 1, stringsAsFactors = FALSE)
files <- file.path(input_dir, rownames(metadata), "quant.sf")
names(files) <- rownames(metadata)

# Create transcript-to-gene mapping (tx2gene) from gtf file
setwd(annotation_dir)
txdb <- txdbmaker::makeTxDbFromGFF("annotation.gtf", format = "gtf")
tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")

#Import
txi <- tximport(files, type="salmon", tx2gene=tx2gene)


## Differential Expression ## 
# Build DESeq dataset and pre-filter lowly expressed genes
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)

# Shrink log2 fold changes
res <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")


## QC & Model Diagnostics ##
# Dispersion plot - check model fit
png("dispersion_plot_salmon.png", width=800, height=600)
disp_plot <- plotDispEsts(dds, main="Dispersion Plot")
dev.off()

# PCA Plot - check sample clustering by condition
png("PCA_plot_salmon.png", width=800, height=600)
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# Sample-to-sample distance heatmap - confirm replicates cluster together
# Extract matrix of counts, compute Euclidean distance, convert to matrix, plot
png("sample_distance_salmon.png", width=800, height=600)
vsd_mat <- assay(vsd)
sample_dists <- dist(t(vsd_mat))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste0(colnames(vst(dds)), " (", vst(dds)$condition, ")")
colnames(sample_dist_matrix) <- paste0(colnames(vst(dds)), " (", vst(dds)$condition, ")")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         main = "Sample-to-Sample Distance Heatmap")
dev.off()

## Differential Expression Results ##
# MA Plot (Raw & Shrunk)
png("MA_plot_raw_salmon.png", width=800, height=600)
plotMA(dds, ylim = c(-5, 5), main="MA Plot - Raw")
dev.off()

png("MA_plot_shrunk_salmon.png", width=800, height=600)
plotMA(res, ylim = c(-5, 5), main="MA Plot - LFC Shrunk")
dev.off()

# Volcano Plot - highlights significant genes
png("volcano_plot_salmon.png", width=800, height=600)
top_genes <- head(res[order(res$padj), ], 5)
# Padding so gene names are not cut off
x_pad <- 0.1 * (max(res$log2FoldChange, na.rm=TRUE) - min(res$log2FoldChange, na.rm=TRUE))
y_pad <- 0.1 * (max(-log10(res$padj), na.rm=TRUE) - 0)
with(res, plot(log2FoldChange, -log10(padj), 
               pch=20, main="Volcano Plot", 
               xlim=c(min(log2FoldChange, na.rm=TRUE)-x_pad, max(log2FoldChange, na.rm=TRUE)+x_pad),
               ylim=c(0, max(-log10(padj), na.rm=TRUE)+y_pad)))
# Red points for significant DEGs
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 1), 
     points(log2FoldChange, -log10(padj), pch=20, col="red"))
# Blue dashed threshold lines
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1,1), col="blue", lty=2)
# Label top genes
text(top_genes$log2FoldChange, -log10(top_genes$padj),
     labels=rownames(top_genes), pos=3, cex=0.6, col="black", xpd=TRUE)
dev.off()

# Counts Summary

# Results Summary


## Gene-Level Exploration ##
# Heatmap of top DEGs

# Normalized counts plot


## Biological Interpretation ##

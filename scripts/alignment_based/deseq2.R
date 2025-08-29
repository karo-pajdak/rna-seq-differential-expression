# =======================================================================================
# DESeq2 Script
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
BiocManager::install("gglot2")
BiocManager::install("pheatmap")

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Define working directory paths
input_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/results/counts/"
meta_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/config/"
output_dir <- file.path(Sys.getenv("HOME"), "BioinformaticsPortfolio/rna-seq-differential-expression/results/deseq")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set working directory to input directory
setwd(input_dir)

# Read the counts from featureCounts
countdata <- read.table("alignment_counts.txt", row.names = 1, header = TRUE)

# Clean up the data to only include gene name and counts per sample; sample names must match metadata
countdata <- countdata[,c(6,7,8,9)]
colnames(countdata) <- gsub(".*(SRR[0-9]+).*", "\\1", colnames(countdata))

# Set working directory to meta data directory
setwd(meta_dir)

# Read meta data with sample conditions
metadata <- read.csv("metadata.csv", row.names = 1, header = TRUE)

# Confirm metadata row names match countdata column names
if (!all(colnames(countdata) %in% rownames(metadata))) {
  stop("ERROR: Some column names in countdata are missing in metadata!")
}
if (!all(rownames(metadata) %in% colnames(countdata))) {
  stop("ERROR: Some row names in metadata are missing in countdata!")
}

# Filter by total reads per gene (at least 10 total reads across all samples)
keep_genes <- rowSums(countdata) >= 10
filtered_counts <- countdata[keep_genes, ]

# Create the DESeq2 object
de <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = ~condition)

# Run DESeq2
de_results <- DESeq(de)

# Extract results
res <- results(de_results)

# Print results to csv
write.csv(res, file.path(output_dir, "deseq_results.csv"))

# Transform data for visualization (using rlog() for smaller data set)
rlog_counts <- rlog(de_results, blind=TRUE)

# Create diagnostic plots
# PCA Plot - see sample clustering
# Enhance PCA plot with ggplot
png("PCA_plot.png", width=800, height=600)
pcaData <- plotPCA(rlog_counts, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# Dispersion plot - check model fit
png("dispersion_plot.png", width=800, height=600)
disp_plot <- plotDispEsts(de_results, main="Dispersion Plot")
dev.off()

# MA plot - fold changes vs expression
png("MA_plot.png", width=800, height=600)
plotMA(rlog_counts, ylim = c(-5, 5), main="MA Plot")
dev.off()

# Sample distance heatmap - check sample clustering
sample_dists <- dist(t(assay(rlog_counts)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste0(colnames(rlog_counts), " (", rlog_counts$condition, ")")
colnames(sample_dist_matrix) <- paste0(colnames(rlog_counts), " (", rlog_counts$condition, ")")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         main="Sample Distance Heatmap")

# Filter by significant genes using padj and log2FoldChange thresholds
padj_threshold <- 0.05
lfc_threshold <- 1    # 2x fold significance
sig_genes <- subset(res, !is.na(padj) & padj < padj_threshold & abs(log2FoldChange) > lfc_threshold)

# Count significant genes
cat("Total genes tested:", nrow(res), "\n")
cat("Significant genes (padj < 0.05, |LFC| > 1):", nrow(sig_genes), "\n")
cat("Upregulated genes:", sum(sig_genes$log2FoldChange > lfc_threshold, na.rm = TRUE), "\n")
cat("Downregulated genes:", sum(sig_genes$log2FoldChange < -lfc_threshold, na.rm = TRUE), "\n")

# Order significant and all results by adjusted p-value
sig_genes <- sig_genes[order(sig_genes$padj),]
res_ordered <- res[order(res$padj),]


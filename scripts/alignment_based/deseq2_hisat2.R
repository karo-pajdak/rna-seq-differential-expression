# =======================================================================================
# DESeq2 Script (HISAT2 Version)
# =======================================================================================
# This script performs differential expression analysis on gene-level counts generated
# from HISAT2 + featureCounts. It includes QC diagnostics, shrinkage of log2 fold 
# changes, visualization (PCA, volcano, MA, heatmaps, normalized counts), 
# and optional biological interpretation (GO/KEGG enrichment).
# =======================================================================================

# Install necessary packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
BiocManager::install("apeglm")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("dplyr")
BiocManager::install("tidyr")
BiocManager::install("ggrepel")

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggrepel)


## Set Up ##
# Define working directory paths
input_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/results/counts/"
meta_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/config/"
output_dir <- file.path(Sys.getenv("HOME"), "BioinformaticsPortfolio/rna-seq-differential-expression/results/deseq_hisat2")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set working directory to input directory
setwd(input_dir)


## Import and Gene-level Summary ##
# Read the counts from featureCounts
countdata <- read.table("alignment_counts.txt", row.names = 1, header = TRUE)

# Clean up the data to only include gene name and counts per sample; sample names must match metadata
countdata <- countdata[,c(6,7,8,9)]
colnames(countdata) <- gsub(".*(SRR[0-9]+).*", "\\1", colnames(countdata))

# Set working directory to meta data directory and read
setwd(meta_dir)
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


## Differential Expression ## 
# Create the DESeq2 object and run DESeq2
de <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = ~condition)
de_results <- DESeq(de)

# Shrink log2 fold changes for stable visualization
res <- lfcShrink(de_results, coef="condition_treated_vs_control", type="apeglm")


## QC & Model Diagnostics ##
# Dispersion plot - check model fit
png("dispersion_plot_hisat2.png", width=800, height=600)
plotDispEsts(dds, main="Dispersion Plot")
dev.off()

# VST transformation for PCA and heatmap
vsd <- vst(de_results, blind=TRUE)

# PCA Plot - check sample clustering by condition
png("PCA_plot_hisat2.png", width=800, height=600)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# Sample-to-sample distance heatmap - confirm replicates cluster together
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste0(colnames(vsd), " (", vsd$condition, ")")
colnames(sample_dist_matrix) <- paste0(colnames(vsd), " (", vsd$condition, ")")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png("sample_distance_heatmap.png", width=800, height=600)
pheatmap(sample_dist_matrix,
         clustering_distance_rows=sample_dists,
         clustering_distance_cols=sample_dists,
         col=colors,
         main="Sample-to-Sample Distance Heatmap")
dev.off()


## Differential Expression Results ##
# MA Plot (Raw & Shrunk)
png("MA_plot_raw_hisat2.png", width=800, height=600)
plotMA(de_results, ylim = c(-5, 5), main="MA Plot - Raw")
dev.off()

png("MA_plot_shrunk_hisat2.png", width=800, height=600)
plotMA(res, ylim = c(-5, 5), main="MA Plot - LFC Shrunk")
dev.off()

# Volcano Plot - highlights significant genes
png("volcano_plot_hisat2.png", width=800, height=600)
top_genes <- head(results_full[order(results_full$padj),], 5)
x_pad <- 0.1 * (max(res$log2FoldChange, na.rm=TRUE) - min(res$log2FoldChange, na.rm=TRUE))
y_pad <- 0.1 * (max(-log10(res$padj), na.rm=TRUE) - 0)
with(res, plot(log2FoldChange, -log10(padj),
               pch=20, main="Volcano Plot",
               xlim=c(min(log2FoldChange, na.rm=TRUE)-x_pad, max(log2FoldChange, na.rm=TRUE)+x_pad),
               ylim=c(0, max(-log10(padj), na.rm=TRUE)+y_pad)))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 1),
     points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1,1), col="blue", lty=2)
text(top_genes$log2FoldChange, -log10(top_genes$padj),
     labels=top_genes$gene_id, pos=3, cex=0.6, col="black", xpd=TRUE)
dev.off()

# Counts Summary (significant genes padj < 0.05 and upregulated vs downregulated)
res_filtered <- res[!is.na(res$padj), ]
total_genes <- nrow(res_filtered)
sig_genes <- sum(res_filtered$padj < 0.05)
up_reg <- sum(res_filtered$padj < 0.05 & res_filtered$log2FoldChange > 0)
down_reg <- sum(res_filtered$padj < 0.05 & res_filtered$log2FoldChange < 0)

sink("counts_summary_deseq_hisat2.txt")
cat("Counts Summary - DESeq (HISAT2)\n")
cat("Total genes tested: ", total_genes, "\n")
cat("Number of significant genes (padj < 0.05): ", sig_genes, "\n")
cat("Number of significant upregulated genes: ", up_reg, "\n")
cat("Number of significant downregulated genes: ", down_reg, "\n")
sink()

# Export full results, significant genes, top 10 DEGs
results_full <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id") %>% arrange(padj)
results_sig <- results_full %>% filter(padj < 0.05)
results_top10 <- results_full %>% arrange(desc(abs(log2FoldChange))) %>% slice_head(n=10)

write.csv(results_full, "results_full_deseq_hisat2.csv")
write.csv(results_sig, "results_sig_deseq_hisat2.csv")
write.csv(results_top10, "results_top10_deseq_hisat2.csv")


## Gene-Level Exploration ##
# Heatmap of top 50 DEGs
top50 <- results_sig %>% slice_head(n=50)
png("top50_heatmap_hisat2.png", width=1000, height=800)
pheatmap(assay(vsd)[top50$gene_id,],
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=metadata["condition"],
         scale="row",
         fontsize_row=10,
         fontsize_col=10,
         color=colorRampPalette(c("navy","white","firebrick3"))(50),
         main="Top 50 Differentially Expressed Genes")
dev.off()

# Normalized counts plots for top 3 genes
top3_genes <- results_top10$gene_id[1:3]
norm_counts <- counts(dds, normalized=TRUE)
top3_counts <- norm_counts[top3_genes,]

df_top3 <- as.data.frame(top3_counts) %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to="sample", values_to="count") %>%
  left_join(metadata %>% tibble::rownames_to_column("sample"), by="sample")

png("norm_counts_top3_hisat2.png", width=1000, height=800)
ggplot(df_top3, aes(x=condition, y=count, fill=condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width=0.1, size=2, alpha=0.9) +
  facet_wrap(~gene_id, scales="free_y") +
  theme_minimal(base_size=12) +
  labs(title="Normalized Counts for Top 3 DEGs", x="Condition", y="Normalized Counts") +
  theme(legend.position="none")
dev.off()


## Biological Interpretation ##
sig_genes_vector <- results_sig$gene_id
sig_genes_clean <- gsub("\\..*$", "", sig_genes_vector)

entrez <- mapIds(org.Hs.eg.db, keys=sig_genes_clean,
                 column="ENTREZID", keytype="ENSEMBL", multiVals="first")

ego <- enrichGO(gene=entrez,
                OrgDb=org.Hs.eg.db,
                keyType="ENTREZID",
                ont="BP",
                pAdjustMethod="BH",
                qvalueCutoff=0.05)

ekegg <- enrichKEGG(gene=entrez,
                    organism="hsa",
                    pvalueCutoff=0.05)

png("go_enrichment_hisat2.png", width=800, height=600)
barplot(ego, showCategory=15, title="Top GO Terms")
dev.off()

png("kegg_enrichment_hisat2.png", width=800, height=600)
dotplot(ekegg, showCategory=15, title="KEGG Pathways")
dev.off()
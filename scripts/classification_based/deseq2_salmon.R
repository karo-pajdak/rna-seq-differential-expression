# =======================================================================================
# DESeq2 Script (Salmon Version)
# =======================================================================================
# This script performs RNA-seq differential expression analysis using Salmon quantifications. 
# It includes data import, gene-level summarization, DE testing with fold-change shrinkage, 
# quality control (dispersion, PCA, sample distances), and visualization of DE results 
# (MA plots, volcano plots, heatmaps, normalized counts). 
# Functional enrichment (GO and KEGG) is performed on significant genes, and all results 
# are exported for reproducibility and portfolio presentation.
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
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(tximport)
library(txdbmaker)
library(DESeq2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db) 


## Set Up ##
# Define working directory paths
input_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/results/salmon"
meta_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/config/"
annotation_dir <- "~/BioinformaticsPortfolio/rna-seq-differential-expression/data/reference/annotation"
output_dir <- file.path(Sys.getenv("HOME"), "BioinformaticsPortfolio/rna-seq-differential-expression/results/deseq_salmon")

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
tx_table <- transcripts(txdb, columns = c("tx_name", "gene_id"))
tx2gene <- tx_table %>%
  as.data.frame() %>%
  dplyr::select(TXNAME = tx_name, GENEID = gene_id) %>%
  mutate(GENEID = sapply(GENEID, `[`, 1))  # extract first element from each list

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

# VST transformation for PCA and heatmap
vsd <- vst(dds, blind = TRUE)

# PCA Plot - check sample clustering by condition
png("PCA_plot_salmon.png", width=800, height=600)
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

# Counts Summary (significant genes padj < 0.05 and upregulated vs downregulated)
res_raw <- results(dds)
res_filtered <- res_raw[!is.na(res_raw$padj),]
total_genes <- nrow(res_filtered)
sig_genes <- sum(res_filtered$padj < 0.05)
up_reg <- sum(res_filtered$padj < 0.05 & res_filtered$log2FoldChange > 0)
down_reg <- sum(res_filtered$padj < 0.05 & res_filtered$log2FoldChange < 0)

setwd(output_dir)
sink("counts_summary_deseq_salmon.txt")
cat("Counts Summary - DESeq (Salmon)\n")
cat("Total genes tested: ", total_genes, "\n")
cat("Number of significant genes (padj < 0.05): ", sig_genes, "\n")
cat("Number of significant upregulated genes: ", up_reg, "\n")
cat("Number of significant downregulated genes: ", down_reg, "\n")
sink()

# Export full results, significant genes, top 10 DEGs
results_full <- as.data.frame(res_raw) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  arrange(padj)

results_sig <- as.data.frame(res_raw) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  filter(padj < 0.05) %>%
  arrange(padj)

results_top_genes <- as.data.frame(res) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 10)

write.csv(results_full, file = "results_full_deseq_salmon.csv")
write.csv(results_sig, file = "results_sig_deseq_salmon.csv")
write.csv(results_top_genes, file = "results_top10_deseq_salmon.csv")


## Gene-Level Exploration ##
# Heatmap of top DEGs
top50 <- as.data.frame(res_raw) %>%
  tibble::rownames_to_column(var = "gene_id") %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50)

png("top_genes_heatmap_salmon.png", width=1000, height=800)
pheatmap(assay(vsd)[top50$gene_id,], 
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=metadata["condition"],
         scale="row", 
         fontsize_row = 10,
         fontsize_col = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main="Top 50 Differentially Expressed Genes")
dev.off()

# Normalized counts plot for top 3 genes
top3_genes <- head(results_top_genes$gene_id,3)
norm_counts <-counts(dds, normalized=TRUE)
top3_counts <- norm_counts[top3_genes,]
# Reshape data to long format for ggplot2
metadata2 <- metadata %>% tibble::rownames_to_column(var = "samples")
df_top3 <- as.data.frame(top3_counts) %>%
  tibble::rownames_to_column(var="gene_id") %>%
  pivot_longer(-gene_id, names_to = "samples", values_to = "count")
df_top3 <- df_top3 %>% left_join(metadata2, by = "samples")

png("normalized_counts_top3_salmon.png", width=1000, height=800)
ggplot(df_top3, aes(x = condition, y = count, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.9) +
  facet_wrap(~gene_id, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Normalized Counts for Top 3 DEGs",
       x = "Condition",
       y = "Normalized Counts") +
  theme(legend.position = "none")
dev.off()

## Biological Interpretation ##
sig_genes_vector <- as.data.frame(res_filtered) %>% 
  filter(padj < 0.05) %>%
  rownames()
sig_genes_clean <- gsub("\\..*$", "", sig_genes_vector)

entrez <- mapIds(org.Hs.eg.db, keys = sig_genes_clean,
                 column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

ego <- enrichGO(gene = entrez,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

ekegg <- enrichKEGG(gene = entrez,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)

png("go_enrichment_salmon.png", width=800, height=600)
barplot(ego, showCategory = 15, title = "Top GO Terms")
dev.off()

png("kegg_enrichment_salmon.png", width=800, height=600)
dotplot(ekegg, showCategory = 15, title = "KEGG Pathways")
dev.off()

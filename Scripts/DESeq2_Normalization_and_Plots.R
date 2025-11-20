# -----------------------------------------
# Load libraries
# -----------------------------------------
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(org.Dr.eg.db)  
library(biomaRt)

# -----------------------------------------
# Load Data
# -----------------------------------------
counts <- read.csv("C:/Users/Dell/Downloads/gene_counts_matrix.csv", row.names = 1)
meta   <- read.csv("C:/Users/Dell/Downloads/metadata (1).csv", row.names = 1)

# Add condition column
meta$condition <- factor(ifelse(meta$tissue == "Tumor", "Tumor", "Normal"))

# Ensure sample names match
meta$Run <- rownames(meta)
rownames(meta) <- meta$Run
meta <- meta[colnames(counts), ]

# -----------------------------------------
# DESeq2 Dataset
# -----------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta,
                              design    = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Set reference level
dds$condition <- relevel(dds$condition, ref = "Normal")

# -----------------------------------------
# Run DESeq2
# -----------------------------------------
dds <- DESeq(dds, fitType = "mean")

# -----------------------------------------
# Extract Results
# -----------------------------------------
res <- results(dds, alpha = 0.05)  # padj < 0.05
res <- res[order(res$pvalue), ] # order by p-value

# Subset significant DEGs |log2FC| > 1
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# -----------------------------------------
# Convert Ensembl IDs to Gene Symbols
# -----------------------------------------
# Using biomaRt
ensembl <- useEnsembl(biomart="genes", dataset="drerio_gene_ensembl")  # zebrafish
gene_map <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = rownames(res),
                  mart = ensembl)

# Map symbols to results
res$gene_name <- gene_map$external_gene_name[match(rownames(res), gene_map$ensembl_gene_id)]
sig_res$gene_name <- gene_map$external_gene_name[match(rownames(sig_res), gene_map$ensembl_gene_id)]

# Replace NA with original ID if no mapping
res$gene_name[is.na(res$gene_name)] <- rownames(res)[is.na(res$gene_name)]
sig_res$gene_name[is.na(sig_res$gene_name)] <- rownames(sig_res)[is.na(sig_res$gene_name)]

# -----------------------------------------
# Save Results
# -----------------------------------------
write.csv(as.data.frame(res), "DEG_results_with_gene_names.csv")
write.csv(as.data.frame(sig_res), "DEG_significant_log2FC1_with_gene_names.csv")

# -----------------------------------------
# Normalized Counts
# -----------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
rownames(norm_counts) <- gene_map$external_gene_name[match(rownames(norm_counts), gene_map$ensembl_gene_id)]
write.csv(norm_counts, "normalized_counts_with_gene_names.csv")

# -----------------------------------------
# PCA Plot
# -----------------------------------------
rld <- vst(dds)
pdf("PCA_plot.pdf")
plotPCA(rld, intgroup = "condition")
dev.off()

png("PCA_plot.png", width = 3000, height = 3000, res = 300)
plotPCA(rld, intgroup = "condition")
dev.off()

# -----------------------------------------
# Volcano Plot
# -----------------------------------------
data <- as.data.frame(res)
data$Expression <- case_when(
  data$log2FoldChange >= 1 & data$padj <= 0.05 ~ "Up-regulated",
  data$log2FoldChange <= -1 & data$padj <= 0.05 ~ "Down-regulated",
  TRUE ~ "Unchanged"
)

# Top 10 significant genes
Top_Hits <- head(arrange(data, padj), 10)
data$label <- if_else(data$gene_name %in% Top_Hits$gene_name, data$gene_name, "")

pdf("volcano_plot.pdf")
ggplot(data, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color = Expression), size = 1.5) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  geom_text_repel(aes(label = label), size = 3) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  theme_minimal()
dev.off()

png("volcano_plot.png", width = 3000, height = 3000, res = 300)
ggplot(data, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color = Expression), size = 1.5) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  geom_text_repel(aes(label = label), size = 3) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  theme_minimal()
dev.off()

# -----------------------------------------
# MA Plot
# -----------------------------------------
pdf("MA_plot.pdf")
plotMA(res, ylim = c(-5,5))
dev.off()

png("MA_plot.png", width = 3000, height = 3000, res = 300)
plotMA(res, ylim = c(-5,5))
dev.off()

# -----------------------------------------
# Heatmap of top 20 DEGs
# -----------------------------------------
top20_genes <- head(sig_res[order(sig_res$padj), "gene_name"], 20)
heat_data <- norm_counts[top20_genes, ]
annotation_col <- data.frame(Condition = meta$condition)
rownames(annotation_col) <- colnames(heat_data)

pdf("heatmap_top20_DEGs.pdf", width = 8, height = 6)
pheatmap(heat_data, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col, scale = "row", fontsize_row = 8)
dev.off()

png("heatmap_top20_DEGs.png", width = 3000, height = 3000, res = 300)
pheatmap(heat_data, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col, scale = "row", fontsize_row = 8)
dev.off()


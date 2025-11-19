# -----------------------------------------
# Load libraries
# -----------------------------------------
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)

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
res <- res[order(res$pvalue), ]    # order by p-value

# Subset significant DEGs |log2FC| > 1
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Save results
write.csv(as.data.frame(res), "DEG_results.csv")
write.csv(as.data.frame(sig_res), "DEG_significant_log2FC1.csv")

# -----------------------------------------
# Normalized Counts
# -----------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "normalized_counts.csv")

# -----------------------------------------
# PCA Plot
# -----------------------------------------
rld <- vst(dds)
pdf("PCA_plot.pdf")  # vector PDF
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
data$label <- if_else(rownames(data) %in% rownames(Top_Hits), rownames(data), "")

# Volcano Plot PDF
pdf("volcano_plot.pdf")
ggplot(data, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color = Expression), size = 1.5) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  geom_text_repel(aes(label = label), size = 3) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  theme_minimal()
dev.off()

# Volcano Plot PNG
png("volcano_plot.png", width = 3000, height = 3000, res = 300)
ggplot(data, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color = Expression), size = 1.5) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  geom_text_repel(aes(label = label), size = 3) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  theme_minimal()
dev.off()


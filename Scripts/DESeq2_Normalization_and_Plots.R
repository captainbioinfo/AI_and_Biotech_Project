# =========================================
# Differential Expression Analysis - DESeq2
# =========================================

# ----------------------------
# Load libraries
# ----------------------------
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyverse)
library(dplyr)
# ----------------------------
# Load Data
# ----------------------------
counts <- read.csv("C:/Users/Dell/Downloads/gene_counts_cleaned.csv",
                   header = TRUE, row.names = 1, check.names = FALSE)

meta <- read.csv("C:/Users/Dell/Downloads/metadata (1).csv", 
                 header = TRUE, row.names = 1)

# Add condition column
meta$condition <- factor(ifelse(meta$tissue == "Tumor", "Tumor", "Normal"))

# Reorder meta to match counts columns
meta <- meta[colnames(counts), ]

# ----------------------------
# DESeq2 Dataset
# ----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta,
                              design    = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Set reference level
dds$condition <- relevel(dds$condition, ref = "Normal")

# ----------------------------
# Run DESeq2
# ----------------------------
dds <- DESeq(dds)

# ----------------------------
# Extract Results
# ----------------------------
res <- results(dds, alpha = 0.05)
res <- res[order(res$pvalue), ]

# Subset significant DEGs |log2FC| > 1
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Save results
write.csv(as.data.frame(res), "DEG_results.csv")
write.csv(as.data.frame(sig_res), "DEG_significant_log2FC1.csv")

# ----------------------------
# Normalized Counts
# ----------------------------
norm_counts <- counts(dds, normalized = TRUE)
rownames(norm_counts) <- rownames(dds)  # ensure rownames match
write.csv(norm_counts, "normalized_counts.csv")

# ----------------------------
# PCA Plots
# ----------------------------
rld <- vst(dds)

# PCA PDF
pdf("PCA_plot.pdf", width = 6, height = 6)
plotPCA(rld, intgroup = "condition")
dev.off()

# PCA PNG
png("PCA_plot.png", width = 3000, height = 3000, res = 300)
plotPCA(rld, intgroup = "condition")
dev.off()

# ----------------------------
# Volcano Plot
# ----------------------------
volcano_data <- as.data.frame(res)
volcano_data <- volcano_data %>%
  mutate(Expression = case_when(
    log2FoldChange >= 1  & padj <= 0.05 ~ "Up-regulated",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
    TRUE                                  ~ "Not Significant"
  ))

pdf("volcano_plot.pdf", width = 7, height = 6)

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Expression), alpha = 0.8, size = 1.6) +
  scale_color_manual(
    values = c(
      "Up-regulated"     = "firebrick3",
      "Down-regulated"   = "dodgerblue3",
      "Not Significant"  = "grey40"
    )
  ) +
geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", color = "grey50") +
  
  xlab(expression("log"[2] * " Fold Change")) +
  ylab(expression("-log"[10] * " P-value")) +
  theme_minimal(base_size = 13) +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  ggtitle("Volcano Plot of Differentially Expressed Genes")

dev.off()

# Volcano PNG
png("volcano_plot.png", width = 3000, height = 3000, res = 300)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Expression), alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c("Down-regulated" = "dodgerblue3",
                                "Not Significant" = "grey40",
                                "Up-regulated" = "firebrick3")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  xlab(expression("log"[2]*" Fold Change")) +
  ylab(expression("-log"[10]*" P-value")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(face = "bold")) +
  ggtitle("Volcano Plot of Differentially Expressed Genes")
dev.off()

# ----------------------------
# MA Plot
# ----------------------------
ma_data <- as.data.frame(res)

ma_data$Expression <- dplyr::case_when(
  ma_data$log2FoldChange >= 1  & ma_data$padj <= 0.05 ~ "Up-regulated",
  ma_data$log2FoldChange <= -1 & ma_data$padj <= 0.05 ~ "Down-regulated",
  TRUE ~ "Not Significant"
)

pdf("MA_plot.pdf", width = 7, height = 6)

ggplot(ma_data, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = Expression), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "Down-regulated" = "dodgerblue3",
    "Not Significant" = "black",
    "Up-regulated" = "firebrick3"
  )) +
  scale_x_log10() +
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "grey40") +
  xlab("Mean Expression (baseMean, log scale)") +
  ylab(expression(log[2]~"Fold Change")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )

dev.off()


# MA plot PNG
png("MA_plot.png", width = 3000, height = 2400, res = 300)

ggplot(ma_data, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = Expression), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "Down-regulated" = "dodgerblue3",
    "Not Significant" = "black",
    "Up-regulated" = "firebrick3"
  )) +
  scale_x_log10() +
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "grey40") +
  xlab("Mean Expression (baseMean, log scale)") +
  ylab(expression(log[2]~"Fold Change")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )

dev.off()

# -------------------------------
# Heatmap of Top 50 DEGs
# -------------------------------
# Significant DEGs (padj <= 0.05)
sig_genes <- rownames(res[res$padj <= 0.05 & !is.na(res$padj), ])

# Transform normalized counts for heatmap
vsd <- vst(dds, blind = FALSE)   # or rlog(dds)
heat_data <- assay(vsd)[sig_genes, ]

#Annotation for columns (samples)
annot <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
colnames(annot) <- "Condition"

# Filter Significant DEGs
sig_df <- as.data.frame(res) %>% 
  filter(padj <= 0.05) %>% 
  arrange(desc(abs(log2FoldChange)))

top50_genes <- rownames(sig_df)[1:50]   # Top 50 by |log2FC|
heat_top50 <- assay(vsd)[top50_genes, ]

top50_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)


# PDF
pdf("heatmap_top50_DEGs.pdf", width = 10, height = 12)
pheatmap(
  heat_top50,
  scale = "row",
  annotation_col = annot,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  color = top50_colors
)
dev.off()

# PNG
png("heatmap_top50_DEGs.png", width = 4000, height = 4000, res = 300)
pheatmap(
  heat_top50,
  scale = "row",
  annotation_col = annot,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 6,
  color = top50_colors
)
dev.off()



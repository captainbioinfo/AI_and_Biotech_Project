# Functional Enrichment Analysis - Colorectal Cancer Project
# Fahad - November 2024

# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
install.packages(c("ggplot2", "dplyr", "RColorBrewer", "viridis", "ggnewscale"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("C:/Users/DELL/Documents/Ai and Biotechnology project")

# Read the DEG data
deg_data <- read.csv("DEG_results.csv", header = TRUE)
head(deg_data)

# Get significant genes (padj < 0.05 and |log2FC| > 1)
sig_genes <- deg_data %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

# Extract Ensembl IDs (column X has the IDs)
genes <- unique(sig_genes$X)
print(paste("Found", length(genes), "significant genes"))

# Convert Ensembl IDs to Entrez IDs
entrez_ids <- bitr(genes, 
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

print(paste("Converted", nrow(entrez_ids), "genes to Entrez IDs"))

# GO enrichment - Biological Process
go_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# GO enrichment - Molecular Function
go_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# GO enrichment - Cellular Component
go_cc <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# KEGG pathway analysis
kegg_result <- enrichKEGG(gene = entrez_ids$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Save results to CSV
write.csv(as.data.frame(go_bp), "GO_BP_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "GO_MF_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "GO_CC_results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_result), "KEGG_results.csv", row.names = FALSE)

# Create publication-quality plots with dark colors - PDF format
# GO BP Dot Plot
pdf("Figure_GO_BP_dotplot.pdf", width = 14, height = 10)
dotplot(go_bp, showCategory = 20, font.size = 10) +
  scale_color_gradientn(colours = c("#00264d", "#0066cc", "#cc0000"),
                        name = "p.adjust") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"))
dev.off()

# GO BP Bar Plot
pdf("Figure_GO_BP_barplot.pdf", width = 12, height = 9)
barplot(go_bp, showCategory = 15, font.size = 10) +
  scale_fill_gradient(low = "#003366", high = "#990000", name = "p.adjust") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.position = "right")
dev.off()

# Gene network plot
pdf("Figure_GO_BP_network.pdf", width = 16, height = 13)
cnetplot(go_bp, showCategory = 5, node_label = "all", 
         cex_label_category = 0.9, cex_label_gene = 0.6, colorEdge = TRUE) +
  scale_color_gradient(low = "#003366", high = "#990000") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9))
dev.off()

# Enrichment map
go_sim <- pairwise_termsim(go_bp)
pdf("Figure_GO_BP_emap.pdf", width = 14, height = 12)
emapplot(go_sim, showCategory = 25, layout = "nicely") +
  scale_fill_gradientn(colours = c("#00264d", "#004080", "#cc0000"),
                       name = "p.adjust") +
  theme(legend.position = "right",
        legend.text = element_text(size = 9))
dev.off()

# KEGG dot plot
pdf("Figure_KEGG_dotplot.pdf", width = 14, height = 10)
dotplot(kegg_result, showCategory = 20, font.size = 10) +
  scale_color_gradientn(colours = c("#00264d", "#0066cc", "#cc0000"),
                        name = "p.adjust") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.position = "right")
dev.off()

# KEGG bar plot
pdf("Figure_KEGG_barplot.pdf", width = 12, height = 9)
barplot(kegg_result, showCategory = 15, font.size = 10) +
  scale_fill_gradient(low = "#003366", high = "#990000", name = "p.adjust") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11, face = "bold"))
dev.off()

# Separate up and down regulated genes
up_genes <- deg_data %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  pull(X) %>%
  unique()

down_genes <- deg_data %>%
  filter(padj < 0.05, log2FoldChange < -1) %>%
  pull(X) %>%
  unique()

print(paste("Upregulated genes:", length(up_genes)))
print(paste("Downregulated genes:", length(down_genes)))

# Convert to Entrez IDs
up_entrez <- bitr(up_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment for upregulated genes
go_up <- enrichGO(gene = up_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# GO enrichment for downregulated genes
go_down <- enrichGO(gene = down_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    readable = TRUE)

# Save separate results
write.csv(as.data.frame(go_up), "GO_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), "GO_downregulated.csv", row.names = FALSE)

# Compare up vs down regulated genes
png("Figure_UP_vs_DOWN_comparison.png", width = 14, height = 10, units = "in", res = 600)
gene_lists <- list(Upregulated = up_entrez$ENTREZID,
                   Downregulated = down_entrez$ENTREZID)
comparison <- compareCluster(geneCluster = gene_lists, 
                             fun = enrichGO,
                             OrgDb = org.Hs.eg.db,
                             ont = "BP")
dotplot(comparison, font.size = 10) +
  scale_color_gradientn(colours = c("#00264d", "#0066cc", "#cc0000")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"))
dev.off()

print("Analysis complete! Check your folder for results and figures.")

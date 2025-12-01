# =====================================================================
# DIRECT GO AND KEGG ENRICHMENT ANALYSIS
# Direct gene enrichment
# =====================================================================

# Set user library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}
.libPaths(c(user_lib, .libPaths()))

# Load required libraries
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", 
                       "ggplot2", "dplyr", "RColorBrewer")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Set working directory and create results folder
results_folder <- "C:/Users/DELL/Downloads/GO_KEGG_Enrichment_Results"
if (!dir.exists(results_folder)) {
  dir.create(results_folder, recursive = TRUE)
}
setwd(results_folder)

cat("\n========== GO AND KEGG ENRICHMENT ANALYSIS ==========\n")
cat("(No normalization, No DEG creation)\n\n")

# =====================================================================
# STEP 1: READ GENE LIST
# =====================================================================
cat("STEP 1: Reading gene list...\n")

# Read from gene_counts_ml.csv
gene_file <- "C:/Users/DELL/Downloads/gene_counts_ml.csv"

if (file.exists(gene_file)) {
  cat("  Reading from gene_counts_ml.csv\n")
  count_data <- read.csv(gene_file, header = TRUE)
  
  # Get gene IDs (first column)
  gene_ids <- count_data[, 1]
  
  # Remove duplicates if any
  gene_ids <- unique(gene_ids)
  
  cat("✓ Total genes loaded:", length(gene_ids), "\n\n")
} else {
  stop("gene_counts_ml.csv not found at C:/Users/DELL/Downloads/")
}

# =====================================================================
# STEP 2: CONVERT ENSEMBL IDS TO ENTREZ IDS
# =====================================================================
cat("STEP 2: Converting Ensembl IDs to Entrez IDs...\n")

entrez_ids <- bitr(gene_ids, 
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

cat("✓ Successfully converted:", nrow(entrez_ids), "genes to Entrez IDs\n")
cat("  (Original genes: ", length(gene_ids), ")\n")
cat("  (Conversion rate: ", round(nrow(entrez_ids)/length(gene_ids)*100, 1), "%)\n\n")

if (nrow(entrez_ids) == 0) {
  stop("No genes could be converted to Entrez IDs")
}

# =====================================================================
# STEP 3: GO ENRICHMENT ANALYSIS
# =====================================================================
cat("STEP 3: Performing GO enrichment analysis...\n")

# GO - Biological Process
cat("  • GO Biological Process enrichment...\n")
go_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

cat("    Found:", nrow(as.data.frame(go_bp)), "BP terms\n")

# GO - Molecular Function
cat("  • GO Molecular Function enrichment...\n")
go_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

cat("    Found:", nrow(as.data.frame(go_mf)), "MF terms\n")

# GO - Cellular Component
cat("  • GO Cellular Component enrichment...\n")
go_cc <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

cat("    Found:", nrow(as.data.frame(go_cc)), "CC terms\n\n")

# =====================================================================
# STEP 4: KEGG PATHWAY ENRICHMENT (ROBUST VERSION)
# =====================================================================
cat("STEP 4: Performing KEGG pathway enrichment...\n")

# Increase timeout for KEGG
options(timeout = 600)

kegg_result <- tryCatch({
  cat("  Attempting KEGG enrichment (Attempt 1/3)...\n")
  enrichKEGG(gene = entrez_ids$ENTREZID,
             organism = "hsa",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
}, error = function(e) {
  cat("  KEGG Attempt 1 failed. Retrying with relaxed parameters...\n")
  
  tryCatch({
    cat("  Attempting KEGG enrichment (Attempt 2/3)...\n")
    enrichKEGG(gene = entrez_ids$ENTREZID,
               organism = "hsa",
               pvalueCutoff = 0.1,        # Relaxed p-value
               qvalueCutoff = 0.3)        # Relaxed q-value
  }, error = function(e2) {
    cat("  KEGG Attempt 2 failed. Trying alternative approach...\n")
    
    tryCatch({
      cat("  Attempting KEGG enrichment (Attempt 3/3)...\n")
      enrichKEGG(gene = entrez_ids$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.5,
                 minGSSize = 1)
    }, error = function(e3) {
      cat("  ✗ KEGG enrichment failed after 3 attempts\n")
      cat("  Error:", e3$message, "\n")
      return(NULL)
    })
  })
})

if (!is.null(kegg_result)) {
  kegg_result_df <- as.data.frame(kegg_result)
  
  if (nrow(kegg_result_df) > 0) {
    kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    cat("✓ KEGG enrichment SUCCESSFUL!\n")
    cat("✓ Found:", nrow(kegg_result), "KEGG pathways\n\n")
  } else {
    cat("! KEGG returned results object but with 0 pathways\n")
    cat("  This may mean genes don't significantly match KEGG pathways\n\n")
    kegg_result <- NULL
  }
} else {
  cat("! KEGG enrichment did not produce results\n")
  cat("  Possible causes:\n")
  cat("    - Network connectivity issues\n")
  cat("    - KEGG API timeout\n")
  cat("    - Genes don't match KEGG database\n\n")
  kegg_result <- NULL
}

# =====================================================================
# STEP 5: SAVE RESULTS TO CSV
# =====================================================================
cat("STEP 5: Saving results to CSV...\n")

write.csv(as.data.frame(go_bp), "GO_Biological_Process.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "GO_Molecular_Function.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "GO_Cellular_Component.csv", row.names = FALSE)

cat("✓ GO results saved\n")

if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
  write.csv(as.data.frame(kegg_result), "KEGG_Pathways.csv", row.names = FALSE)
  cat("✓ KEGG results saved\n")
}

cat("\n")

# =====================================================================
# STEP 6: CREATE VISUALIZATIONS
# =====================================================================
cat("STEP 6: Creating visualizations...\n")

# Color scheme - Red and Blue gradient
neat_gradient <- c("#2E5090", "#5B9BD5", "#ED7D31", "#C5504C")
neat_blue <- "#5B9BD5"
neat_red <- "#E60000"

# ===== GO BIOLOGICAL PROCESS =====
cat("  • Creating GO BP dot plot...\n")
pdf("01_GO_BP_DotPlot.pdf", width = 18, height = 11)
p1 <- dotplot(go_bp, showCategory = 20, font.size = 12) +
  scale_color_gradientn(colours = neat_gradient, name = "Adjusted\np-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.text.x = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major.x = element_line(colour = "lightgrey", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
print(p1)
dev.off()

cat("  • Creating GO BP bar plot...\n")
pdf("02_GO_BP_BarPlot.pdf", width = 16, height = 12)
p2 <- barplot(go_bp, showCategory = 15, font.size = 12) +
  scale_fill_gradient(low = neat_blue, high = neat_red, name = "Adjusted\np-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.text.x = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid.major.x = element_line(colour = "lightgrey", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
print(p2)
dev.off()

cat("  • Creating GO BP network plot...\n")
pdf("03_GO_BP_Network.pdf", width = 18, height = 14)
p3 <- cnetplot(go_bp, showCategory = 8, node_label = "category",
               cex_label_category = 1.3, cex_label_gene = 0.9, colorEdge = TRUE) +
  scale_color_gradient(low = neat_blue, high = neat_red) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
print(p3)
dev.off()

cat("  • Creating GO BP enrichment map...\n")
go_sim <- pairwise_termsim(go_bp)
pdf("04_GO_BP_EnrichmentMap.pdf", width = 18, height = 14)
p4 <- emapplot(go_sim, showCategory = 25, layout = "nicely") +
  scale_fill_gradientn(colours = neat_gradient, name = "Adjusted\np-value") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"))
print(p4)
dev.off()

# ===== GO MOLECULAR FUNCTION =====
cat("  • Creating GO MF dot plot...\n")
pdf("05_GO_MF_DotPlot.pdf", width = 18, height = 11)
p5 <- dotplot(go_mf, showCategory = 15, font.size = 12) +
  scale_color_gradientn(colours = neat_gradient, name = "Adjusted\np-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.text.x = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major.x = element_line(colour = "lightgrey", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
print(p5)
dev.off()

# ===== GO CELLULAR COMPONENT =====
cat("  • Creating GO CC dot plot...\n")
pdf("06_GO_CC_DotPlot.pdf", width = 18, height = 11)
p6 <- dotplot(go_cc, showCategory = 15, font.size = 12) +
  scale_color_gradientn(colours = neat_gradient, name = "Adjusted\np-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.text.x = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major.x = element_line(colour = "lightgrey", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
print(p6)
dev.off()

# ===== KEGG PATHWAYS =====
if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
  
  cat("  • Creating KEGG dot plot...\n")
  pdf("07_KEGG_DotPlot.pdf", width = 18, height = 11)
  p7 <- dotplot(kegg_result, showCategory = 15, font.size = 12) +
    scale_color_gradientn(colours = neat_gradient, name = "Adjusted\np-value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),
          axis.text.x = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 11),
          axis.title = element_text(size = 14, face = "bold"),
          panel.grid.major.x = element_line(colour = "lightgrey", size = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  print(p7)
  dev.off()
  
  cat("  • Creating KEGG bar plot...\n")
  pdf("08_KEGG_BarPlot.pdf", width = 16, height = 12)
  p8 <- barplot(kegg_result, showCategory = 12, font.size = 12) +
    scale_fill_gradient(low = neat_blue, high = neat_red, name = "Adjusted\np-value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),
          axis.text.x = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 11),
          panel.grid.major.x = element_line(colour = "lightgrey", size = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
  print(p8)
  dev.off()
}

cat("✓ Visualizations created\n\n")

# =====================================================================
# STEP 7: FINAL SUMMARY
# =====================================================================
cat("\n")
cat("=========================================\n")
cat("===  ENRICHMENT ANALYSIS COMPLETE  ===\n")
cat("=========================================\n")
cat("Total genes analyzed:", length(gene_ids), "\n")
cat("Genes converted to Entrez ID:", nrow(entrez_ids), "\n")
cat("\nGO Results:\n")
cat("  • Biological Process terms:", nrow(as.data.frame(go_bp)), "\n")
cat("  • Molecular Function terms:", nrow(as.data.frame(go_mf)), "\n")
cat("  • Cellular Component terms:", nrow(as.data.frame(go_cc)), "\n")

if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
  cat("\nKEGG Results:\n")
  cat("  • Pathways found:", nrow(as.data.frame(kegg_result)), "\n")
}

cat("\n=========================================\n")
cat("Output directory:", getwd(), "\n")
cat("=========================================\n\n")

cat("CSV Files Generated:\n")
cat("  ✓ GO_Biological_Process.csv\n")
cat("  ✓ GO_Molecular_Function.csv\n")
cat("  ✓ GO_Cellular_Component.csv\n")
if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
  cat("  ✓ KEGG_Pathways.csv\n")
}

cat("\nPDF Visualizations Generated:\n")
cat("  ✓ 01_GO_BP_DotPlot.pdf\n")
cat("  ✓ 02_GO_BP_BarPlot.pdf\n")
cat("  ✓ 03_GO_BP_Network.pdf\n")
cat("  ✓ 04_GO_BP_EnrichmentMap.pdf\n")
cat("  ✓ 05_GO_MF_DotPlot.pdf\n")
cat("  ✓ 06_GO_CC_DotPlot.pdf\n")
if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
  cat("  ✓ 07_KEGG_DotPlot.pdf\n")
  cat("  ✓ 08_KEGG_BarPlot.pdf\n")
}

cat("\n=========================================\n")
cat("✓ ALL TASKS COMPLETED SUCCESSFULLY!\n")
cat("=========================================\n\n")

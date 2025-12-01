# AI_and_Biotech_Project
### *Exploring RNA-Seq derived Biomarkers for colorectal cancer classification via Machine Learning.*

---


# 👥 Project Members

- **Yusuf Munir Aliyu**  
- **Saniya Khurshid**  
- **Fahad Sajjad**  
- **Farha Tarique**
---



# 📌 Project Overview

Colorectal cancer (CRC) is the third most prevalent malignancy worldwide and accounts for nearly **10% of all cancer-related deaths**. Accurate early detection remains a challenge due to tumour heterogeneity and lack of reliable biomarkers.

This project integrates:

- **RNA-seq preprocessing and normalisation**
- **Differential gene expression (DEG) analysis**
- **Feature selection (LASSO, ROC-AUC screening)**
- **Functional enrichment (GO, KEGG)**
- **Machine-learning classification (SVM)**

The pipeline leverages both **Galaxy** (for raw FASTQ analysis) and **R/Python** for computational downstream analysis.

### 🔬 Hybrid Workflow Used:
✔ **Galaxy** → QC → trimming → alignment → quantification  
✔ **R + Python** → DEG → biomarker selection → ML classification  

---

# 🎯 Objectives

1. **Identify differentially expressed genes (DEGs)** between CRC tumour and normal tissues.  
2. **Screen biomarker genes** using:
   - LASSO regression  
   - ROC analysis  
   - ML feature ranking  
3. **Build supervised ML models**:
      - Support Vector Machine (SVM)  
    
4.  
5. **Provide a reproducible bioinformatics workflow** for future CRC biomarker research.

---

# 🧬 Workflow

<img width="600" height="337" alt="image" src="https://github.com/user-attachments/assets/425fd34d-c7cd-49ca-a440-da8e55251b3c" />


# 📁 Dataset Information

**Accession ID:** GSE156451  
**Samples:** 144 total  
- 72 colorectal cancer tumour tissues  
- 72 matched adjacent normal tissues  

**Platform:** Illumina RNA-seq  
**Data type:** Raw FASTQ + processed gene counts  

Raw dataset available from GEO:  
🔗 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156451

---

# 📂 Repository Structure


```
CRC-Biomarker-Discovery/
│
├── 📁 Data/
│   ├── gene_counts_cleaned.csv
│   ├── metadata.csv
│   └── README.md
│
├── 📁 References/
│   └── Tool_References.md
│
├── 📁 Scripts/
│   ├── DESeq2_Normalization_and_Plotting.R
│   ├── GO_and_KEGG.R
│   ├── ML_script.ipynb
│   ├── biomarker_figures_script.Rmd
│   ├── to_merge_featurecounts.py
│   └── README.md
│
├── 📁 Results/
│   ├── 📁 DEG/
│   │   ├── volcano_plot.png
│   │   ├── heatmap_top50.png
│   │   └── DESeq2_results.csv
│   │
│   ├── 📁 Enrichment_Analysis/
│   │   ├── GO_BP_MF_CC.csv
│   │   ├── KEGG_pathways.csv
│   │   └── enrichment_plots/
│   │
│   ├── 📁 ML/
│   │   ├── LASSO_results.csv
│   │   ├── SVM_RFE_results.csv
│   │   ├── core_genes.csv
│   │   ├── ROC_curves.png
│   │   └── stability_scores.csv
│   │
│   ├── 📁 QC/
│   │   ├── fastqc_reports/
│   │   └── multiqc_report.html
│
├── 📁 Figures/
│   ├── workflow_diagram.png
│   ├── PCA_UMAP.png
│   └── biomarker_violin_density.png
│
├── 📁 Docs/
│   ├── Project_Overview.md
│   ├── Pipeline_Workflow.md
│   ├── QC_Guidelines.md
│   └── References.md
│
├── LICENSE
└── README.md

```
---

# 🧰 Tools Used

### 🧪 **RNA-seq Processing (Galaxy platform)**  
- **FastQC** — Quality control  
- **Fastp** — Trimming and filtering  
- **MultiQC** — QC report summary  
- **Hisat2** — Genome alignment  
- **FeatureCounts** — Gene-level quantification  

### 📊 **R Packages**
- **DESeq2** — Differential expression  
- **ClusterProfiler** — GO/KEGG enrichment  
- **pROC** — ROC curve analysis  
- **glmnet** — LASSO regression  

### 🤖 **Python Packages**
- numpy  
- pandas  
- scipy  
- scikit-learn  
- matplotlib / seaborn  

---
📄 License
This project is released under the MIT License.
You are free to use, modify, and distribute with attribution.

📣 Citation
Farha, T., Munir, Y. A., Saniya, K., Fahad, S. 
AI_and_Biotech_Project: ----

📬 Contact
Email: 
        
    

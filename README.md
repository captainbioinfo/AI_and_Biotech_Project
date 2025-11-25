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

- **RNA-seq preprocessing**
- **Differential gene expression (DEG) analysis**
- **Feature selection (LASSO, ROC-AUC screening)**
- **Machine-learning classification (SVM)**
- **Functional enrichment (GO, KEGG)**

The pipeline leverages both **Galaxy** (for raw FASTQ analysis) and **R/Python** for computational downstream analysis.

### 🔬 Hybrid Workflow Used:
✔ **Galaxy** → QC → trimming → alignment → quantification  
✔ **R + Python** → DEG → biomarker selection → ML classification → immune analysis  

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

<img width="800" height="450" alt="workflow" src="https://github.com/user-attachments/assets/7fcbff94-fa59-4e86-b7bd-09ecefd83c1b" />

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
├── README.md
│
├── galaxy_history_export.html/
│   ├── fastp_results
│   ├── feature_count_results
│   ├── gene_counts_matrix
│   ├── alignments/
│   └── hisat2_results
│
│       
├── data/
│   ├── metadata.csv
│   ├── processed/
│   │   ├── normalized_counts.csv
│   │   ├── DEGs_results.csv
│   │   └── biomarker_genes.csv
│   └── raw/ (empty – raw FASTQs stored in Galaxy)
│
├── scripts/
│   ├── 01_to_merge_featurecounts.py
│   ├── 02_LASSO_biomarker_selection.ipynb
│   ├── 03_ML_classification_models.ipynb
│   ├── 04_ROC_evaluation.ipynb
    ├── 01_deseq2_DEG_analysis.R
│   └── 05_immune_infiltration_analysis.ipynb
│
├── results/
│   ├── DEG_plots/
│   ├── ML_performance/
│   ├── ROC_curves/
│   ├── biomarker_analysis/
│   └── immune_infiltration/
│
└── environment/
    ├── requirements.txt
    └── environment.yml
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
Farha T, Yusuf S, Saniya, Fahad. 
AI_and_Biotech_Project: ----

📬 Contact
Email: 
        
    

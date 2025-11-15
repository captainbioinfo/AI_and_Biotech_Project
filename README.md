# Project Overview

Colorectal cancer (CRC) is characterised by extensive transcriptomic dysregulation that can be effectively profiled using RNA-sequencing (RNA-seq). This project integrates RNA-seq data processing, differential gene expression analysis, functional enrichment, and machine-learning classification to identify biomarker genes capable of distinguishing colorectal tumor tissues from matched normal samples. Using the publicly available **GSE156451** dataset, this project establishes a reproducible workflow for precision oncology research and biomarker discovery.
In this study, RNA-seq data from GSE156451 (72 tumor and 72 matched normal tissues) were processed using a hybrid workflow:

✔ Galaxy for RNA-seq QC → trimming → alignment → quantification
✔ R & Python for downstream analysis (DEG, biomarkers, ML models, immune infiltration)

# Objectives
Identify differentially expressed genes (DEGs)

Screen for biomarker genes using LASSO, ROC, and ML-driven feature ranking

Build cancer classification models (RF, SVM, ANN, GBM)

Evaluate immune infiltration patterns

Provide a reproducible analysis pipeline for future CRC biomarker research
#
<img width="800" height="450" alt="image" src="https://github.com/user-attachments/assets/7fcbff94-fa59-4e86-b7bd-09ecefd83c1b" />


# Dataset Information

Accession ID: GSE156451

Samples: 144 total

72 colorectal cancer tumor samples

72 matched adjacent normal tissues

Platform: Illumina RNA-seq

Type: Raw FASTQ files + processed counts

Raw data can be downloaded from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156451
GSE156451 – NCBI Gene Expression Omnibus.


# Project Members-
Yusuf, Bilal, Saniya, Fahad and Farha 

# Repository Structure

( ``` ) 
├── README.md
│
├── galaxy_workflow/
│   ├── galaxy_history_export.html
│   ├── qc_reports/
│   │   ├── fastqc_raw/
│   │   └── fastqc_trimmed/
│   ├── trimmed_reads/
│   ├── alignments/
│   │   ├── *.bam
│   ├── counts/
│       ├── gene_counts.tsv
│       └── featurecounts_summary.txt
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
│   ├── 01_deseq2_DEG_analysis.R
│   ├── 02_LASSO_biomarker_selection.ipynb
│   ├── 03_ML_classification_models.ipynb
│   ├── 04_ROC_evaluation.ipynb
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
( ``` )
## Tools Used
The following tools and packages were used in this RNA-seq analysis pipeline:

    FastQC – for raw sequence quality control

    Fastp – for trimming and filtering of reads

    MultiQC – for aggregating QC reports

    Hisat2 – for efficient alignment of RNA-seq reads

    FeatureCounts – for efficient read summarisation

    R Packages:

        DESeq2 – differential expression analysis.
        ClusterProfiler-
        
    Python Packages:
        numpy
        pandas
        scipy
        scikit-learn


        
    

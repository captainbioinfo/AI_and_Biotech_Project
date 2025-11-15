# AI_and_Biotech_Project

Colorectal cancer (CRC) is characterised by extensive transcriptomic dysregulation that can be effectively profiled using RNA-sequencing (RNA-seq). This project integrates RNA-seq data processing, differential gene expression analysis, functional enrichment, and machine-learning classification to identify biomarker genes capable of distinguishing colorectal tumor tissues from matched normal samples. Using the publicly available **GSE156451** dataset, this project establishes a reproducible workflow for precision oncology research and biomarker discovery.

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

# Project Members-
Yusuf, Bilal, Saniya, Fahad and Farha 

# Requirements

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


        
    

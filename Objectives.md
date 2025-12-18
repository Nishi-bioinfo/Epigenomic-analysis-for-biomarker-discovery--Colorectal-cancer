## Project Objectives

The objectives of this project are designed to align with the implemented whole-exome sequencing (WES) pipeline and the broader goal of colorectal cancer biomarker discovery. The specific objectives are to:

- Perform **quality control and preprocessing** of raw paired-end sequencing data using standard NGS quality assessment and trimming tools.
- Align high-quality reads to the **human reference genome (hg38)** and generate analysis-ready alignment files following **GATK best practices**.
- Apply **post-alignment processing**, including duplicate marking and **base quality score recalibration (BQSR)**, to improve variant-calling accuracy.
- Conduct **chromosome-specific variant discovery (chromosome 12)** using **GATK HaplotypeCaller**, with a targeted focus on the **KRAS gene**, a key biomarker in colorectal cancer.
- Identify and filter **high-confidence single nucleotide variants (SNVs) and small insertions/deletions (INDELs)** relevant to colorectal cancer.
- Prepare variant datasets for **functional annotation and downstream biological interpretation**.
- Establish a **reproducible and scalable bioinformatics workflow** that can be extended to additional cancer-associated genes or integrated with epigenomic and transcriptomic data in future studies.

---

## Methodological Approach

This project primarily implements a **DNA-based whole-exome sequencing (WES) analysis pipeline** following **GATK best practices** to identify high-confidence genetic variants associated with colorectal cancer. The pipeline focuses on chromosome 12 to investigate clinically relevant mutations in the **KRAS gene**.

In addition to genomic variant discovery, the project is framed within a broader **epigenomic and transcriptomic context**. RNA-based high-throughput technologies, such as RNA sequencing (RNA-seq), are discussed as complementary approaches for understanding gene expression changes and epigenetic regulation in colorectal cancer.

By integrating variant-level genomic insights with known epigenetic and transcriptomic mechanisms reported in the literature, this project establishes a foundation for future multi-omics analyses.

This integrative framework supports:

- Early cancer detection
- Risk stratification
- Identification of therapeutic targets
- Development of personalized treatment strategies

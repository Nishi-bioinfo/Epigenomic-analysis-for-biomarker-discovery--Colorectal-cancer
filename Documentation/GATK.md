# Colorectal Cancer Biomarker Discovery using KRAS (Chromosome 12)

## 1. Introduction to Biomarker Discovery in Colorectal Cancer

**Biomarker discovery** is the process of identifying measurable biological indicators (DNA variants, RNA expression, proteins, metabolites) that are associated with disease presence, progression, prognosis, or treatment response.

In **colorectal cancer (CRC)**, genomic biomarkers are crucial for:

* Early diagnosis and risk stratification
* Predicting prognosis
* Guiding targeted therapy (precision medicine)
* Monitoring treatment response and resistance

Somatic variants detected using **next-generation sequencing (NGS)**, particularly **whole-exome sequencing (WES)**, are widely used to identify cancer-associated biomarkers.

---

## 2. What is Colorectal Cancer?

**Colorectal cancer** is a malignancy that arises from the epithelial cells lining the colon or rectum. It is one of the leading causes of cancer-related mortality worldwide.

Key molecular characteristics of CRC include:

* Accumulation of somatic mutations
* Dysregulation of oncogenes and tumor suppressor genes
* Alterations in signaling pathways such as **MAPK**, **PI3K-AKT**, and **WNT**

---

## 3. KRAS Gene and Its Role in Colorectal Cancer

### 3.1 What is KRAS?

* **Gene name:** KRAS (Kirsten rat sarcoma viral oncogene homolog)
* **Chromosomal location:** Chromosome 12 (12p12.1)
* **Gene type:** Proto-oncogene
* **Protein function:** Small GTPase involved in signal transduction

KRAS acts as a molecular switch in the **EGFR–RAS–RAF–MEK–ERK (MAPK) pathway**, regulating:

* Cell proliferation
* Differentiation
* Survival

### 3.2 KRAS as a Biomarker

Mutations in KRAS are among the most frequent driver mutations in CRC.

Clinical relevance:

* KRAS mutations predict **resistance to anti-EGFR therapies** (e.g., cetuximab, panitumumab)
* Used as a **predictive biomarker** in CRC treatment decisions

Common mutation hotspots:

* Codon 12
* Codon 13
* Codon 61

---

## 4. Why Chromosome 12 Was Chosen

* The **KRAS gene is located on chromosome 12**
* Focusing on a single chromosome reduces computational complexity
* Enables targeted analysis while retaining biological relevance
* Suitable for demonstration of a complete GATK-based WES variant discovery pipeline

This project limits alignment, recalibration, and variant calling to **chr12**, allowing efficient biomarker-focused analysis.

---

## 5. Project Objective

The primary objectives of this project are:

1. To process raw NGS data using a standard **GATK best-practices workflow**
2. To identify high-confidence single nucleotide variants (SNVs) and indels
3. To generate analysis-ready variants from **chromosome 12**
4. To prepare variant data for **biomarker discovery in CRC**, focusing on the KRAS gene

---

## 6. Workflow Overview

**Overall pipeline:**

1. Environment and tool setup
2. Raw FASTQ download and subsampling
3. Quality control (FastQC)
4. Reference genome preparation (chr12)
5. Read alignment (BWA-MEM)
6. BAM sorting and indexing
7. Duplicate marking (GATK)
8. Base Quality Score Recalibration (BQSR)
9. Variant calling (HaplotypeCaller)
10. Variant selection (SNPs/INDELs)
11. Variant filtration
12. Annotation preparation (ANNOVAR)

---

## 7. Tools and Technologies

| Category         | Tool       | Purpose                           |
| ---------------- | ---------- | --------------------------------- |
| OS-level         | apt, conda | Package management                |
| QC               | FastQC     | Read quality assessment           |
| Subsampling      | SeqTK      | Downsampling FASTQ reads          |
| Alignment        | BWA-MEM    | Read alignment to reference       |
| BAM processing   | SAMtools   | Sorting and indexing              |
| Variant tools    | BCFtools   | VCF manipulation                  |
| Variant calling  | GATK       | Recalibration & variant discovery |
| Containerization | Docker     | Reproducible GATK execution       |
| Annotation       | ANNOVAR    | Functional annotation             |

---

## 8. Environment Setup

### 8.1 Verify Docker Installation

```bash
docker --version
```

Checks whether Docker is installed and accessible.

### 8.2 System Update and Tool Installation

```bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y fastqc bwa samtools bcftools vcftools docker.io
conda install -c bioconda seqtk -y
```

Installs all required bioinformatics tools.

### 8.3 Pull GATK Docker Image

```bash
sudo docker pull broadinstitute/gatk:latest
```

Downloads the official GATK container.

---

## 9. Data Organization

```bash
mkdir -p Raw_Data Outputs
```

Creates directories for input data and results.

---

## 10. Raw Data Download and Subsampling

### 10.1 FASTQ Download

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/SRR35521082_1.fastq.gz -P Raw_Data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/SRR35521082_2.fastq.gz -P Raw_Data
```

Downloads paired-end sequencing reads.

### 10.2 Subsampling Reads

```bash
seqtk sample -s123 Raw_Data/SRR35521082_1.fastq.gz 300000 > Raw_Data/SRR35521082_1_300k.fastq
seqtk sample -s123 Raw_Data/SRR35521082_2.fastq.gz 300000 > Raw_Data/SRR35521082_2_300k.fastq
```

Reduces dataset size for faster processing.

**Output:**

* `SRR35521082_1_300k.fastq`
* `SRR35521082_2_300k.fastq`

---

## 11. Quality Control

```bash
fastqc -o Outputs Raw_Data/SRR35521082_1_300k.fastq Raw_Data/SRR35521082_2_300k.fastq
```

Generates quality reports.

**Output:**

* HTML and ZIP FastQC reports in `Outputs/`

---

## 12. Reference Genome Preparation (Chromosome 12)

### 12.1 Download Reference

```bash
wget https://ftp.ncbi.nlm.nih.gov/.../chr12.fna.gz
gunzip chr12.fna.gz
mv chr12.fna Outputs/chr12.fa
```

Downloads and prepares chromosome 12 reference sequence.

### 12.2 Indexing

```bash
bwa index chr12.fa
samtools faidx chr12.fa
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk CreateSequenceDictionary -R /data/chr12.fa -O /data/chr12.dict
```

Creates all required reference indices.

**Output:**

* `chr12.fa`
* `chr12.fa.fai`
* `chr12.dict`

---

## 13. Alignment

```bash
bwa mem -t 4 -R "@RG\tID:SRR35521082\tPL:ILLUMINA\tSM:SRR35521082" \
  chr12.fa ../Raw_Data/SRR35521082_1_300k.fastq \
  ../Raw_Data/SRR35521082_2_300k.fastq | \
  samtools sort -o SRR35521082.sorted.bam
```

Aligns reads and sorts alignments.

**Output:**

* `SRR35521082.sorted.bam`

---

## 14. Duplicate Marking

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk MarkDuplicates \
  -I /data/SRR35521082.sorted.bam \
  -O /data/SRR35521082.sorted.dedup.bam \
  -M /data/SRR35521082.metrics.txt
```

Marks PCR duplicates to avoid biased variant calls.

**Output:**

* `SRR35521082.sorted.dedup.bam`
* `SRR35521082.metrics.txt`

---

## 15. Base Quality Score Recalibration (BQSR)

### 15.1 Known Sites Preparation

```bash
bcftools view -r chr12 Homo_sapiens_assembly38.dbsnp138.vcf.gz -Oz -o known_sites_chr12.vcf.gz
```

Extracts known variants for chr12.

### 15.2 Recalibration

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk BaseRecalibrator \
  -I /data/SRR35521082.sorted.dedup.bam \
  -R /data/chr12.fa \
  --known-sites /data/known_sites_chr12_NC.vcf.gz \
  -O /data/SRR35521082.recal_data.table
```

Generates recalibration model.

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk ApplyBQSR \
  -I /data/SRR35521082.sorted.dedup.bam \
  -R /data/chr12.fa \
  --bqsr-recal-file /data/SRR35521082.recal_data.table \
  -O /data/SRR35521082.sorted.dedup.BQSR.bam
```

**Output:**

* `SRR35521082.sorted.dedup.BQSR.bam`

---

## 16. Variant Calling

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk HaplotypeCaller \
  -R /data/chr12.fa \
  -I /data/SRR35521082.sorted.dedup.BQSR.bam \
  -O /data/SRR35521082.raw.vcf
```

**Output:**

* `SRR35521082.raw.vcf`

---

## 17. Variant Filtering

### 17.1 SNP Selection and Filtration

```bash
gatk SelectVariants --select-type SNP
gatk VariantFiltration
```

Applies hard filters to remove low-quality variants.

**Output:**

* `SRR35521082.analysis_ready_snps.vcf`

---

## 18. Annotation Preparation

```bash
convert2annovar.pl -format vcf4 SRR35521082.analysis_ready_snps.vcf
```

Prepares variants for functional annotation.

**Output:**

* `SRR35521082.analysis_ready_snps.avinput`

---

## 19. Final Outcome

This pipeline produces **high-confidence variants on chromosome 12**, suitable for:

* Identifying KRAS mutations
* Biomarker discovery in colorectal cancer
* Downstream annotation and clinical interpretation
* Deleted the data due to storage issue
---

## 20. Conclusion

This project demonstrates a **reproducible, chromosome-focused GATK pipeline** for CRC biomarker discovery, highlighting the clinical importance of **KRAS mutations** and providing a strong foundation for translational cancer genomics research.

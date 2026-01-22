# 🧬 SRR35521082 Variant Annotation Pipeline

**Complete Step-by-Step pipeline for Colorectal Cancer Biomarker Discovery**

---

## 📋 Table of Contents

1. [Project Overview](#1-project-overview)
2. [Input Files](#2-input-files)
3. [Step 1: VCF Header Inspection](#3-step-1--vcf-header-inspection)
4. [Step 2: Variant Record Examination](#4-step-2--variant-record-examination)
5. [Step 3: BGZIP Re-compression](#5-step-3--bgzip-re-compression)
6. [Step 4: VCF Indexing](#6-step-4--vcf-indexing)
7. [Step 5: Ensembl VEP Annotation](#7-step-5--variant-annotation-with-ensembl-vep)
8. [Step 6: Cancer Variant Classification](#8-step-6--cancer-variant-classification-script)
9. [Step 7: Final Variant Annotation](#9-step-7--final-variant-annotation-script)
10. [Pipeline Summary](#11-pipeline-summary-flow)

---

## 1. Project Overview

This pipeline processes **GATK-called SNP variants** from Chromosome 12 (GRCh38) and performs:

* Variant quality validation
* BGZIP re-indexing
* Ensembl VEP annotation (gene, protein, clinical significance)
* Cancer gene extraction (KRAS, TP53, PIK3CA)
* Pathway classification and biomarker prioritization

---

## 2. Input Files

### Primary Input File

```
Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

### VCF Column Structure

| Column | Description                       |
| ------ | --------------------------------- |
| CHROM  | Chromosome (NC_000012.12 = chr12) |
| POS    | Genomic position                  |
| REF    | Reference allele                  |
| ALT    | Alternate allele                  |
| QUAL   | Variant confidence score          |
| FILTER | PASS or filter reason             |
| INFO   | Depth, AF, strand bias, etc.      |
| FORMAT | Genotype fields (GT:AD:DP:GQ:PL)  |
| SAMPLE | Sample-specific values            |

---

## 3. Step 1 — VCF Header Inspection

### Command

```bash
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | head -20
```

### Expected Output (Key Lines)

```text
##fileformat=VCFv4.2
##source="GATK HaplotypeCaller 4.6.2.0"
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
```

### Key Learnings

* GATK v4.6.2.0 with GRCh38 reference
* Hard filters applied: FS, MQ, QD, SOR
* Clinical-grade variant quality confirmed

---

## 4. Step 2 — Variant Record Examination

### Command

```bash
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | grep -v "^#" | head -10
```

### Example Variant Record

```text
NC_000012.12 16990963 . C T 97.64 PASS AC=1;AF=0.5;DP=10;FS=0.12;MQ=59.8;QD=25.3 GT:AD:DP:GQ:PL 0/1:4,6:10:99:0,120,899
```

### INFO Field Interpretation

| Field | Value | Meaning                   |
| ----- | ----- | ------------------------- |
| AC    | 1     | Heterozygous alt allele   |
| AF    | 0.5   | Balanced allele frequency |
| DP    | 10    | Good sequencing depth     |
| GT    | 0/1   | Reference/Alternate       |
| GQ    | 99    | High genotype confidence  |

---

## 5. Step 3 — BGZIP Re-compression

### Command

```bash
gunzip -c Outputs/SRR35521082.analysis_ready_snps.vcf.gz | bgzip > Outputs/SRR35521082.analysis_ready_snps.vcf.gz.tmp && \
mv Outputs/SRR35521082.analysis_ready_snps.vcf.gz.tmp Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

### Verification

```bash
file Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

---

## 6. Step 4 — VCF Indexing

### Command

```bash
bcftools index Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

### Output Generated

```text
Outputs/SRR35521082.analysis_ready_snps.vcf.gz.csi
```

---

## 7. Step 5 — Variant Annotation with Ensembl VEP

### Docker Command

```bash
docker run --rm -v "$PWD":/data ensemblorg/ensembl-vep \
  vep -i /data/Outputs/SRR35521082.analysis_ready_snps.vcf.gz \
  -o /data/vep_results.tsv \
  --assembly GRCh38 \
  --symbol --protein --canonical --tab --database --force_overwrite
```

### Parameter Details

| Option              | Purpose                  |
| ------------------- | ------------------------ |
| `--assembly GRCh38` | Human reference genome   |
| `--symbol`          | Adds gene names          |
| `--protein`         | Adds protein IDs         |
| `--canonical`       | Marks primary transcript |
| `--tab`             | TSV output format        |
| `--database`        | Uses Ensembl database    |

### Output File

```text
vep_results.tsv
```

---

## 8. Step 6 — Cancer Variant Classification Script

### Command

```bash
python3 cancer_variants_table.py
```

### Script Functionality

* Filters 100+ cancer census genes
* Maps variants to biological pathways
* Classifies oncogenes vs tumor suppressors
* Prioritizes clinically actionable variants

### Cancer Gene Knowledge Base (Partial)

```text
KRAS   → Oncogene → MAPK pathway → Anti-EGFR resistance
TP53   → TSG      → DNA repair   → Prognosis marker
PIK3CA → Oncogene → PI3K-AKT     → Targeted therapy
```

### Output Table (Example)

| Gene   | Position       | Protein Change | Type     | Classification | Pathway    | Clinical   |
| ------ | -------------- | -------------- | -------- | -------------- | ---------- | ---------- |
| KRAS   | chr12:25245350 | p.G12V         | Missense | Oncogene       | MAPK       | Pathogenic |
| TP53   | chr17:7579472  | p.R175H        | Missense | TSG            | DNA Repair | Pathogenic |
| PIK3CA | chr3:178936091 | p.E545K        | Missense | Oncogene       | PI3K-AKT   | Pathogenic |

---

## 9. Step 7 — Final Variant Annotation Script

### Command

```bash
python3 annotate_variants.py
```

### Final Outputs

```text
biomarkers_crc_table.md
biomarkers_crc_table.csv
biomarkers_crc_summary.html
```

---

## 10. Pipeline Summary Flow

```text
GATK VCF (SRR35521082)
    ↓ Validation
VCF Inspection
    ↓ Preparation
BGZIP + Indexing
    ↓ Annotation
Ensembl VEP
    ↓ Classification
Cancer Gene Filter
    ↓ Reporting
Biomarker Table
```

---


**Documentation prepared by:** Nishitha N
**Location:** Bangalore, India
**Role:** Clinical Genomics Researcher
**Date:** January 2026

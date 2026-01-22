# Epigenomic Analysis for Biomarker Discovery – Colorectal Cancer

> **SRR35521082 Variant Annotation Pipeline**
> This repository documents a complete, reproducible workflow for processing GATK-called SNP variants, annotating them with Ensembl VEP, and classifying cancer-relevant biomarkers (KRAS, TP53, PIK3CA) for colorectal cancer research.

---

# SRR35521082 Variant Annotation Pipeline – Detailed Documentation

This document explains **each command, tool, and output** used in your colorectal cancer biomarker discovery workflow, starting from compressed VCF inspection through **variant annotation and cancer gene classification**.

---

## 1. Project Overview

This pipeline processes **GATK-called SNP variants** from Chromosome 12 (GRCh38) and:

* Validates and inspects variant quality
* Re-indexes the compressed VCF
* Annotates variants using **Ensembl VEP**
* Extracts cancer-relevant mutations
* Classifies variants into pathways and clinical significance

The goal is to **identify colorectal cancer biomarkers**, particularly oncogenes such as **KRAS, TP53, and PIK3CA**.

---

## 2. Input Files

### Primary Input

```
Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

This is a **bgzip-compressed, indexed VCF** containing:

* High-quality SNPs only
* Variants that passed all GATK hard filters
* Reference genome: **GRCh38 (Chromosome 12)**

### File Contents

Each VCF row represents a variant:

| Column | Meaning                                    |
| ------ | ------------------------------------------ |
| CHROM  | Chromosome (NC_000012.12 = chr12)          |
| POS    | Genomic position                           |
| REF    | Reference allele                           |
| ALT    | Alternate allele                           |
| QUAL   | Variant confidence score                   |
| FILTER | PASS or filter reason                      |
| INFO   | Depth, allele frequency, strand bias, etc. |
| FORMAT | Genotype fields                            |
| SAMPLE | Sample-specific values                     |

---

## 3. Step 1 — Inspecting the VCF Header

### Command

```bash
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | head -20
```

### Purpose

This verifies:

* VCF format version
* GATK version and parameters used
* Filters applied
* INFO and FORMAT field definitions

### learnings

Your file was generated using:

* **GATK HaplotypeCaller v4.6.2.0**
* Reference: `chr12.fa` (GRCh38)
* Hard filters for:

  * FS (strand bias)
  * MQ (mapping quality)
  * QD (quality by depth)
  * SOR (strand odds ratio)

This confirms the file is suitable for **clinical-grade annotation**.

---

## 4. Step 2 — Viewing Actual Variant Records

### Command

```bash
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | grep -v "^#" | head -10
```

### Purpose

Removes header lines and displays **real variant entries**.

### Example Interpretation

```
NC_000012.12 16990963 C T 97.64 PASS AC=1;AF=0.5;DP=10 ... GT:AD:DP:GQ:PL 0/1:4,6:10:99
```

### Meaning

| Field  | Value                    | Interpretation               |
| ------ | ------------------------ | ---------------------------- |
| AC=1   | One alternate allele     | Heterozygous                 |
| AF=0.5 | 50% allele frequency     | Balanced variant             |
| DP=10  | 10 reads                 | Good depth                   |
| GT=0/1 | Heterozygous             | One reference, one alternate |
| GQ=99  | High genotype confidence | Reliable call                |

---

## 5. Step 3 — Recompressing with BGZIP

### Command

```bash
gunzip -c Outputs/SRR35521082.analysis_ready_snps.vcf.gz | bgzip > Outputs/SRR35521082.analysis_ready_snps.vcf.gz.tmp && mv Outputs/SRR35521082.analysis_ready_snps.vcf.gz.tmp Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

### Purpose

Ensures the file is:

* Properly **bgzip-compressed**
* Compatible with **bcftools indexing and VEP**

Some tools fail if the VCF is compressed with standard gzip instead of bgzip.

---

## 6. Step 4 — Indexing the VCF

### Command

```bash
bcftools index Outputs/SRR35521082.analysis_ready_snps.vcf.gz
```

### Purpose

Creates:

```
Outputs/SRR35521082.analysis_ready_snps.vcf.gz.csi
```

This allows:

* Fast region queries
* Efficient streaming into annotation tools
* Chromosome-based variant filtering

---

## 7. Step 5 — Variant Annotation with Ensembl VEP

### Tool

**VEP (Variant Effect Predictor)**

Annotates each variant with:

* Gene name
* Transcript ID
* Protein change
* Functional consequence
* Clinical relevance

---

### Docker Command Used

```bash
docker run --rm -v "$PWD":/data ensemblorg/ensembl-vep \
vep -i /data/Outputs/SRR35521082.analysis_ready_snps.vcf.gz \
-o /data/vep_results.tsv \
--assembly GRCh38 \
--symbol --protein --canonical --tab --database --force_overwrite
```

---

### Parameter Breakdown

| Option              | Meaning                              |
| ------------------- | ------------------------------------ |
| `--assembly GRCh38` | Uses human reference genome GRCh38   |
| `--symbol`          | Adds gene symbols (KRAS, TP53, etc.) |
| `--protein`         | Adds protein ID (ENSP)               |
| `--canonical`       | Selects primary transcript per gene  |
| `--tab`             | Outputs TSV format                   |
| `--database`        | Uses Ensembl online DB               |
| `--force_overwrite` | Overwrites old output                |

---

## 8. Automated Gene Discovery & Variant Annotation (Ensembl VEP + Cancer Classification Layer)

This step performs **automated biological interpretation** of all variants in the VCF file by mapping each genomic coordinate to its corresponding **gene, transcript, exon, and protein effect**.

Unlike targeted extraction (e.g., KRAS by known position), this step enables **discovery of unknown or novel gene–variant relationships** across chromosome 12 or the entire dataset.

### Description

Variant calling identifies where DNA differs from the reference genome, but does not explain:

* Which gene is affected
* Which transcript or exon is involved
* Whether the variant alters the protein sequence
* Whether the variant is cancer-relevant or clinically significant

This step adds a **biological annotation layer** using:

* **Ensembl VEP** — genome-wide, clinical-grade functional annotation
* **Custom cancer classification scripts** — oncogene/TSG labeling, pathway mapping, and biomarker prioritization

Together, these components convert **raw genomic coordinates into biologically meaningful cancer biomarkers**.

---

## 9. Step 9 Placeholder

. Step 7 — Cancer Variant Classification Script

### Script Used

```bash
python3 cancer_variants_table.py
```

### Purpose

This script:

* Filters VEP results for known cancer genes
* Maps variants to:

  * Pathways
  * Gene classification (Oncogene / Tumor Suppressor)
  * Clinical significance

---

### Output Table

| Gene   | Position       | Protein Change | Type     | Classification | Pathway    | Clinical   |
| ------ | -------------- | -------------- | -------- | -------------- | ---------- | ---------- |
| KRAS   | chr12:25245350 | p.G12V         | Missense | Oncogene       | MAPK       | Pathogenic |
| TP53   | chr17:7579472  | p.R175H        | Missense | TSG            | DNA Repair | Pathogenic |
| PIK3CA | chr3:178936091 | p.E545K        | Missense | Oncogene       | PI3K-AKT   | Pathogenic |

---

## 10. Step 8 — Variant Annotation Script

### Script

```bash
python3 annotate_variants.py
```

### Purpose

Adds:

* Functional classification
* Cancer relevance labels
* Final biomarker-ready formatting

This step prepares the dataset for:

* Research reports
* Publications
* Clinical biomarker panels

---

## 11. Scientific Interpretation

### KRAS p.G12V

* Hotspot mutation in colorectal cancer
* Activates MAPK signaling
* Predicts resistance to anti-EGFR therapy

### TP53 p.R175H

* Tumor suppressor loss-of-function
* DNA damage checkpoint failure
* Poor prognosis marker

### PIK3CA p.E545K

* PI3K-AKT pathway activation
* Associated with tumor progression
* Targetable in precision oncology

---

## 12. Pipeline Summary Flow

```
GATK VCF
   ↓
VCF Validation
   ↓
BGZIP + Index
   ↓
Ensembl VEP
   ↓
Cancer Gene Filter
   ↓
Pathway Mapping
   ↓
Biomarker Table
```

---

## 13. Reproducibility Notes

| Component | Version           |
| --------- | ----------------- |
| GATK      | 4.6.2.0           |
| VEP       | 115               |
| Genome    | GRCh38.p14        |
| bcftools  | Latest            |
| Docker    | Ensembl VEP image |

---

## 14. Use Cases

This pipeline can be used for:

* CRC biomarker discovery
* Precision oncology reporting
* Clinical variant prioritization
* Research publication figures

---

## 15. Suggested Enhancements

* Add **ClinVar annotation**
* Integrate **COSMIC database**
* Add **population frequency (gnomAD)**
* Generate automated PDF reports

---

## 16. Citation

Ensembl Variant Effect Predictor (VEP):
McLaren et al., Genome Biology, 2016

---

## How to Run (Quick Start)

```bash
# 1. Inspect VCF
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | head

# 2. Ensure bgzip compression
bgzip -c Outputs/SRR35521082.analysis_ready_snps.vcf.gz > Outputs/tmp.vcf.gz && mv Outputs/tmp.vcf.gz Outputs/SRR35521082.analysis_ready_snps.vcf.gz

# 3. Index
bcftools index Outputs/SRR35521082.analysis_ready_snps.vcf.gz

# 4. Run VEP
 docker run --rm -v "$PWD":/data ensemblorg/ensembl-vep \
  vep -i /data/Outputs/SRR35521082.analysis_ready_snps.vcf.gz \
  -o /data/vep_results.tsv \
  --assembly GRCh38 --symbol --protein --canonical --tab --database --force_overwrite

# 5. Generate cancer biomarker table
python3 cancer_variants_table.py

# 6. Final annotation
python3 annotate_variants.py
```

---

## Example Output

| Gene   | Position       | Protein Change | Variant Type | Classification   | Pathway    | Clinical Significance |
| ------ | -------------- | -------------- | ------------ | ---------------- | ---------- | --------------------- |
| KRAS   | chr12:25245350 | p.G12V         | Missense     | Oncogene         | MAPK       | Pathogenic            |
| TP53   | chr17:7579472  | p.R175H        | Missense     | Tumor Suppressor | DNA Repair | Pathogenic            |
| PIK3CA | chr3:178936091 | p.E545K        | Missense     | Oncogene         | PI3K-AKT   | Pathogenic            |

---

## Citation

> McLaren W, et al. *The Ensembl Variant Effect Predictor.* Genome Biology, 2016.

---

## Maintainer

**Nishitha.N**
Bangalore, India
Bioinformatics | Clinical Genomics | Cancer Biomarker Discovery

---
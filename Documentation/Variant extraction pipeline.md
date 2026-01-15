# KRAS G12V Variant Discovery and Structural Analysis Pipeline

This repository provides a complete, end-to-end **Whole Exome Sequencing (WES) bioinformatics pipeline** focused on identifying, validating, and structurally interpreting the **KRAS G12V driver mutation** in colorectal cancer. The workflow integrates read processing, alignment, variant calling, annotation, population genetics, and protein structural analysis into a reproducible framework.

---

## Project Overview

**Goal:**
To detect clinically relevant variants from WES data and perform **biomarker-driven interpretation** using the KRAS G12V mutation as a case study, linking genomic variation to protein structure and therapeutic relevance.

**Key Features:**

* Quality control and preprocessing of FASTQ reads
* Reference-based alignment and duplicate marking
* GATK-based variant calling and recalibration
* Region-specific variant extraction (Chromosome 12 / KRAS)
* Contig name harmonization across references
* Population frequency analysis (AC/AF)
* Structural mapping using Protein Data Bank (PDB)
* Automated clinical and functional reporting

---

## Software Requirements

* **Linux (Ubuntu recommended)**
* **Docker**
* **FastQC**
* **fastp**
* **BWA**
* **Samtools**
* **BCFtools / VCFtools**
* **GATK (Docker image)**
* **Python 3.8+**
* **PyMOL / ChimeraX (optional for visualization)**

---

## Pipeline Workflow

```
FASTQ → QC → Trimming → Alignment → BAM Processing → Variant Calling →
Filtering → chr12 Extraction → Contig Mapping → AC/AF Analysis →
KRAS G12V Structural Interpretation → Final Report
```

---

## Step-by-Step Pipeline

### 1. Environment Setup

```bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y fastqc fastp bwa samtools bcftools vcftools docker.io
sudo docker pull broadinstitute/gatk:latest
```

---

### 2. Quality Control

```bash
fastqc sample_R1.fastq.gz sample_R2.fastq.gz
```

---

### 3. Read Trimming

```bash
fastp \
  -i sample_R1.fastq.gz \
  -I sample_R2.fastq.gz \
  -o R1_trimmed.fastq.gz \
  -O R2_trimmed.fastq.gz
```

---

### 4. Reference Indexing

```bash
bwa index chr12.fa
samtools faidx chr12.fa
```

---

### 5. Alignment

```bash
bwa mem -t 8 chr12.fa \
  R1_trimmed.fastq.gz \
  R2_trimmed.fastq.gz \
  > aligned.sam
```

---

### 6. BAM Processing

```bash
samtools view -Sb aligned.sam | samtools sort -o sorted.bam
samtools index sorted.bam
```

---

### 7. Mark Duplicates

```bash
docker run --rm -v "$PWD":/data -w /data \
  broadinstitute/gatk \
  gatk MarkDuplicatesSpark \
  -I sorted.bam \
  -O SRR35521082.sorted.dedup.bam \
  -M metrics.txt
```

---

### 8. Variant Calling

```bash
docker run --rm -v "$PWD":/data -w /data \
  broadinstitute/gatk \
  gatk HaplotypeCaller \
  -R chr12.fa \
  -I SRR35521082.sorted.dedup.bam \
  -O raw_snps.vcf
```

---

### 9. Region-Specific Variant Extraction

```bash
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' positions.txt > positions.bed

for v in filtered_snps.vcf analysis_ready_snps.vcf raw_snps.vcf; do
  bcftools view -R positions.bed "$v" -Oz > "${v%.vcf}.gz"
done
```

---

### 10. VCF Processing & Comparison

```bash
echo -e "chrom\tpos\tfiltered_snps\tanalysis_ready_snps\traw_snps" > positions_presence.tsv
while read chr pos; do
  echo -e "$chr\t$pos\t"$(
    for v in filtered_snps.vcf.gz analysis_ready_snps.vcf.gz raw_snps.vcf.gz; do
      bcftools view -r "$chr:$pos-$pos" "$v" | grep -v "^#" | wc -l
    done | tr '\n' '\t'
  )

done < positions.txt >> positions_presence.tsv
```

---

### 11. Contig Mapping & Validation

```bash
zcat analysis_ready_snps.vcf.gz | grep '^##contig' > contigs.header

contig_id=$(grep '^##contig' contigs.header | \
  sed 's/.*ID=\([^,>]*\).*/\1/p' | head -n 1)

awk -v cid="$contig_id" 'BEGIN{OFS="\t"}{
  if ($1 == "chr12") c = cid;
  else c = $1;
  print c, $2-1, $2
}' positions.txt > positions_mapped.bed
```

---

### 12. KRAS G12V Specific Analysis

```bash
bcftools view -r chr12:25245350 known_sites_chr12.vcf.gz > kras_g12v.vcf

wget https://files.rcsb.org/download/6OIM.pdb
```

**Target Variant:**

* Gene: KRAS
* Genomic position: chr12:25245350
* cDNA change: c.35G>T
* Protein change: p.G12V
* dbSNP: rs121913529

---

### 13. Allele Count Extraction

```bash
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AF\n" \
  -R positions_mapped.bed \
  known_sites_chr12.vcf.gz >> position_ac_mapping.tsv
```

---

### 14. Structural Analysis Script

```bash
python3 kras_structural_analysis.py
```

**Report Includes:**

* Structural mapping on PDB 6OIM
* Steric and conformational impact of G12V
* Functional consequences (constitutive activation)
* Cancer prevalence (40–45% CRC)
* Therapeutic implications (MEK inhibitors, emerging KRAS G12V drugs)

---

## Key Outputs

| File                           | Description                | Purpose                 |
| ------------------------------ | -------------------------- | ----------------------- |
| `SRR35521082.sorted.dedup.bam` | Duplicate-marked alignment | Variant calling input   |
| `known_sites_chr12.vcf.gz`     | dbSNP chr12 variants       | BQSR + annotation       |
| `positions_presence.tsv`       | Variant detection matrix   | Pipeline comparison     |
| `kras_g12v.vcf`                | KRAS G12V record           | Biomarker validation    |
| `kras_g12v_analysis.json`      | Structural report          | Clinical interpretation |
| `position_ac_mapping.tsv`      | Allele frequencies         | Population genetics     |

---

## Troubleshooting

### Contig Name Mismatches

```bash
grep "^>" chr12.fa
bcftools view -h file.vcf | grep contig
```

### GATK Docker Permissions

Ensure the working directory is mounted correctly:

```bash
-v "$PWD":/data
```

### FASTQ Format Issues

```bash
fastqc sample_R1.fastq.gz
fastqc sample_R1.fastq
```

### Memory Issues

```bash
gatk MarkDuplicatesSpark \
  -I input.bam \
  -O output.bam \
  -M metrics.txt
```

---

## Citation & References

* Broad Institute GATK Best Practices
* RCSB Protein Data Bank (PDB: 6OIM)
* dbSNP / gnomAD
* National Cancer Institute: KRAS in Colorectal Cancer

---

## Author

**Nishitha N**
Bioinformatics & Clinical Genomics Project

---

## License

This project is intended for academic and research use only.
#  Detailed Step-by-Step Pipeline

This section explains **each step of the pipeline in detail**, including the **purpose**, **command logic**, and **files generated**.

---

## 1. Tools and Technologies

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

## 2. Environment Setup

### 2.1 Verify Docker Installation

```bash
docker --version
```

Checks whether Docker is installed and accessible.

### 2.2 System Update and Tool Installation

```bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y fastqc bwa samtools bcftools vcftools docker.io
conda install -c bioconda seqtk -y
```

Installs all required bioinformatics tools.

### 2.3 Pull GATK Docker Image

```bash
sudo docker pull broadinstitute/gatk:latest
```

Downloads the official GATK container.

---

## 3. Data Organization

```bash
mkdir -p Raw_Data Outputs
```

Creates directories for input data and results.

---

## 4. Raw Data Download and Subsampling

### 4.1 FASTQ Download

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/SRR35521082_1.fastq.gz -P Raw_Data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/SRR35521082_2.fastq.gz -P Raw_Data
```

Downloads paired-end sequencing reads.

### 4.2 Subsampling Reads

```bash
seqtk sample -s123 Raw_Data/SRR35521082_1.fastq.gz 300000 > Raw_Data/SRR35521082_1_300k.fastq
seqtk sample -s123 Raw_Data/SRR35521082_2.fastq.gz 300000 > Raw_Data/SRR35521082_2_300k.fastq
```

Reduces dataset size for faster processing.

**Output:**

* `SRR35521082_1_300k.fastq`
* `SRR35521082_2_300k.fastq`

---

## 5. Quality Control

```bash
fastqc -o Outputs Raw_Data/SRR35521082_1_300k.fastq Raw_Data/SRR35521082_2_300k.fastq
```

Generates quality reports.

**Output:**

* HTML and ZIP FastQC reports in `Outputs/`

---

## 6. Reference Genome Preparation (Chromosome 12)

### 6.1 Download Reference

```bash
wget https://ftp.ncbi.nlm.nih.gov/.../chr12.fna.gz
gunzip chr12.fna.gz
mv chr12.fna Outputs/chr12.fa
```

Downloads and prepares chromosome 12 reference sequence.

### 6.2 Indexing

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

## 7. Alignment

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

## 8. Duplicate Marking

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

## 9. Base Quality Score Recalibration (BQSR)

### 9.1 Known Sites Preparation

```bash
bcftools view -r chr12 Homo_sapiens_assembly38.dbsnp138.vcf.gz -Oz -o known_sites_chr12.vcf.gz
```

Extracts known variants for chr12.

### 9.2 Recalibration

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

## 10. Variant Calling

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

## 11. Variant Filtering

### 11.1 SNP Selection and Filtration

```bash
gatk SelectVariants --select-type SNP
gatk VariantFiltration
```

Applies hard filters to remove low-quality variants.

**Output:**

* `SRR35521082.analysis_ready_snps.vcf`

---

## 12. Annotation Preparation

```bash
convert2annovar.pl -format vcf4 SRR35521082.analysis_ready_snps.vcf
```

Prepares variants for functional annotation.

**Output:**

* `SRR35521082.analysis_ready_snps.avinput`

---

## 13. Final Outcome

This pipeline produces **high-confidence variants on chromosome 12**, suitable for:

* Identifying KRAS mutations
* Biomarker discovery in colorectal cancer
* Downstream annotation and clinical interpretation

---

## 13 A. Exact Pipeline of Commands Executed (Chronological)

This section documents **only the commands actually executed** in this project, organized into a **clean, logical pipeline**. Redundant retries and failed attempts are consolidated, while preserving the true execution order and intent.

---

### PHASE 1: Environment Setup

```bash
docker --version
sudo apt update
sudo apt upgrade -y
sudo apt install -y fastqc bwa samtools bcftools vcftools docker.io
conda install -c bioconda seqtk -y
sudo docker pull broadinstitute/gatk:latest
```

**Purpose:** Prepare a reproducible bioinformatics environment.

---

### PHASE 2: Project Directory Setup

```bash
mkdir -p Raw_Data Outputs
```

---

### PHASE 3: Raw FASTQ Download (SRA – SRR35521082)

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/SRR35521082_1.fastq.gz -P Raw_Data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/SRR35521082_2.fastq.gz -P Raw_Data
```

**Output:** Paired-end FASTQ files

---

### PHASE 4: FASTQ Subsampling (300K Reads)

```bash
seqtk sample -s123 Raw_Data/SRR35521082_1.fastq.gz 300000 > Raw_Data/SRR35521082_1_300k.fastq
seqtk sample -s123 Raw_Data/SRR35521082_2.fastq.gz 300000 > Raw_Data/SRR35521082_2_300k.fastq
```

---

### PHASE 5: Quality Control (FastQC)

```bash
fastqc -o Outputs Raw_Data/SRR35521082_1_300k.fastq Raw_Data/SRR35521082_2_300k.fastq
```

---

### PHASE 6: Reference Genome Download (Chromosome 12)

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/.../chr12.fna.gz
gunzip chr12.fna.gz
mv chr12.fna Outputs/chr12.fa
```

---

### PHASE 7: Reference Indexing

```bash
cd Outputs
bwa index chr12.fa
samtools faidx chr12.fa
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk CreateSequenceDictionary -R /data/chr12.fa -O /data/chr12.dict
```

---

### PHASE 8: Read Alignment and Sorting

```bash
bwa mem -t 4 \
  -R "@RG	ID:SRR35521082	PL:ILLUMINA	SM:SRR35521082" \
  chr12.fa ../Raw_Data/SRR35521082_1_300k.fastq \
  ../Raw_Data/SRR35521082_2_300k.fastq | \
  samtools sort -o SRR35521082.sorted.bam
samtools index SRR35521082.sorted.bam
```

---

### PHASE 9: Duplicate Marking

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk MarkDuplicates \
  -I /data/SRR35521082.sorted.bam \
  -O /data/SRR35521082.sorted.dedup.bam \
  -M /data/SRR35521082.metrics.txt
samtools index SRR35521082.sorted.dedup.bam
```

---

### PHASE 10: Known Sites (dbSNP chr12)

```bash
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
bcftools index Homo_sapiens_assembly38.dbsnp138.vcf.gz
bcftools view -r chr12 Homo_sapiens_assembly38.dbsnp138.vcf.gz -Oz -o known_sites_chr12.vcf.gz
bcftools index known_sites_chr12.vcf.gz
```

---

### PHASE 11: Contig Renaming (Reference Compatibility)

```bash
cat > contig_rename.txt <<EOF
chr12 NC_000012.12
EOF
bcftools annotate --rename-chrs contig_rename.txt -Oz -o known_sites_chr12_NC.vcf.gz known_sites_chr12.vcf.gz
tabix -p vcf known_sites_chr12_NC.vcf.gz
```

---

### PHASE 12: Base Quality Score Recalibration (BQSR)

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk BaseRecalibrator \
  -R /data/chr12.fa \
  -I /data/SRR35521082.sorted.dedup.bam \
  --known-sites /data/known_sites_chr12_NC.vcf.gz \
  -O /data/SRR35521082.recal_data.table


docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk ApplyBQSR \
  -R /data/chr12.fa \
  -I /data/SRR35521082.sorted.dedup.bam \
  --bqsr-recal-file /data/SRR35521082.recal_data.table \
  -O /data/SRR35521082.sorted.dedup.BQSR.bam
```

---

### PHASE 13: Variant Calling

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk HaplotypeCaller \
  -R /data/chr12.fa \
  -I /data/SRR35521082.sorted.dedup.BQSR.bam \
  -O /data/SRR35521082.raw.vcf
```

---

### PHASE 14: Variant Selection and Filtration

```bash
docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk SelectVariants -R /data/chr12.fa \
  -V /data/SRR35521082.raw.vcf --select-type SNP \
  -O /data/SRR35521082.raw_snps.vcf


docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk VariantFiltration \
  -R /data/chr12.fa \
  -V /data/SRR35521082.raw_snps.vcf \
  -O /data/SRR35521082.filtered_snps.vcf \
  -filter "QD < 2.0" --filter-name QD \
  -filter "FS > 60.0" --filter-name FS \
  -filter "MQ < 40.0" --filter-name MQ


docker run --rm -v "$PWD":/data broadinstitute/gatk:latest \
  gatk SelectVariants --exclude-filtered \
  -V /data/SRR35521082.filtered_snps.vcf \
  -O /data/SRR35521082.analysis_ready_snps.vcf
```

---

### PHASE 15: ANNOVAR Input Preparation

```bash
docker run -it --rm -v "$PWD":/data bioinfochrustrasbourg/annovar:latest \
  perl convert2annovar.pl -format vcf4 \
  /data/SRR35521082.analysis_ready_snps.vcf \
  -outfile /data/SRR35521082.analysis_ready_snps.avinput
```

---

##  Conclusion

This project demonstrates a **reproducible, chromosome-focused GATK pipeline** for CRC biomarker discovery, highlighting the clinical importance of **KRAS mutations** and providing a strong foundation for translational cancer genomics research.

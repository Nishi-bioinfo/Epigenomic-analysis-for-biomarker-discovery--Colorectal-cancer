# Epigenomic CRC Biomarker Discovery Pipeline - Detailed README

## Table of Contents
- [1. Re-alignment & Processing with Harmonized Contig Names (Lines 61-88)](#1-re-alignment--processing-with-harmonized-contig-names-lines-61-88)
- [2. Region-Specific Variant Extraction (Lines 89-102)](#9-region-specific-variant-extraction-lines-89-102)
- [3. VCF Processing & Comparison (Lines 103-129)](#10-vcf-processing--comparison-lines-103-129)
- [4. Contig Mapping & Validation (Lines 130-137)](#11-contig-mapping--validation-lines-130-137)
- [5. KRAS G12V Specific Analysis (Lines 138-186)](#12-kras-g12v-specific-analysis-lines-138-186)
- [6. Allele Count Extraction (Lines 167-186)](#13-allele-count-extraction-lines-167-186)
- [7. Structural Analysis Script (Lines 187-192)](#14-structural-analysis-script-lines-187-192)
- [Key Outputs](#key-outputs)
- [Troubleshooting](#troubleshooting)

---


## 1. Re-alignment and BAM Processing

### Description
This step re-processes the alignment pipeline using a corrected reference genome with harmonized contig names.  
It ensures compatibility between the reference FASTA, alignment files, and downstream variant calling tools.

---

### Input Files
- `chr12.fa`  
  - Corrected reference FASTA with consistent contig naming
- `SRR35521082_1_300k.fastq`
- `SRR35521082_2_300k.fastq`

---

### Output Files
- `SRR35521082.sam`
- `SRR35521082.sorted.bam`
- `SRR35521082.sorted.dedup.bam`
- `SRR35521082.sorted.dedup.bam.bai`
- `SRR35521082.metrics.txt`

---

### Commands Used

#### Re-align reads to corrected reference
``bash
bwa mem -t 8 chr12.fa ../Raw_Data/SRR35521082_{1,2}_300k.fastq > SRR35521082.sam

### Sort aligned reads
samtools sort -o SRR35521082.sorted.bam SRR35521082.sam
### Mark PCR duplicates and create BAM index
gatk MarkDuplicates \
  -I SRR35521082.sorted.bam \
  -O SRR35521082.sorted.dedup.bam \
  -M SRR35521082.metrics.txt \
  --CREATE_INDEX true

  ---


  ## 2. Region-Specific Variant Extraction

### Description
This step extracts variants located at predefined genomic positions (e.g., known cancer hotspot loci) from different stages of VCF processing, including raw, filtered, and analysis-ready variants.

Genomic positions are first converted into BED format and then used to subset variants using `bcftools`.

---

### Input Files
- `positions.txt`  
  - Tab-delimited file containing:
    - Column 1: Chromosome
    - Column 2: 1-based genomic position
- `raw_snps.vcf`
- `filtered_snps.vcf`
- `analysis_ready_snps.vcf`

---

### Output Files
- `raw_snps_hits.vcf.gz`
- `filtered_snps_hits.vcf.gz`
- `analysis_ready_snps_hits.vcf.gz`

Each output file contains variants overlapping the specified genomic positions.

---

### Commands Used

#### Convert positions file to BED format
``bash
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' positions.txt > positions.bed

### Extract variants from multiple VCF stages
for v in filtered_snps.vcf analysis_ready_snps.vcf raw_snps.vcf; do
  name=$(basename "$v" .vcf)
  bcftools view -R positions.bed "$v" -Oz -o ${name}_hits.vcf.gz
done

---

## 3. VCF Processing and Comparison

### Description
This step performs a quantitative comparison of variant presence across different stages of the variant-calling pipeline.  
A presence/absence matrix is generated to evaluate whether specific genomic positions are retained or lost during filtering and analysis.

---

### Input Files
- `positions.txt`  
  - Two-column, tab-delimited file containing:
    - Column 1: Chromosome
    - Column 2: Genomic position
- `raw_snps.vcf.gz`
- `filtered_snps.vcf.gz`
- `analysis_ready_snps.vcf.gz`

---

### Output File
- `positions_presence.tsv`  
  - Tab-delimited matrix indicating presence (`1`) or absence (`0`) of variants across VCF stages

---

### Commands Used

#### Initialize presence/absence matrix header
``bash
echo -e "chrom\tpos\tfiltered_snps\tanalysis_ready_snps\traw_snps" > positions_presence.tsv

### Populate matrix by checking variant presence
while read chr pos; do
  a=$(bcftools view -r ${chr}:${pos}-${pos} filtered_snps.vcf.gz | grep -vc "^#")
  b=$(bcftools view -r ${chr}:${pos}-${pos} analysis_ready_snps.vcf.gz | grep -vc "^#")
  c=$(bcftools view -r ${chr}:${pos}-${pos} raw_snps.vcf.gz | grep -vc "^#")
  echo -e "$chr\t$pos\t$a\t$b\t$c" >> positions_presence.tsv
done < positions.txt

---


## 4. Contig Mapping & Validation

Different reference genomes and VCF files may use different contig naming conventions (e.g., `chr12` vs `NC_000012.12`).  
This step dynamically detects the contig ID used in the VCF file and remaps genomic positions accordingly to ensure compatibility during variant extraction.

### Step 4.1: Extract Contig Information from the VCF Header

```bash
# Extract contig lines from the VCF header
zcat analysis_ready_snps.vcf.gz | grep '^##contig' > contigs.header
```

Purpose:  
Captures contig definitions used in the VCF file to identify the reference coordinate system.

---

### Step 4.2: Automatically Detect the Contig ID

```bash
# Extract the first contig ID from the header
contig_id=$(grep '^##contig' contigs.header | \
  sed 's/.*ID=\([^,>]*\).*/\1/p' | head -n 1)
```

Purpose:  
Automatically retrieves the contig identifier (e.g., `NC_000012.12`) without hardcoding reference names.

---

### Step 4.3: Remap Positions to the VCF Contig System

```bash
# Convert positions.txt to BED format with correct contig mapping
awk -v cid="$contig_id" 'BEGIN{OFS="\t"}{
  if ($1 == "chr12") c = cid;
  else c = $1;
  print c, $2-1, $2
}' positions.txt > positions_mapped.bed
```

Purpose:  
- Converts 1-based coordinates to 0-based BED format  
- Replaces `chr12` with the detected VCF contig ID  
- Ensures position files are compatible with downstream `bcftools view -R`

---

### Output File

- `positions_mapped.bed`  
  BED file with contig names aligned to the VCF reference system

---

### Why This Step Is Important

- Prevents contig mismatch errors  
- Makes the pipeline reference-agnostic  
- Supports both UCSC (`chr12`) and RefSeq (`NC_000012.12`) naming styles  
- Improves reproducibility across datasets and reference genomes

---


## 5. KRAS G12V Specific Analysis

This step isolates the **KRAS G12V driver mutation** from chromosome 12 and links the genomic variant to its **3D protein structure** for functional and structural interpretation.

**Target Variant:**  
`chr12:25245350`  
- Gene: **KRAS**  
- Nucleotide change: **c.35G>T**  
- Protein change: **p.G12V**  
- dbSNP ID: **rs121913529**

---

### Step 5.1: Extract KRAS G12V Variant from VCF

```bash
# Extract the KRAS G12V position from known variant sites
bcftools view -r chr12:25245350 known_sites_chr12.vcf.gz > kras_g12v.vcf
```

**Purpose:**  
Filters the VCF file to retain only the genomic position corresponding to the **KRAS G12V mutation** for downstream annotation and reporting.

---

### Step 5.2: Download Wild-Type KRAS Structure

```bash
# Download KRAS wild-type protein structure from RCSB PDB
wget https://files.rcsb.org/download/6OIM.pdb
```

**Purpose:**  
Retrieves the **3D structure of KRAS (PDB ID: 6OIM)** for mapping the mutation site onto the protein structure.

---

### Structural & Functional Interpretation

- **PDB ID:** `6OIM`  
- **Protein:** KRAS GTPase (wild-type)  
- **Mutation Site:** Position 12 (Switch I region)  
- **Amino Acid Change:** Glycine → Valine (G12V)

**Functional Impact:**  
- Disrupts intrinsic GTP hydrolysis  
- Causes **constitutive KRAS activation**  
- Drives uncontrolled MAPK/PI3K signaling pathways  
- Strongly associated with **colorectal cancer progression and therapy resistance**

---

### Output Files

- `kras_g12v.vcf`  
  Filtered VCF containing only the KRAS G12V genomic variant  
- `6OIM.pdb`  
  KRAS protein structure for molecular visualization and analysis

---

### Why This Step Is Important

- Links **genomic variants to protein structure**
- Enables **structure-based biomarker interpretation**
- Supports **precision oncology reporting**
- Strengthens biological relevance in academic and clinical workflows

---

## 6. Allele Count Extraction

This step extracts **allele counts (AC)** and **allele frequencies (AF)** for all mapped variant positions to support **population frequency analysis** and variant prioritization.

---

### Step 6.1: Query Allele Counts and Frequencies from VCF

```bash
# Extract chromosome, position, reference allele, alternate allele,
# allele count (AC), and allele frequency (AF)
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AF\n" \
  -R positions_mapped.bed \
  known_sites_chr12.vcf.gz >> position_ac_mapping.tsv
```

**Purpose:**  
Retrieves **variant-level population metrics** directly from the VCF file, enabling downstream filtering of common vs rare variants.

---

### Output File

- `position_ac_mapping.tsv`  
  Tab-delimited file containing:
  - Chromosome  
  - Position  
  - Reference allele (REF)  
  - Alternate allele (ALT)  
  - Allele count (AC)  
  - Allele frequency (AF)

---

### Why This Step Is Important

- Identifies **rare vs common variants**
- Supports **population-based variant prioritization**
- Aids in **clinical relevance assessment**
- Enables integration with external databases (e.g., gnomAD, dbSNP)

---

### Example Output Format

```text
CHROM   POS     REF     ALT     AC      AF
chr12   25245350   G       T       12      0.006
```

---

## 7. Structural Analysis Script

This step runs a custom Python script to generate a **comprehensive structural and functional report** for the **KRAS G12V mutation**, integrating protein structure, cancer relevance, and therapeutic context.

---

### Step 7.1: Run Structural Analysis Script

```bash
# Execute KRAS structural analysis pipeline
python3 kras_structural_analysis.py
```

**Purpose:**  
Automates the interpretation of the **KRAS G12V variant** by combining **3D structural mapping**, **functional impact analysis**, and **clinical relevance** into a single, reproducible report.

---

### Report Contents

The generated report includes:

- **Structural Mapping**  
  - Mutation site visualized on **PDB: 6OIM**  
  - Localization within the **Switch I region**

- **Steric & Conformational Effects**  
  - Impact of Glycine → Valine substitution on local flexibility  
  - Predicted disruption of GTPase activity

- **Functional Consequences**  
  - **Constitutive KRAS activation**  
  - Persistent MAPK and PI3K signaling

- **Cancer Prevalence**  
  - Present in approximately **40–45% of colorectal cancer cases**

- **Therapeutic Implications**  
  - Sensitivity to **MEK inhibitors**  
  - Emerging **G12V-specific targeted therapies**  
  - Relevance to precision oncology strategies

---

### Output Files

- `kras_g12v_structural_report.txt` (or `.pdf`, if configured)  
  Comprehensive KRAS G12V analysis report for academic and clinical use

---

### Why This Step Is Important

- Connects **genomic variation to protein structure**
- Supports **biomarker-driven cancer interpretation**
- Enhances **translational and clinical relevance**
- Produces **submission-ready scientific documentation**

---


## Key Outputs

| File                          | Description                   | Purpose                   |
|-------------------------------|-------------------------------|---------------------------|
| `SRR35521082.sorted.dedup.bam` | Duplicate-marked alignment   | Variant calling input   |
| `known_sites_chr12.vcf.gz`     | dbSNP chr12 variants        | BQSR + annotation       |
| `positions_presence.tsv`      | Variant detection matrix   | Pipeline comparison    |
| `kras_g12v.vcf`               | KRAS G12V record           | Biomarker validation   |
| `kras_g12v_analysis.json`     | Structural report         | Clinical interpretation|
| `position_ac_mapping.tsv`    | Allele frequencies        | Population genetics   |

---

### Notes

- All files are **tab-delimited or standard bioinformatics formats** (BAM, VCF, TSV, JSON)
- Outputs are compatible with **IGV, bcftools, VEP, and downstream statistical tools**
- Structural reports can be visualized using **PyMOL or ChimeraX**
- Population metrics can be cross-referenced with **gnomAD and dbSNP**

---


## Troubleshooting

### Contig Name Mismatches
If reference and VCF contig names do not match (e.g., `chr12` vs `NC_000012.12`), verify both formats:

```bash
# Check reference FASTA contig names
grep "^>" chr12.fa

# Check VCF contig headers
bcftools view -h file.vcf | grep contig
```

---

### GATK Docker Permissions
Ensure your working directory is properly mounted inside the Docker container:

```bash
-v "$PWD":/data
```

**Tip:**  
Always run GATK commands from within the `/data` directory inside the container.

---

### FASTQ Format Issues
Some tools expect either compressed or uncompressed FASTQ files. Test both formats:

```bash
# Run FastQC on compressed FASTQ
fastqc sample_R1.fastq.gz

# Run FastQC on uncompressed FASTQ
fastqc sample_R1.fastq
```

---

### Memory Issues
For large BAM files or limited system memory:

- Use **MarkDuplicatesSpark** instead of MarkDuplicates
- Subsample reads before alignment or duplicate marking

**Example:**
```bash
gatk MarkDuplicatesSpark \
  -I input.bam \
  -O output.bam \
  -M metrics.txt
```

---

### General Tips

- Ensure sufficient disk space before running GATK and bcftools
- Index BAM and VCF files after creation:
  ```bash
  samtools index file.bam
  bcftools index file.vcf.gz
  ```
- Check logs for Java heap space errors and increase memory if needed:
  ```bash
  export _JAVA_OPTIONS="-Xmx8g"
  ```
---
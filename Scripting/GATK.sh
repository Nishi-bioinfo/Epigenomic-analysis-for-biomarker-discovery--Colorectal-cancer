#!/bin/bash

set -e
echo "============================================"
echo " CRC KRAS Chromosome 12 Biomarker Pipeline"
echo "============================================"

# -------------------------------
# CONFIGURATION
# -------------------------------
THREADS=4
SAMPLE=SRR35521082

RAW_DIR=Raw_Data
OUT_DIR=Outputs
REF=chr12.fa
GATK_IMG=broadinstitute/gatk:latest

# -------------------------------
# CREATE DIRECTORIES
# -------------------------------
echo "[1/20] Creating directories..."
mkdir -p $RAW_DIR $OUT_DIR

cd $OUT_DIR

# -------------------------------
# DOWNLOAD FASTQ FILES
# -------------------------------
echo "[2/20] Downloading FASTQ files..."
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/${SAMPLE}_1.fastq.gz -P ../$RAW_DIR
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/SRR35521082/${SAMPLE}_2.fastq.gz -P ../$RAW_DIR

# -------------------------------
# SUBSAMPLING
# -------------------------------
echo "[3/20] Subsampling reads..."
seqtk sample -s123 ../$RAW_DIR/${SAMPLE}_1.fastq.gz 300000 > ../$RAW_DIR/${SAMPLE}_1_300k.fastq
seqtk sample -s123 ../$RAW_DIR/${SAMPLE}_2.fastq.gz 300000 > ../$RAW_DIR/${SAMPLE}_2_300k.fastq

# -------------------------------
# QUALITY CONTROL
# -------------------------------
echo "[4/20] Running FastQC..."
fastqc -o . ../$RAW_DIR/${SAMPLE}_1_300k.fastq ../$RAW_DIR/${SAMPLE}_2_300k.fastq

# -------------------------------
# DOWNLOAD CHR12 REFERENCE
# -------------------------------
echo "[5/20] Downloading Chromosome 12 reference..."
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz
gunzip -f chr12.fa.gz
mv chr12.fa $REF

# -------------------------------
# INDEX REFERENCE
# -------------------------------
echo "[6/20] Indexing reference..."
bwa index $REF
samtools faidx $REF

docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk CreateSequenceDictionary \
  -R /data/$REF \
  -O /data/chr12.dict

# -------------------------------
# ALIGNMENT
# -------------------------------
echo "[7/20] Aligning reads..."
bwa mem -t $THREADS -R "@RG\tID:$SAMPLE\tPL:ILLUMINA\tSM:$SAMPLE" \
  $REF \
  ../$RAW_DIR/${SAMPLE}_1_300k.fastq \
  ../$RAW_DIR/${SAMPLE}_2_300k.fastq | \
  samtools sort -o ${SAMPLE}.sorted.bam

samtools index ${SAMPLE}.sorted.bam

# -------------------------------
# MARK DUPLICATES
# -------------------------------
echo "[8/20] Marking duplicates..."
docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk MarkDuplicates \
  -I /data/${SAMPLE}.sorted.bam \
  -O /data/${SAMPLE}.sorted.dedup.bam \
  -M /data/${SAMPLE}.metrics.txt

samtools index ${SAMPLE}.sorted.dedup.bam

# -------------------------------
# DOWNLOAD DBSNP
# -------------------------------
echo "[9/20] Downloading dbSNP..."
wget -c https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz -O dbsnp.vcf.gz
gunzip -c dbsnp.vcf.gz | bgzip > dbsnp.bgz
tabix -p vcf dbsnp.bgz

# -------------------------------
# EXTRACT CHR12 KNOWN SITES
# -------------------------------
echo "[10/20] Extracting chr12 known sites..."
bcftools view -r chr12 dbsnp.bgz -Oz -o known_sites_chr12.vcf.gz
tabix -p vcf known_sites_chr12.vcf.gz

# -------------------------------
# BASE RECALIBRATION
# -------------------------------
echo "[11/20] Base recalibration..."
docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk BaseRecalibrator \
  -I /data/${SAMPLE}.sorted.dedup.bam \
  -R /data/$REF \
  --known-sites /data/known_sites_chr12.vcf.gz \
  -O /data/${SAMPLE}.recal.table

# -------------------------------
# APPLY BQSR
# -------------------------------
echo "[12/20] Applying BQSR..."
docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk ApplyBQSR \
  -I /data/${SAMPLE}.sorted.dedup.bam \
  -R /data/$REF \
  --bqsr-recal-file /data/${SAMPLE}.recal.table \
  -O /data/${SAMPLE}.BQSR.bam

samtools index ${SAMPLE}.BQSR.bam

# -------------------------------
# VARIANT CALLING
# -------------------------------
echo "[13/20] Calling variants..."
docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk HaplotypeCaller \
  -R /data/$REF \
  -I /data/${SAMPLE}.BQSR.bam \
  -O /data/${SAMPLE}.raw.vcf

# -------------------------------
# SELECT SNPS
# -------------------------------
echo "[14/20] Selecting SNPs..."
docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk SelectVariants \
  -R /data/$REF \
  -V /data/${SAMPLE}.raw.vcf \
  --select-type SNP \
  -O /data/${SAMPLE}.snps.vcf

# -------------------------------
# FILTER SNPS
# -------------------------------
echo "[15/20] Filtering SNPs..."
docker run --rm -v "$PWD":/data $GATK_IMG \
  gatk VariantFiltration \
  -R /data/$REF \
  -V /data/${SAMPLE}.snps.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "LOW_QUALITY" \
  -O /data/${SAMPLE}.analysis_ready_snps.vcf

# -------------------------------
# ANNOVAR PREPARATION
# -------------------------------
echo "[16/20] Preparing ANNOVAR input..."
convert2annovar.pl -format vcf4 ${SAMPLE}.analysis_ready_snps.vcf \
  > ${SAMPLE}.analysis_ready_snps.avinput

# -------------------------------
# KRAS EXTRACTION
# -------------------------------
echo "[17/20] Extracting KRAS variants (chr12:25205246-25250929)..."
bcftools view -r chr12:25205246-25250929 \
  ${SAMPLE}.analysis_ready_snps.vcf \
  > ${SAMPLE}.KRAS_variants.vcf

# -------------------------------
# FINAL REPORT
# -------------------------------
echo "[18/20] Pipeline complete!"
echo "============================================"
echo "Key Outputs:"
echo " - Final SNPs: ${SAMPLE}.analysis_ready_snps.vcf"
echo " - ANNOVAR Input: ${SAMPLE}.analysis_ready_snps.avinput"
echo " - KRAS Variants: ${SAMPLE}.KRAS_variants.vcf"
echo "============================================"


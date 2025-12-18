#!/bin/bash
set -euo pipefail

############################################
# VARIABLES
############################################

SAMPLE="SRR35521082"
THREADS=4

BASE_DIR="$PWD"
RAW_DIR="${BASE_DIR}/Raw_Data"
OUT_DIR="${BASE_DIR}/Outputs"
REF_DIR="${OUT_DIR}"

FASTQ1="${RAW_DIR}/${SAMPLE}_1_300k.fastq"
FASTQ2="${RAW_DIR}/${SAMPLE}_2_300k.fastq"

REF_FA="${REF_DIR}/chr12.fa"
DBSNP_VCF="${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
KNOWN_CHR12_VCF="${REF_DIR}/known_sites_chr12.vcf.gz"
KNOWN_CHR12_NC_VCF="${REF_DIR}/known_sites_chr12_NC.vcf.gz"

############################################
# DIRECTORY SETUP
############################################

mkdir -p "${RAW_DIR}" "${OUT_DIR}"

############################################
# STEP 1: ENVIRONMENT CHECK
############################################

docker --version
fastqc --version
bwa 2>&1 | head -n 1
samtools --version | head -n 1

############################################
# STEP 2: FASTQ DOWNLOAD
############################################

wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/${SAMPLE}/${SAMPLE}_1.fastq.gz -P "${RAW_DIR}"
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR355/082/${SAMPLE}/${SAMPLE}_2.fastq.gz -P "${RAW_DIR}"

############################################
# STEP 3: SUBSAMPLING (300K READS)
############################################

seqtk sample -s123 "${RAW_DIR}/${SAMPLE}_1.fastq.gz" 300000 > "${FASTQ1}"
seqtk sample -s123 "${RAW_DIR}/${SAMPLE}_2.fastq.gz" 300000 > "${FASTQ2}"

############################################
# STEP 4: QUALITY CONTROL
############################################

fastqc -o "${OUT_DIR}" "${FASTQ1}" "${FASTQ2}"

############################################
# STEP 5: CHR12 REFERENCE PREPARATION
############################################

wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/seqs_for_alignment_pipelines.ucsc_ids/chr12.fna.gz
gunzip -f chr12.fna.gz
mv chr12.fna "${REF_FA}"

bwa index "${REF_FA}"
samtools faidx "${REF_FA}"

docker run --rm -v "${REF_DIR}":/data broadinstitute/gatk:latest \
  gatk CreateSequenceDictionary \
  -R /data/chr12.fa \
  -O /data/chr12.dict

############################################
# STEP 6: ALIGNMENT
############################################

bwa mem -t "${THREADS}" \
  -R "@RG\tID:${SAMPLE}\tPL:ILLUMINA\tSM:${SAMPLE}" \
  "${REF_FA}" "${FASTQ1}" "${FASTQ2}" | \
  samtools sort -o "${OUT_DIR}/${SAMPLE}.sorted.bam"

samtools index "${OUT_DIR}/${SAMPLE}.sorted.bam"

############################################
# STEP 7: MARK DUPLICATES
############################################

docker run --rm -v "${OUT_DIR}":/data broadinstitute/gatk:latest \
  gatk MarkDuplicates \
  -I /data/${SAMPLE}.sorted.bam \
  -O /data/${SAMPLE}.sorted.dedup.bam \
  -M /data/${SAMPLE}.metrics.txt

samtools index "${OUT_DIR}/${SAMPLE}.sorted.dedup.bam"

############################################
# STEP 8: dbSNP → chr12 KNOWN SITES
############################################

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz -P "${REF_DIR}"
bcftools index "${DBSNP_VCF}"

bcftools view -r chr12 "${DBSNP_VCF}" -Oz -o "${KNOWN_CHR12_VCF}"
bcftools index "${KNOWN_CHR12_VCF}"

############################################
# STEP 9: CONTIG RENAMING (chr12 → NC_)
############################################

cat > "${REF_DIR}/contig_rename.txt" <<EOF
chr12 NC_000012.12
EOF

bcftools annotate \
  --rename-chrs "${REF_DIR}/contig_rename.txt" \
  -Oz -o "${KNOWN_CHR12_NC_VCF}" "${KNOWN_CHR12_VCF}"

tabix -p vcf "${KNOWN_CHR12_NC_VCF}"

############################################
# STEP 10: BQSR
############################################

docker run --rm -v "${OUT_DIR}":/data -v "${REF_DIR}":/ref broadinstitute/gatk:latest \
  gatk BaseRecalibrator \
  -R /ref/chr12.fa \
  -I /data/${SAMPLE}.sorted.dedup.bam \
  --known-sites /ref/known_sites_chr12_NC.vcf.gz \
  -O /data/${SAMPLE}.recal_data.table

docker run --rm -v "${OUT_DIR}":/data -v "${REF_DIR}":/ref broadinstitute/gatk:latest \
  gatk ApplyBQSR \
  -R /ref/chr12.fa \
  -I /data/${SAMPLE}.sorted.dedup.bam \
  --bqsr-recal-file /data/${SAMPLE}.recal_data.table \
  -O /data/${SAMPLE}.sorted.dedup.BQSR.bam

############################################
# STEP 11: VARIANT CALLING
############################################

docker run --rm -v "${OUT_DIR}":/data -v "${REF_DIR}":/ref broadinstitute/gatk:latest \
  gatk HaplotypeCaller \
  -R /ref/chr12.fa \
  -I /data/${SAMPLE}.sorted.dedup.BQSR.bam \
  -O /data/${SAMPLE}.raw.vcf

############################################
# STEP 12: SNP SELECTION & FILTERING
############################################

docker run --rm -v "${OUT_DIR}":/data -v "${REF_DIR}":/ref broadinstitute/gatk:latest \
  gatk SelectVariants \
  -R /ref/chr12.fa \
  -V /data/${SAMPLE}.raw.vcf \
  --select-type SNP \
  -O /data/${SAMPLE}.raw_snps.vcf

docker run --rm -v "${OUT_DIR}":/data -v "${REF_DIR}":/ref broadinstitute/gatk:latest \
  gatk VariantFiltration \
  -R /ref/chr12.fa \
  -V /data/${SAMPLE}.raw_snps.vcf \
  -O /data/${SAMPLE}.filtered_snps.vcf \
  -filter "QD < 2.0" --filter-name QD \
  -filter "FS > 60.0" --filter-name FS \
  -filter "MQ < 40.0" --filter-name MQ

docker run --rm -v "${OUT_DIR}":/data broadinstitute/gatk:latest \
  gatk SelectVariants \
  --exclude-filtered \
  -V /data/${SAMPLE}.filtered_snps.vcf \
  -O /data/${SAMPLE}.analysis_ready_snps.vcf

############################################
# STEP 13: ANNOVAR INPUT
############################################

docker run --rm -v "${OUT_DIR}":/data bioinfochrustrasbourg/annovar:latest \
  perl convert2annovar.pl -format vcf4 \
  /data/${SAMPLE}.analysis_ready_snps.vcf \
  -outfile /data/${SAMPLE}.analysis_ready_snps.avinput

############################################
# COMPLETION
############################################

echo "CRC–KRAS chr12 pipeline completed successfully"

#!/bin/bash
set -euo pipefail

############################################
# KRAS G12V Variant Discovery Pipeline
# Author: Nishitha N
############################################

# -------- CONFIGURATION --------
THREADS=8
REF="chr12.fa"
R1="sample_R1.fastq.gz"
R2="sample_R2.fastq.gz"

WORKDIR="$PWD"
OUTDIR="$WORKDIR/output"
GATK_IMAGE="broadinstitute/gatk:latest"

POSITIONS_FILE="positions.txt"
KNOWN_SITES="known_sites_chr12.vcf.gz"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== KRAS G12V PIPELINE STARTED ==="

# -------- STEP 1: QC --------
echo "[1/14] Running FastQC..."
fastqc "$R1" "$R2"

# -------- STEP 2: TRIMMING --------
echo "[2/14] Trimming reads..."
fastp \
  -i "$R1" \
  -I "$R2" \
  -o R1_trimmed.fastq.gz \
  -O R2_trimmed.fastq.gz

# -------- STEP 3: REFERENCE INDEXING --------
echo "[3/14] Indexing reference..."
bwa index "$REF"
samtools faidx "$REF"

# -------- STEP 4: ALIGNMENT --------
echo "[4/14] Aligning reads..."
bwa mem -t "$THREADS" "$REF" \
  R1_trimmed.fastq.gz \
  R2_trimmed.fastq.gz \
  > aligned.sam

# -------- STEP 5: BAM PROCESSING --------
echo "[5/14] Converting and sorting BAM..."
samtools view -Sb aligned.sam | samtools sort -o sorted.bam
samtools index sorted.bam

# -------- STEP 6: MARK DUPLICATES --------
echo "[6/14] Marking duplicates with GATK..."
docker run --rm -v "$PWD":/data -w /data \
  "$GATK_IMAGE" \
  gatk MarkDuplicatesSpark \
  -I sorted.bam \
  -O SRR35521082.sorted.dedup.bam \
  -M metrics.txt

samtools index SRR35521082.sorted.dedup.bam

# -------- STEP 7: VARIANT CALLING --------
echo "[7/14] Calling variants..."
docker run --rm -v "$PWD":/data -w /data \
  "$GATK_IMAGE" \
  gatk HaplotypeCaller \
  -R "$REF" \
  -I SRR35521082.sorted.dedup.bam \
  -O raw_snps.vcf

# -------- STEP 8: REGION-SPECIFIC EXTRACTION --------
echo "[8/14] Creating BED file from positions.txt..."
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' "$POSITIONS_FILE" > positions.bed

echo "[8/14] Extracting region-specific variants..."
for v in raw_snps.vcf; do
  bcftools view -R positions.bed "$v" -Oz > "${v%.vcf}.gz"
done

# -------- STEP 9: VCF COMPARISON --------
echo "[9/14] Generating presence/absence matrix..."
echo -e "chrom\tpos\traw_snps" > positions_presence.tsv

while read chr pos; do
  echo -e "$chr\t$pos\t"$(
    bcftools view -r "$chr:$pos-$pos" raw_snps.vcf.gz | grep -v "^#" | wc -l
  )
done < "$POSITIONS_FILE" >> positions_presence.tsv

# -------- STEP 10: CONTIG MAPPING --------
echo "[10/14] Detecting contig naming..."
zcat raw_snps.vcf.gz | grep '^##contig' > contigs.header

contig_id=$(grep '^##contig' contigs.header | \
  sed 's/.*ID=\([^,>]*\).*/\1/p' | head -n 1)

echo "[10/14] Mapping positions to contig: $contig_id"

awk -v cid="$contig_id" 'BEGIN{OFS="\t"}{
  if ($1 == "chr12") c = cid;
  else c = $1;
  print c, $2-1, $2
}' "$POSITIONS_FILE" > positions_mapped.bed

# -------- STEP 11: KRAS G12V EXTRACTION --------
echo "[11/14] Extracting KRAS G12V variant..."
bcftools view -r chr12:25245350 "$KNOWN_SITES" > kras_g12v.vcf

# -------- STEP 12: DOWNLOAD STRUCTURE --------
echo "[12/14] Downloading KRAS PDB structure..."
wget -q https://files.rcsb.org/download/6OIM.pdb

# -------- STEP 13: ALLELE COUNTS --------
echo "[13/14] Extracting allele counts and frequencies..."
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AF\n" \
  -R positions_mapped.bed \
  "$KNOWN_SITES" >> position_ac_mapping.tsv

# -------- STEP 14: STRUCTURAL ANALYSIS --------
echo "[14/14] Running structural analysis script..."
python3 kras_structural_analysis.py

echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "Results saved in: $OUTDIR"
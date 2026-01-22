#!/bin/bash

# Biomarker Discovery Pipeline Script
# This script automates the colorectal cancer biomarker discovery workflow
# Usage: ./biomarker_discovery.sh

set -e  # Exit on error

echo "Starting Biomarker Discovery Pipeline..."

# Step 1: Inspect VCF
echo "Step 1: Inspecting VCF header..."
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | head -20

echo "Step 2: Viewing variant records..."
zcat Outputs/SRR35521082.analysis_ready_snps.vcf.gz | grep -v "^#" | head -10

# Step 3: Recompress with BGZIP
echo "Step 3: Recompression with BGZIP..."
gunzip -c Outputs/SRR35521082.analysis_ready_snps.vcf.gz | bgzip > Outputs/SRR35521082.analysis_ready_snps.vcf.gz.tmp
mv Outputs/SRR35521082.analysis_ready_snps.vcf.gz.tmp Outputs/SRR35521082.analysis_ready_snps.vcf.gz

# Step 4: Index VCF
echo "Step 4: Indexing VCF..."
bcftools index Outputs/SRR35521082.analysis_ready_snps.vcf.gz

# Step 5: Annotate with VEP
echo "Step 5: Annotating with VEP..."
docker run --rm -v "$PWD":/data ensemblorg/ensembl-vep \
  vep -i /data/Outputs/SRR35521082.analysis_ready_snps.vcf.gz \
  -o /data/vep_results.tsv \
  --assembly GRCh38 \
  --symbol --protein --canonical --tab --database --force_overwrite

# Step 6: Generate cancer biomarker table
echo "Step 6: Generating cancer biomarker table..."
python3 cancer_variants_table.py

# Step 7: Final annotation
echo "Step 7: Final annotation..."
python3 annotate_variants.py

echo "Pipeline completed successfully!"
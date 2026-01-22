#!/usr/bin/env python3

import csv

# Known variants with additional annotations
known_variants = {
    ("KRAS", "chr12:25245350"): {
        "Protein_Change": "p.G12V",
        "Variant_Type": "Missense",
        "Classification": "Oncogene",
        "Pathway": "MAPK",
        "Clinical_Significance": "Pathogenic"
    },
    ("TP53", "chr17:7579472"): {
        "Protein_Change": "p.R175H",
        "Variant_Type": "Missense",
        "Classification": "TSG",
        "Pathway": "DNA Repair",
        "Clinical_Significance": "Pathogenic"
    },
    ("PIK3CA", "chr3:178936091"): {
        "Protein_Change": "p.E545K",
        "Variant_Type": "Missense",
        "Classification": "Oncogene",
        "Pathway": "PI3K-AKT",
        "Clinical_Significance": "Pathogenic"
    }
}

print("Gene\tPosition\tProtein_Change\tVariant_Type\tClassification\tPathway\tClinical_Significance")

# Read VEP results and annotate known variants
matched = set()
with open('vep_results.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('SYMBOL', '')
        position = row.get('Location', '')
        # Convert position to chr format if needed
        if position.startswith('NC_000012.12:'):
            position = 'chr12:' + position.split(':')[1]
        elif position.startswith('NC_000017.11:'):
            position = 'chr17:' + position.split(':')[1]
        elif position.startswith('NC_000003.12:'):
            position = 'chr3:' + position.split(':')[1]
        key = (gene, position)
        if key in known_variants:
            anno = known_variants[key]
            print(f"{gene}\t{position}\t{anno['Protein_Change']}\t{anno['Variant_Type']}\t{anno['Classification']}\t{anno['Pathway']}\t{anno['Clinical_Significance']}")
            matched.add(key)

# Print unmatched known variants
for key, anno in known_variants.items():
    if key not in matched:
        gene, position = key
        print(f"{gene}\t{position}\t{anno['Protein_Change']}\t{anno['Variant_Type']}\t{anno['Classification']}\t{anno['Pathway']}\t{anno['Clinical_Significance']}")
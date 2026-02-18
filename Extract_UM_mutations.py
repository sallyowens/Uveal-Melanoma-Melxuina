#!/usr/bin/env python3
"""
Extract Uveal Melanoma Driver Mutations from Novogene Somatic Annotation Files

This script processes Novogene WGS somatic mutation annotation files (.xls.gz) 
and extracts mutations in key uveal melanoma driver genes:
- GNAQ, GNA11 (early drivers)
- BAP1, SF3B1, EIF1AX (metastatic risk markers)
- PRAME, MBD4 (additional prognostic markers)

Outputs:
1. Detailed TSV with all mutations (gene, position, amino acid change, VAF, etc.)
2. Clean matrix TSV ready for Excel/SPSS with separate mutation and VAF columns per gene

Author: Sally Owens
Date: February 2026
"""

import os
import gzip
import pandas as pd
from pathlib import Path
import argparse


# Uveal melanoma driver genes
UM_GENES = ['GNAQ', 'GNA11', 'BAP1', 'PRAME', 'MBD4', 'SF3B1', 'EIF1AX']


def extract_mutations_from_file(filepath, var_type, um_genes):
    """
    Extract mutations in UM genes from a single Novogene annotation file.
    
    Args:
        filepath: Path to .xls.gz annotation file
        var_type: 'SNV' or 'INDEL'
        um_genes: List of gene names to extract
    
    Returns:
        List of mutation dictionaries
    """
    sample = filepath.name.split('.')[0]
    mutations = []
    
    try:
        # Read compressed annotation file
        with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
            df = pd.read_csv(f, sep='\t', low_memory=False)
        
        # Filter for PASS variants only
        df = df[df['FILTER'] == 'PASS']
        
        # Filter for UM genes using GeneName column (not Gene which contains transcript IDs)
        df_um = df[df['GeneName'].isin(um_genes)].copy()
        
        # Keep only exonic and splicing variants
        df_um = df_um[df_um['Func'].isin(['exonic', 'splicing', 'exonic;splicing'])]
        
        if len(df_um) == 0:
            return mutations
        
        # Find tumour column (starts with T_)
        tumor_col = next((c for c in df.columns if c.startswith('T_')), None)
        
        for _, row in df_um.iterrows():
            gene   = row.get('GeneName', 'NA')
            exfunc = row.get('ExonicFunc', 'NA')
            aac    = str(row.get('AAChange', ''))
            chrom  = row.get('CHROM', 'NA')
            pos    = row.get('POS', 'NA')
            ref    = row.get('REF', 'NA')
            alt    = row.get('ALT', 'NA')
            cosmic = row.get('cosmic', '.')
            
            # Extract VAF (Variant Allele Frequency) from tumour FORMAT field
            vaf = 'NA'
            if tumor_col and pd.notna(row.get(tumor_col)):
                fmt_keys = str(row.get('FORMAT', '')).split(':')
                fmt_vals = str(row[tumor_col]).split(':')
                fmt_dict = dict(zip(fmt_keys, fmt_vals))
                vaf = fmt_dict.get('AF', 'NA')
                try:
                    vaf = f"{float(vaf):.1%}"
                except:
                    pass
            
            # Extract amino acid change from AAChange field
            aa = 'NA'
            if aac not in ['.', 'NA', '', 'nan']:
                # AAChange format: GENE:TRANSCRIPT:EXON:cDNA:PROTEIN
                # e.g., GNA11:NM_002067:exon5:c.A626T:p.Q209L
                parts = aac.split(',')
                for p in reversed(parts):
                    if 'p.' in p:
                        aa = p.split(':')[-1]
                        break
                if aa == 'NA':
                    aa = aac[:40]  # Fallback to truncated full string
            
            mutations.append({
                'Sample'   : sample,
                'Type'     : var_type,
                'Gene'     : gene,
                'Function' : exfunc,
                'AA_Change': aa,
                'VAF'      : vaf,
                'Chr'      : chrom,
                'Pos'      : pos,
                'Ref'      : ref,
                'Alt'      : alt,
                'COSMIC'   : cosmic if pd.notna(cosmic) else '.'
            })
    
    except Exception as e:
        print(f"ERROR processing {sample}: {e}")
    
    return mutations


def create_mutation_matrix(df_mutations, um_genes):
    """
    Create a clean wide-format matrix with one row per sample.
    For each gene: separate columns for mutation name and VAF.
    
    Args:
        df_mutations: DataFrame with all mutations
        um_genes: List of gene names
    
    Returns:
        DataFrame in wide format
    """
    all_samples = sorted(df_mutations['Sample'].unique())
    rows = []
    
    for sample in all_samples:
        row = {'Sample': sample}
        sample_data = df_mutations[df_mutations['Sample'] == sample]
        
        for gene in um_genes:
            gene_data = sample_data[sample_data['Gene'] == gene]
            
            if len(gene_data) == 0:
                row[f'{gene}_Mutation'] = 'WT'
                row[f'{gene}_VAF'] = ''
            else:
                # If multiple mutations in same gene, take the one with highest VAF
                gene_data = gene_data.copy()
                gene_data['vaf_num'] = gene_data['VAF'].str.rstrip('%').apply(
                    lambda x: float(x) if x not in ['NA', ''] else 0
                )
                best = gene_data.sort_values('vaf_num', ascending=False).iloc[0]
                row[f'{gene}_Mutation'] = best['AA_Change']
                row[f'{gene}_VAF'] = best['VAF']
        
        rows.append(row)
    
    matrix = pd.DataFrame(rows)
    
    # Reorder columns: Sample, then for each gene: Mutation then VAF
    cols = ['Sample']
    for gene in um_genes:
        cols.append(f'{gene}_Mutation')
        cols.append(f'{gene}_VAF')
    matrix = matrix[cols]
    
    return matrix


def main():
    parser = argparse.ArgumentParser(
        description='Extract UM driver mutations from Novogene somatic annotation files'
    )
    parser.add_argument(
        '--snv-dir',
        required=True,
        help='Directory containing SNV annotation files (*.somatic.snv.annovar*.xls.gz)'
    )
    parser.add_argument(
        '--indel-dir',
        required=True,
        help='Directory containing INDEL annotation files (*.somatic.indel.annovar*.xls.gz)'
    )
    parser.add_argument(
        '--output-prefix',
        default='UM_mutations',
        help='Prefix for output files (default: UM_mutations)'
    )
    
    args = parser.parse_args()
    
    # Collect all annotation files
    snv_dir = Path(args.snv_dir)
    indel_dir = Path(args.indel_dir)
    
    ann_files = []
    for f in snv_dir.glob('*.xls.gz'):
        if not f.name.startswith('._'):
            ann_files.append(('SNV', f))
    for f in indel_dir.glob('*.xls.gz'):
        if not f.name.startswith('._'):
            ann_files.append(('INDEL', f))
    
    print(f"Found {len(ann_files)} annotation files")
    print(f"Processing {len(set(f.name.split('.')[0] for _, f in ann_files))} unique samples")
    print("")
    
    # Extract mutations from all files
    all_mutations = []
    for var_type, filepath in sorted(ann_files, key=lambda x: x[1].name):
        sample = filepath.name.split('.')[0]
        print(f"Processing: {sample} ({var_type})")
        mutations = extract_mutations_from_file(filepath, var_type, UM_GENES)
        
        if mutations:
            for m in mutations:
                print(f"  ✓ {m['Gene']}  {m['AA_Change']}  VAF={m['VAF']}")
        
        all_mutations.extend(mutations)
    
    print("")
    
    if not all_mutations:
        print("No mutations found in UM driver genes")
        return
    
    # Create DataFrames
    df_mutations = pd.DataFrame(all_mutations)
    
    # Save detailed output
    detail_file = f"{args.output_prefix}_detailed.tsv"
    df_mutations.to_csv(detail_file, sep='\t', index=False)
    
    # Create and save mutation matrix
    matrix = create_mutation_matrix(df_mutations, UM_GENES)
    matrix_file = f"{args.output_prefix}_matrix.tsv"
    matrix.to_csv(matrix_file, sep='\t', index=False)
    
    # Print summary
    all_samples = sorted(df_mutations['Sample'].unique())
    
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total samples analysed: {len(all_samples)}")
    print(f"Total mutations found: {len(all_mutations)}")
    print("")
    print("Mutations per gene:")
    for gene in UM_GENES:
        n = len(df_mutations[df_mutations['Gene'] == gene]['Sample'].unique())
        pct = (n / len(all_samples)) * 100
        print(f"  {gene:8s}: {n:2d}/{len(all_samples)} ({pct:.0f}%)")
    print("")
    print("Output files:")
    print(f"  {detail_file}  <- Full details (all mutations)")
    print(f"  {matrix_file}  <- Clean matrix (open in Excel)")
    print("=" * 70)


if __name__ == '__main__':
    main()

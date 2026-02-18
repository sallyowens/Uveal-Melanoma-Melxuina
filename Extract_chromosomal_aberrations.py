#!/usr/bin/env python3
"""
Extract Chromosomal Aberrations from CNVkit and FACETS Results

This script processes CNVkit (.cns files) and FACETS (segments_*.csv) outputs
to create a summary table of key chromosomal changes in uveal melanoma:

Key UM chromosomal aberrations:
- Monosomy 3 (chromosome 3 loss) - strong metastatic predictor
- 8q gain - associated with metastasis
- 6p gain - additional risk marker
- 1p loss - poor prognosis marker
- Tumor purity and ploidy (from FACETS)

Author: Sally Owens
Date: February 2026
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')


# Chromosome lengths (hg38)
CHR_LENGTHS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}


def calculate_chr_metrics(cns_file, chromosome):
    """
    Calculate copy number metrics for a specific chromosome from CNVkit .cns file.
    
    Returns dict with: mean_cn, pct_loss, pct_gain, status
    """
    try:
        df = pd.read_csv(cns_file, sep='\t', comment='#')
        
        # Handle different column names
        chr_col = 'chromosome' if 'chromosome' in df.columns else 'chrom'
        
        # Filter for this chromosome
        chr_data = df[df[chr_col].astype(str).str.replace('chr', '') == 
                      chromosome.replace('chr', '')].copy()
        
        if len(chr_data) == 0:
            return None
        
        # Calculate segment lengths
        chr_data['length'] = chr_data['end'] - chr_data['start']
        chr_length = CHR_LENGTHS.get(chromosome, chr_data['end'].max())
        
        # Determine which CN column to use
        if 'cn' in chr_data.columns:
            cn_col = 'cn'
            loss_threshold = 2
            gain_threshold = 2
        elif 'log2' in chr_data.columns:
            # Convert log2 to copy number
            chr_data['cn_status'] = 2 * (2 ** chr_data['log2'])
            cn_col = 'cn_status'
            loss_threshold = 1.5
            gain_threshold = 2.5
        else:
            return None
        
        # Calculate loss and gain percentages
        loss_data = chr_data[chr_data[cn_col] < loss_threshold]
        gain_data = chr_data[chr_data[cn_col] > gain_threshold]
        
        loss_length = loss_data['length'].sum()
        gain_length = gain_data['length'].sum()
        
        pct_loss = (loss_length / chr_length) * 100
        pct_gain = (gain_length / chr_length) * 100
        
        # Weighted mean CN
        weighted_cn = (chr_data[cn_col] * chr_data['length']).sum() / chr_data['length'].sum()
        
        return {
            'mean_cn': round(weighted_cn, 2),
            'pct_loss': round(pct_loss, 2),
            'pct_gain': round(pct_gain, 2)
        }
    
    except Exception as e:
        print(f"Error processing {cns_file} for {chromosome}: {e}")
        return None


def classify_chr3_status(metrics):
    """Classify chromosome 3 status based on CN metrics."""
    if metrics is None:
        return 'Unknown'
    
    mean_cn = metrics['mean_cn']
    pct_loss = metrics['pct_loss']
    
    if mean_cn < 1.3 and pct_loss > 80:
        return 'Monosomy'
    elif mean_cn > 1.7 and pct_loss < 20:
        return 'Disomy'
    elif pct_loss > 20:
        return 'Partial_Loss'
    else:
        return 'Disomy'


def classify_8q_status(chr8_metrics):
    """Classify 8q status (8q is the long arm, positions ~48M-145M)."""
    if chr8_metrics is None:
        return 'Normal'
    
    # For 8q gain, we care about the mean CN and gain percentage
    mean_cn = chr8_metrics['mean_cn']
    pct_gain = chr8_metrics['pct_gain']
    
    if mean_cn > 2.5 or pct_gain > 50:
        return 'Gain'
    else:
        return 'Normal'


def classify_6p_status(chr6_metrics):
    """Classify 6p status (6p is short arm, positions 0-~60M)."""
    if chr6_metrics is None:
        return 'Normal'
    
    mean_cn = chr6_metrics['mean_cn']
    pct_gain = chr6_metrics['pct_gain']
    
    if mean_cn > 2.5 or pct_gain > 30:
        return 'Gain'
    else:
        return 'Normal'


def classify_1p_status(chr1_metrics):
    """Classify 1p status (1p is short arm, positions 0-~125M)."""
    if chr1_metrics is None:
        return 'Normal'
    
    pct_loss = chr1_metrics['pct_loss']
    
    if pct_loss > 30:
        return 'Loss'
    else:
        return 'Normal'


def get_facets_data(facets_file):
    """Extract purity and ploidy from FACETS segments file."""
    try:
        df = pd.read_csv(facets_file)
        
        # FACETS summary info is typically in the first few rows or separate file
        # Check if we have purity/ploidy columns
        purity = None
        ploidy = None
        
        # Sometimes purity is in a separate summary file
        # For now, we'll try to extract from the segments file
        # This may need adjustment based on your actual FACETS output format
        
        if 'cf' in df.columns or 'cf.em' in df.columns:
            # cf = cellular fraction (purity estimate for each segment)
            cf_col = 'cf.em' if 'cf.em' in df.columns else 'cf'
            purity = df[cf_col].median() * 100  # Convert to percentage
        
        return {
            'purity': round(purity, 1) if purity else None,
            'ploidy': round(ploidy, 1) if ploidy else None
        }
    
    except Exception as e:
        print(f"Warning: Could not extract FACETS data from {facets_file}: {e}")
        return {'purity': None, 'ploidy': None}


def process_cnvkit_file(cns_file):
    """Process a single CNVkit .cns file and extract all chromosomal metrics."""
    sample = Path(cns_file).stem.replace('.cns', '').replace('.call', '').replace('.revised', '')
    
    result = {'Sample': sample}
    
    # Chromosome 3 (monosomy 3)
    chr3_metrics = calculate_chr_metrics(cns_file, 'chr3')
    result['Chr3_Status'] = classify_chr3_status(chr3_metrics)
    result['Chr3_Mean_CN'] = chr3_metrics['mean_cn'] if chr3_metrics else 'NA'
    result['Chr3_Pct_Loss'] = chr3_metrics['pct_loss'] if chr3_metrics else 'NA'
    
    # Chromosome 8 (8q gain)
    chr8_metrics = calculate_chr_metrics(cns_file, 'chr8')
    result['Chr8_Status'] = classify_8q_status(chr8_metrics)
    result['Chr8_Mean_CN'] = chr8_metrics['mean_cn'] if chr8_metrics else 'NA'
    result['Chr8_Pct_Gain'] = chr8_metrics['pct_gain'] if chr8_metrics else 'NA'
    
    # Chromosome 6 (6p gain)
    chr6_metrics = calculate_chr_metrics(cns_file, 'chr6')
    result['Chr6_Status'] = classify_6p_status(chr6_metrics)
    result['Chr6_Mean_CN'] = chr6_metrics['mean_cn'] if chr6_metrics else 'NA'
    
    # Chromosome 1 (1p loss)
    chr1_metrics = calculate_chr_metrics(cns_file, 'chr1')
    result['Chr1_Status'] = classify_1p_status(chr1_metrics)
    result['Chr1_Mean_CN'] = chr1_metrics['mean_cn'] if chr1_metrics else 'NA'
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Extract chromosomal aberrations from CNVkit and FACETS results'
    )
    parser.add_argument(
        '--cnvkit-dir',
        required=True,
        help='Directory containing CNVkit .cns files'
    )
    parser.add_argument(
        '--facets-dir',
        help='Optional directory containing FACETS segments_*.csv files'
    )
    parser.add_argument(
        '--output',
        default='chromosome_aberrations.tsv',
        help='Output file name (default: chromosome_aberrations.tsv)'
    )
    
    args = parser.parse_args()
    
    # Find all CNVkit files
    cnvkit_dir = Path(args.cnvkit_dir)
    cns_files = list(cnvkit_dir.glob('*.cns'))
    
    if not cns_files:
        print(f"ERROR: No .cns files found in {cnvkit_dir}")
        return
    
    print(f"Found {len(cns_files)} CNVkit .cns files")
    print("")
    
    # Process each CNVkit file
    results = []
    for cns_file in sorted(cns_files):
        sample = Path(cns_file).stem.replace('.cns', '').replace('.call', '').replace('.revised', '')
        print(f"Processing: {sample}")
        
        result = process_cnvkit_file(cns_file)
        
        # Add FACETS data if available
        if args.facets_dir:
            facets_dir = Path(args.facets_dir)
            facets_file = facets_dir / f"segments_{sample}.csv"
            
            if facets_file.exists():
                facets_data = get_facets_data(facets_file)
                result['Purity'] = facets_data['purity']
                result['Ploidy'] = facets_data['ploidy']
            else:
                result['Purity'] = None
                result['Ploidy'] = None
        
        results.append(result)
        
        # Print key findings
        if result['Chr3_Status'] == 'Monosomy':
            print(f"  ⚠ Monosomy 3 detected")
        if result['Chr8_Status'] == 'Gain':
            print(f"  ⚠ 8q gain detected")
    
    print("")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Save output
    df.to_csv(args.output, sep='\t', index=False)
    
    # Print summary
    print("=" * 70)
    print("CHROMOSOMAL ABERRATION SUMMARY")
    print("=" * 70)
    print(f"Total samples: {len(results)}")
    print("")
    
    print("Chromosome 3 status:")
    print(df['Chr3_Status'].value_counts().to_string())
    print("")
    
    print("Chromosome 8 status:")
    print(df['Chr8_Status'].value_counts().to_string())
    print("")
    
    if 'Chr6_Status' in df.columns:
        chr6_gain = len(df[df['Chr6_Status'] == 'Gain'])
        if chr6_gain > 0:
            print(f"6p gain: {chr6_gain} samples")
    
    if 'Chr1_Status' in df.columns:
        chr1_loss = len(df[df['Chr1_Status'] == 'Loss'])
        if chr1_loss > 0:
            print(f"1p loss: {chr1_loss} samples")
    
    print("")
    print(f"Output saved to: {args.output}")
    print("=" * 70)


if __name__ == '__main__':
    main()

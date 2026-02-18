# Uveal Melanoma Driver Mutation Extractor

Extract somatic mutations in key uveal melanoma driver genes from Novogene WGS annotation files.

## Overview

This tool processes Novogene's somatic mutation annotation files (`.muTect2.somatic.*.annovar.hg38_multianno.xls.gz`) and extracts mutations in 7 key uveal melanoma driver genes:

- **GNAQ, GNA11** - Early drivers (~90% of UMs)
- **BAP1** - Loss associated with metastatic disease
- **SF3B1, EIF1AX** - Alternative splicing mutations, lower metastatic risk
- **PRAME, MBD4** - Additional prognostic markers

## Features

- Processes both SNV and INDEL annotation files
- Filters for PASS-quality variants only
- Extracts exonic and splicing variants
- Reports amino acid changes and VAF (variant allele frequency)
- Generates two output formats:
  1. Detailed TSV with all mutation information
  2. Clean matrix TSV ready for Excel/statistical analysis

## Requirements

```bash
pip install pandas
```

Python 3.7+

## Usage

```bash
python extract_um_mutations.py \
  --snv-dir /path/to/SNV/annotations \
  --indel-dir /path/to/INDEL/annotations \
  --output-prefix UM_mutations
```

### Input Files

The script expects Novogene annotation files with this naming pattern:
```
T_25RV04.muTect2.somatic.snv.annovar.hg38_multianno.xls.gz
T_25RV04.muTect2.somatic.indel.annovar.hg38_multianno.xls.gz
```

### Output Files

**1. Detailed output (`UM_mutations_detailed.tsv`)**

All mutations with full annotation:

| Sample | Type | Gene | Function | AA_Change | VAF | Chr | Pos | Ref | Alt | COSMIC |
|---|---|---|---|---|---|---|---|---|---|---|
| T_22RV17 | SNV | GNA11 | missense SNV | p.Q209L | 37.9% | chr19 | 3118942 | A | T | COSM... |
| T_25RV08 | SNV | BAP1 | missense SNV | p.L49P | 90.5% | chr3 | 52436823 | T | C | ... |

**2. Matrix output (`UM_mutations_matrix.tsv`)**

Clean wide format ready for Excel/SPSS:

| Sample | GNAQ_Mutation | GNAQ_VAF | GNA11_Mutation | GNA11_VAF | BAP1_Mutation | BAP1_VAF | ... |
|---|---|---|---|---|---|---|---|
| T_22RV17 | WT | | p.Q209L | 37.9% | p.S596Tfs*21 | 77.3% | ... |
| T_25RV08 | WT | | WT | | p.L49P | 90.5% | ... |

## Understanding VAF (Variant Allele Frequency)

**VAF = percentage of DNA reads carrying the mutation**

### Interpretation Guide

- **40-60% VAF**: Clonal heterozygous mutation in reasonably pure tumour
- **70-100% VAF**: Very pure tumour OR loss of other allele (e.g., BAP1 with monosomy 3)
- **25-40% VAF**: Lower purity tumour or emerging subclone
- **<20% VAF**: Subclonal mutation, contamination, or potential artifact

### Example

```
T_22RV17: GNA11 p.Q209L (37.9%)
```

This means:
- 37.9% of reads have the mutant allele
- This is a **clonal driver mutation** (present in all tumour cells)
- VAF ~38% suggests **~76% tumour purity** (38% × 2, accounting for heterozygous state)
- Remaining ~24% is normal stromal/immune cells

```
T_25RV08: BAP1 p.L49P (90.5%)
```

This very high VAF indicates:
- Either very high tumour purity (>90%), OR
- Tumour has **monosomy 3** (lost the other chromosome 3), so only one BAP1 allele remains and it's mutant

## Typical UM Mutation Patterns

**Classic high-risk profile:**
- GNAQ or GNA11 mutation (Q209L/P, R183C)
- BAP1 loss-of-function mutation
- Monosomy 3 (visible in CNV data)
- 8q gain (visible in CNV data)

**Lower-risk alternative:**
- GNAQ or GNA11 mutation
- SF3B1 or EIF1AX mutation (NOT BAP1)
- Disomy 3

## File Structure

Expected Novogene delivery structure:
```
Novogene_Delivery/
└── 03.Result_*/
    └── result/
        └── Somatic/
            ├── Somatic_SNV/
            │   └── Annotation/
            │       ├── T_25RV04.muTect2.somatic.snv.annovar.hg38_multianno.xls.gz
            │       └── ...
            └── Somatic_INDEL/
                └── Annotation/
                    ├── T_25RV04.muTect2.somatic.indel.annovar.hg38_multianno.xls.gz
                    └── ...
```

## Example Workflow

```bash
# 1. Organize Novogene files
mkdir -p UM_annotations/SNV UM_annotations/INDEL

# 2. Copy annotation files
cp /path/to/Novogene/*/Somatic_SNV/Annotation/*.xls.gz UM_annotations/SNV/
cp /path/to/Novogene/*/Somatic_INDEL/Annotation/*.xls.gz UM_annotations/INDEL/

# 3. Run extraction
python extract_um_mutations.py \
  --snv-dir UM_annotations/SNV \
  --indel-dir UM_annotations/INDEL \
  --output-prefix my_cohort

# 4. Open results in Excel
open my_cohort_mutation_matrix.tsv
```

## Troubleshooting

**No mutations found**

Check that:
1. Files are from the `Somatic` folder (not `Mutation` folder - that's germline)
2. File names contain `somatic.snv` or `somatic.indel`
3. Files are the ANNOVAR annotation files (`.annovar.hg38_multianno.xls.gz`)

**Missing samples**

The script only outputs samples that have at least one mutation in the 7 UM genes. Samples with no driver mutations won't appear in the matrix. Check the detailed output for a complete list of processed samples.

## Citation

If you use this tool in your research, please cite:

```
Owens S. (2026). Uveal Melanoma Driver Mutation Extractor. 
GitHub: https://github.com/[your-username]/um-mutation-extractor
```

## License

MIT License - see LICENSE file

## Author

Sally Owens  
Royal College of Surgeons in Ireland / Dublin City University  
Contact: [your-email]

## Acknowledgments

- Sequencing performed by Novogene
- Analysis pipeline developed for uveal melanoma genomics research

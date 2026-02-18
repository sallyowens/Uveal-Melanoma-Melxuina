# Chromosomal Aberrations Extractor for Uveal Melanoma

Extract key chromosomal changes from CNVkit and FACETS copy number analysis results.

## Overview

This tool processes CNVkit `.cns` files (and optionally FACETS segments files) to create a summary table of chromosomal aberrations important in uveal melanoma prognosis:

- **Monosomy 3** - Loss of chromosome 3 (strongest metastatic predictor)
- **8q gain** - Gain of chromosome 8 long arm (associated with metastasis)
- **6p gain** - Additional risk marker
- **1p loss** - Poor prognosis marker
- **Tumor purity and ploidy** (from FACETS if available)

## Requirements

```bash
pip install pandas numpy
```

Python 3.7+

## Usage

### With CNVkit only

```bash
python extract_chromosomal_aberrations.py \
  --cnvkit-dir /path/to/cnvkit/results \
  --output chr_aberrations.tsv
```

### With CNVkit + FACETS

```bash
python extract_chromosomal_aberrations.py \
  --cnvkit-dir /path/to/cnvkit/results \
  --facets-dir /path/to/facets/results \
  --output chr_aberrations.tsv
```

## Input Files

**CNVkit files** (required):
```
T_25RV04.cns
T_25RV05.cns
...
```

**FACETS files** (optional):
```
segments_T_25RV04.csv
segments_T_25RV05.csv
...
```

## Output Format

Tab-separated table with one row per sample:

| Sample | Chr3_Status | Chr3_Mean_CN | Chr3_Pct_Loss | Chr8_Status | Chr8_Mean_CN | Chr8_Pct_Gain | Chr6_Status | Chr6_Mean_CN | Chr1_Status | Chr1_Mean_CN | Purity | Ploidy |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T_25RV04 | Monosomy | 1.2 | 95.3 | Gain | 3.2 | 68.5 | Normal | 2.1 | Normal | 2.0 | 82.5 | 2.1 |
| T_25RV05 | Disomy | 2.0 | 8.2 | Normal | 2.1 | 5.3 | Normal | 2.0 | Normal | 2.0 | 76.3 | 2.0 |

## Understanding the Classifications

### Chromosome 3 Status

- **Monosomy**: Mean CN < 1.3 AND >80% of chromosome lost
- **Disomy** (normal): Mean CN > 1.7 AND <20% loss
- **Partial_Loss**: Some loss but doesn't meet monosomy criteria

**Clinical significance**: Monosomy 3 is the strongest predictor of metastatic disease in UM.

### Chromosome 8 Status

- **Gain**: Mean CN > 2.5 OR >50% of chromosome gained
- **Normal**: Otherwise

**Clinical significance**: 8q gain combined with monosomy 3 indicates very high metastatic risk.

### Chromosome 6p Status

- **Gain**: Mean CN > 2.5 OR >30% gained
- **Normal**: Otherwise

**Clinical significance**: 6p gain is seen more frequently in disomy 3 tumors and associated with favorable prognosis when isolated.

### Chromosome 1p Status

- **Loss**: >30% of chromosome lost
- **Normal**: Otherwise

**Clinical significance**: 1p loss associated with poor prognosis.

## Typical UM Prognostic Profiles

**High-risk (Class 2):**
```
Chr3_Status: Monosomy
Chr8_Status: Gain
Expected metastasis rate: 70-90% by 5 years
```

**Low-risk (Class 1):**
```
Chr3_Status: Disomy
Chr8_Status: Normal
6p may show Gain
Expected metastasis rate: <5% by 5 years
```

**Intermediate (Class 1B):**
```
Chr3_Status: Partial_Loss or Disomy
Chr8_Status: Normal or Gain
SF3B1 or EIF1AX mutation often present
```

## Purity and Ploidy (FACETS)

**Purity**: Percentage of tumor cells in the sample
- >70%: High purity (reliable CNV calls)
- 40-70%: Moderate purity
- <40%: Low purity (CNV calls less reliable)

**Ploidy**: Average number of chromosome copies across the genome
- ~2.0: Diploid (normal)
- >2.5: Polyploid
- <1.8: Hypodiploid

## Example Workflow

```bash
# If you have your MeluXina results organized by group:
cd /Volumes/Expansion/UVM_WGS_Analysis

# Process all groups at once
python extract_chromosomal_aberrations.py \
  --cnvkit-dir grp1_results/cnvkit \
  --facets-dir grp1_results/FACETS/results \
  --output grp1_chr_aberrations.tsv

# Or combine all groups by first collecting all .cns files in one place:
mkdir -p all_cnvkit_results
cp grp*/cnvkit/*.cns all_cnvkit_results/

python extract_chromosomal_aberrations.py \
  --cnvkit-dir all_cnvkit_results \
  --output all_samples_chr_aberrations.tsv
```

## Combining with Mutation Data

Merge this chromosomal data with your mutation matrix in Excel/R:

```R
# In R
chr_data <- read.delim("chr_aberrations.tsv")
mut_data <- read.delim("UM_mutation_matrix.tsv")

combined <- merge(chr_data, mut_data, by="Sample")

# Cox regression for survival
library(survival)
coxph(Surv(OS_months, OS_status) ~ 
        Chr3_Status + BAP1_Mutation + Chr8_Status, 
      data=combined)
```

## Troubleshooting

**No .cns files found**

Make sure you're pointing to a directory that contains `.cns` files directly, not a parent directory. CNVkit typically outputs files like:
- `sample.cns` (segmented calls)
- `sample.call.cns` (thresholded calls)  
- `sample.revised.call.cns` (revised calls)

The script will process any of these.

**FACETS purity/ploidy not working**

FACETS output formats can vary. If purity/ploidy aren't being extracted, check your FACETS segments file format. You may need to adjust the `get_facets_data()` function to match your specific FACETS output. Alternatively, run without `--facets-dir` and add purity/ploidy manually later.

## Citation

If you use this tool, please cite:

```
Owens S. (2026). Chromosomal Aberrations Extractor for Uveal Melanoma.
GitHub: https://github.com/[your-username]/um-cnv-extractor
```

## License

MIT License

## Author

Sally Owens  
Royal College of Surgeons in Ireland / Dublin City University

## Acknowledgments

- CNV analysis performed using CNVkit and FACETS
- Pipeline developed for uveal melanoma genomics research

# Whole Genome Sequencing Analysis Pipeline for Uveal Melanoma

Comprehensive bioinformatics pipeline for copy number variation (CNV) analysis and tumor characterization from paired tumor-normal whole genome sequencing data.

## Overview

This pipeline processes WGS BAM files to perform:
- Quality control assessment
- Coverage analysis
- Copy number variation detection (CNVkit)
- Tumor purity and ploidy estimation (FACETS)

**Institution**: MeluXina HPC  
**Reference Genome**: GRCh38 (hg38)  
**Sample Type**: Paired tumor-normal WGS samples

---

## Table of Contents

- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Input Data](#input-data)
- [Usage](#usage)
- [Output Files](#output-files)
- [Workflow Details](#workflow-details)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contact](#contact)

---

## Pipeline Overview
```
BAM Files (*.final.bam)
    ↓
[1. BAM QC] → Quality metrics, coverage statistics, duplication rates
    ↓
[2. Coverage Analysis] → Genome-wide depth analysis
    ↓
[3. CNVkit] → Copy number variation detection
    ↓
[4. SNP Pileup] → Extract SNP information for FACETS
    ↓
[5. FACETS] → Tumor purity, ploidy, and refined CNV calls
```

---

## Requirements

### Software Dependencies
- **Samtools** v1.17
- **Picard** v2.27.5
- **FastQC** v0.12.1
- **MultiQC** v1.14
- **CNVkit** v0.9.10
- **snp-pileup** (from FACETS suite)
- **BWA** v0.7.18
- **Cutadapt** v5.2
- **R** v4.2.3
  - facets v0.6.2
  - ggplot2, dplyr, tidyr, readr

### Compute Resources
- **HPC System**: SLURM job scheduler
- **Storage**: ~500 GB for intermediate and output files
- **Reference Genome**: GRCh38 (~3 GB)

---

## Installation

### 1. Clone Repository
```bash
git clone https://github.com/YOUR_USERNAME/wgs-uvm-analysis.git
cd wgs-uvm-analysis
```

### 2. Set Up Conda Environment
```bash
# Install Miniconda (if not already installed)
cd ~/Project
mkdir -p software
cd software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/Project/software/miniconda3
rm Miniconda3-latest-Linux-x86_64.sh

# Initialize conda
~/Project/software/miniconda3/bin/conda init bash
source ~/.bashrc
```

### 3. Create Analysis Environment
```bash
# Create environment with all tools
conda create -n wgs_analysis python=3.10 -y
conda activate wgs_analysis

# Install bioinformatics tools
conda install -c bioconda -c conda-forge \
    samtools=1.17 \
    picard=2.27.5 \
    fastqc=0.12.1 \
    multiqc=1.14 \
    cnvkit=0.9.10 \
    snp-pileup \
    htslib \
    bwa \
    cutadapt \
    r-base=4.2 \
    r-ggplot2 \
    r-dplyr \
    r-tidyr \
    r-readr \
    r-facets \
    -y
```

### 4. Download Reference Genome
```bash
mkdir -p ~/Project/reference
cd ~/Project/reference

# Download GRCh38 reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Download index
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
```

### 5. Create Directory Structure
```bash
mkdir -p ~/Project/{logs,scripts,bam_qc,cnvkit,coverage,FACETS/{VCF,pileups,results,plots,reports},tmp}
```

---

## Directory Structure
```
Project/
├── grp1/                          # Input BAM files
│   ├── T_24RV18.final.bam        # Tumor samples (T_*)
│   ├── B_24RV18.final.bam        # Normal samples (B_*)
│   └── *.bai                      # BAM index files
├── scripts/                       # Analysis scripts
│   ├── bam_qc_check.sh
│   ├── coverage_analysis.sh
│   ├── cnvkit_analysis.sh
│   ├── snp_pileup.sh
│   └── facets_analysis.R
├── logs/                          # SLURM job logs
├── reference/                     # Reference genome
│   └── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
├── bam_qc/                        # QC results
│   ├── stats/
│   └── multiqc_report/
├── coverage/                      # Coverage analysis
│   ├── coverage_stats/
│   ├── depth_files/
│   └── plots_data/
├── cnvkit/                        # CNVkit results
│   ├── raw/
│   └── revised/
└── FACETS/                        # FACETS analysis
    ├── VCF/
    ├── pileups/
    ├── results/
    ├── plots/
    └── reports/
```

---

## Input Data

### File Naming Convention
- **Tumor samples**: `T_[SAMPLE_ID].final.bam`
- **Normal samples**: `B_[SAMPLE_ID].final.bam`
- Index files: `*.final.bam.bai`

### Required File Properties
- Coordinate-sorted BAM files
- Read groups present
- Duplicates marked/removed
- hg38/GRCh38 alignment

### Sample Data Structure
```
grp1/
├── T_24RV18.final.bam (53 GB)  →  B_24RV18.final.bam (25 GB)
├── T_24RV21.final.bam (58 GB)  →  B_24RV21.final.bam (30 GB)
├── T_24RV22.final.bam (55 GB)  →  B_24RV22.final.bam (26 GB)
└── T_24RV24.final.bam (45 GB)  →  B_24RV24.final.bam (27 GB)
```

---

## Usage

### Update SLURM Account in Scripts

Before running, update the account in all scripts:
```bash
cd ~/Project/scripts

# Replace p201093 with your account number in all scripts
sed -i 's/p201093/YOUR_ACCOUNT/g' *.sh
```

### Step 1: BAM Quality Control
```bash
cd ~/Project/scripts
sbatch bam_qc_check.sh

# Monitor progress
squeue -u $USER
tail -f ~/Project/logs/bam_qc_*.out
```

**Runtime**: ~2-4 hours  
**Output**: QC metrics, flagstat, coverage statistics, MultiQC report

### Step 2: Coverage Analysis (Optional)
```bash
sbatch coverage_analysis.sh
```

**Runtime**: ~6-12 hours  
**Output**: Genome-wide depth files, coverage statistics, sample metadata

### Step 3: CNVkit Analysis
```bash
sbatch cnvkit_analysis.sh
```

**Runtime**: ~12-24 hours  
**Output**: Copy number segments, scatter plots, chromosome diagrams

### Step 4: SNP Pileup for FACETS
```bash
sbatch snp_pileup.sh
```

**Runtime**: ~24-48 hours  
**Output**: Sorted pileup files for FACETS

### Step 5: FACETS Analysis
```bash
sbatch --wrap="Rscript ~/Project/scripts/facets_analysis.R" \
       --job-name=facets_R \
       --account=p201093 \
       --partition=cpu \
       --qos=default \
       --output=~/Project/logs/facets_R_%j.out \
       --error=~/Project/logs/facets_R_%j.err \
       --cpus-per-task=4 \
       --mem=32G \
       --time=06:00:00
```

**Runtime**: ~2-6 hours  
**Output**: Purity/ploidy estimates, refined CNV calls, FACETS plots

---

## Output Files

### BAM QC
- `BAM_QC_Summary_Report.txt` - Comprehensive QC summary
- `multiqc_report.html` - Interactive QC dashboard
- `*_flagstat.txt` - Alignment statistics
- `*_duplication_metrics.txt` - Duplication rates

### Coverage Analysis
- `all_samples_coverage.txt` - Combined coverage data
- `sample_metadata.txt` - Sample type classification
- `*_depth.txt.gz` - Compressed depth files
- `coverage_summary.txt` - Analysis summary

### CNVkit
- `*.revised.call.cns` - Final copy number segments
- `*.scatter.pdf` - Genome-wide CNV plots
- `*.diagram.pdf` - Chromosome-level diagrams
- `cnvkit_summary.txt` - Analysis parameters

### FACETS
- `segments_[SAMPLE].csv` - Copy number segments with purity/ploidy
- `facets_[SAMPLE].png` - Main FACETS plot
- `spider_[SAMPLE].png` - Purity-ploidy spider plot
- `facets_summary_all_samples.csv` - Combined results table
- `fit_[SAMPLE].rds` - R object for reanalysis

---

## Workflow Details

### 1. BAM Quality Control
- Validates BAM integrity
- Checks for read groups
- Calculates alignment statistics
- Assesses duplication rates
- Generates per-chromosome coverage
- Creates MultiQC aggregate report

**Key Metrics**:
- Mapping rate (expect >95%)
- Duplication rate (expect <30%)
- Mean coverage (expect 30-60x for WGS)

### 2. Coverage Analysis
- Generates genome-wide depth profiles
- Creates coverage histograms
- Calculates per-chromosome statistics
- Classifies samples as tumor/normal

**Output Format**: Tab-delimited files compatible with R/Python

### 3. CNVkit - Copy Number Detection
**Method**: WGS mode (bin-based approach)

**Parameters**:
- Bin size: Automatically determined for WGS
- Segmentation: CBS algorithm
- Copy number calling: Clonal method

**Workflow**:
1. Create genome access file
2. Batch process tumor-normal pairs
3. Segment copy number changes
4. Call integer copy numbers
5. Generate visualizations

### 4. SNP Pileup
**Purpose**: Extract allele-specific read counts for FACETS

**Parameters**:
- Minimum base quality: 25
- Minimum mapping quality: 20
- Minimum read depth: 10

**VCF Source**: dbSNP common variants (GRCh38)

### 5. FACETS - Purity & Ploidy
**Purpose**: Estimate tumor purity, ploidy, and refine CNV calls

**Parameters**:
- Critical value (cval): 50 (more stringent for WGS)
- Minimum depth: 25
- Minimum heterozygous SNPs: 15
- EM iterations: 30

**Outputs**:
- Tumor purity (fraction of tumor cells)
- Tumor ploidy (average DNA content)
- Allele-specific copy number
- Loss of heterozygosity (LOH) regions

---

## Troubleshooting

### Common Issues

#### Job Submission Fails
```bash
# Error: Invalid account or partition
# Fix: Check your account
sacctmgr show associations user=$USER format=account,partition

# Update scripts with correct account
sed -i 's/p201093/YOUR_ACCOUNT/g' ~/Project/scripts/*.sh
```

#### Out of Memory
```bash
# Increase memory in SLURM header
#SBATCH --mem=128G  # or higher
```

#### Missing Conda Environment
```bash
# Reactivate environment
conda activate wgs_analysis

# If missing, recreate
conda env list
conda create -n wgs_analysis python=3.10 -y
```

#### SNP Pileup Takes Too Long
```bash
# Check if running
squeue -u $USER

# Increase time limit
#SBATCH --time=72:00:00
```

#### FACETS Fails
```bash
# Check pileup files exist
ls -lh ~/Project/FACETS/pileups/*.sorted.pileup

# Check R packages
R -e "library(facets); packageVersion('facets')"
```

---

## Performance Benchmarks

| Analysis Step | Runtime | Memory | CPUs | Output Size |
|--------------|---------|--------|------|-------------|
| BAM QC | 2-4 hrs | 32 GB | 8 | ~500 MB |
| Coverage | 6-12 hrs | 64 GB | 8 | ~50-100 GB |
| CNVkit | 12-24 hrs | 128 GB | 16 | ~2-5 GB |
| SNP Pileup | 24-48 hrs | 64 GB | 4 | ~20-40 GB |
| FACETS | 2-6 hrs | 32 GB | 4 | ~2 GB |

*Based on 4 tumor-normal pairs, ~40-60x coverage*

---

## Key Publications

### Methods
- **CNVkit**: Talevich et al. (2016) *PLOS Comput Biol* - [DOI: 10.1371/journal.pcbi.1004873](https://doi.org/10.1371/journal.pcbi.1004873)
- **FACETS**: Shen & Seshan (2016) *Nucleic Acids Res* - [DOI: 10.1093/nar/gkw520](https://doi.org/10.1093/nar/gkw520)

### Reference
- **GRCh38**: Genome Reference Consortium - [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)

---

## Repository Structure
```
.
├── README.md                     # This file
├── scripts/
│   ├── bam_qc_check.sh          # BAM quality control
│   ├── coverage_analysis.sh     # Coverage analysis
│   ├── cnvkit_analysis.sh       # CNVkit CNV detection
│   ├── snp_pileup.sh            # FACETS preprocessing
│   └── facets_analysis.R        # FACETS purity/ploidy
├── docs/
│   ├── installation.md          # Detailed installation
│   ├── parameters.md            # Parameter descriptions
│   └── interpretation.md        # Result interpretation
└── LICENSE                       # MIT License
```

---

## Citation

If you use this pipeline, please cite:
```bibtex
@software{wgs_uvm_pipeline,
  author = {Your Name},
  title = {WGS Analysis Pipeline for Uveal Melanoma},
  year = {2025},
  url = {https://github.com/YOUR_USERNAME/wgs-uvm-analysis}
}
```

And cite the individual tools used (CNVkit, FACETS, etc.)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Researcher**: Sally Owens  
**Institution**: Dublin City University  
**Email**: sally.owens7@mail.dcu.ie 
**Project**: Uveal Melanoma WGS Analysis  
**HPC**: MeluXina

---

## Acknowledgments

- MeluXina HPC for computational resources
- CNVkit and FACETS development teams
- The Bioconda community

---

## Version History

- **v1.0.0** (2025-01-XX) - Initial release
  - BAM QC pipeline
  - CNVkit WGS analysis
  - FACETS integration
  - 4 paired samples processed

---

## Future Development

- [ ] Somatic variant calling integration
- [ ] Structural variant detection
- [ ] Automated report generation
- [ ] Integration with additional CNV callers
- [ ] Visualization dashboard

---

*Last updated: January 2025*

# WGS Analysis Pipeline for Meluxina

Adapted scripts for whole genome sequencing (WGS) analysis on Meluxina HPC, starting from BAM/BAI files.

## 📋 Overview

This pipeline processes matched tumor-normal WGS BAM files (T_*.final.bam and B_*.final.bam format) through:
1. **Quality Control** - Assess BAM file quality
2. **BAM Processing** - Add read groups, filter, deduplicate
3. **Coverage Analysis** - Calculate genome-wide coverage metrics
4. **CNVkit** - Copy number variation calling
5. **FACETS** - Purity and ploidy estimation, CNV analysis

## 📁 Your Current Setup

**Working Directory:** `/home/users/u103499/Project`

**Current BAM files (grp1):**
- Tumor files: `T_24RV18.final.bam`, `T_24RV21.final.bam`, `T_24RV22.final.bam`, `T_24RV24.final.bam`
- Normal files: `B_24RV18.final.bam`, `B_24RV21.final.bam`, `B_24RV22.final.bam`, `B_24RV24.final.bam`

**File naming:** T_ prefix = tumor, B_ prefix = blood/normal

## 🚀 Quick Start

### 1. Setup Directory Structure

```bash
# Navigate to your Project directory
cd /home/users/u103499/Project

# Create necessary directories
mkdir -p scripts
mkdir -p logs
mkdir -p reference
mkdir -p processed_bams
mkdir -p bam_qc
mkdir -p coverage
mkdir -p CNVkit
mkdir -p FACETS

# Copy all scripts to scripts directory
cd scripts
# Upload the 6 script files here

# Make scripts executable
chmod +x *.sh
```

### 2. Download Reference Genome (One-time Setup)

```bash
cd /home/users/u103499/Project/reference

# Download GRCh38 reference (same as in your original pipeline)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Download index file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

# Index with samtools (after loading SAMtools module)
module load SAMtools
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

### 3. Check Available Modules on Meluxina

```bash
# Check what modules are available
module avail SAMtools
module avail Picard
module avail Python
module avail R

# Load modules to test
module load SAMtools
module load Picard

# Check if tools work
samtools --version
picard -h
```

### 4. Submit Your First Job (BAM QC)

```bash
cd /home/users/u103499/Project/scripts

# Submit BAM quality check
sbatch bam_qc.sh

# Check job status
squeue -u u103499

# View log in real-time
tail -f ../logs/bam_qc_*.out
```

### 5. Process BAMs After QC

```bash
# Submit BAM processing (will process all 4 pairs)
sbatch bam_processing.sh

# Monitor progress
watch -n 30 'squeue -u u103499'
```

### 6. Run Complete Pipeline

```bash
# Option A: Submit all steps with dependencies (automatic)
# This ensures each step waits for the previous one to complete

JOB_QC=$(sbatch bam_qc.sh | awk '{print $4}')
JOB_PROC=$(sbatch --dependency=afterok:${JOB_QC} bam_processing.sh | awk '{print $4}')
JOB_COV=$(sbatch --dependency=afterok:${JOB_PROC} coverage.sh | awk '{print $4}')
JOB_CNV=$(sbatch --dependency=afterok:${JOB_PROC} cnvkit.sh | awk '{print $4}')
JOB_FACETS=$(sbatch --dependency=afterok:${JOB_PROC} snp_pileup.sh | awk '{print $4}')

echo "Jobs submitted:"
echo "QC: $JOB_QC"
echo "Processing: $JOB_PROC (waits for QC)"
echo "Coverage: $JOB_COV (waits for Processing)"
echo "CNVkit: $JOB_CNV (waits for Processing)"
echo "FACETS: $JOB_FACETS (waits for Processing)"
```

## 📂 Directory Structure After Setup

```
/home/users/u103499/Project/
├── grp1/                          # Your original BAM files
│   ├── T_24RV18.final.bam
│   ├── T_24RV18.final.bam.bai
│   ├── B_24RV18.final.bam
│   ├── B_24RV18.final.bam.bai
│   └── ... (other pairs)
├── scripts/                       # All SLURM scripts
│   ├── bam_qc.sh
│   ├── bam_processing.sh
│   ├── coverage.sh
│   ├── cnvkit.sh
│   └── snp_pileup.sh
├── logs/                          # SLURM output logs
├── reference/                     # Reference genome
│   └── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
├── processed_bams/               # Processed BAM files
│   ├── final_bams/
│   │   ├── 24RV18_tumor_final.bam
│   │   ├── 24RV18_normal_final.bam
│   │   └── ... (other pairs)
│   ├── metrics/
│   └── intermediate_files/
├── bam_qc/                       # QC results
│   ├── stats/
│   ├── metrics/
│   ├── plots/
│   └── BAM_QC_Report.txt
├── coverage/                     # Coverage analysis
│   ├── coverage_stats/
│   ├── plots_data/
│   └── combined_coverage_summary.txt
├── CNVkit/                       # CNVkit results
│   ├── bedfiles/
│   └── results/
│       └── revised/              # Final CNV calls
└── FACETS/                       # FACETS results
    ├── VCF/
    ├── pileups/
    └── logs/
```

## 🔧 Script Details

### 1. bam_qc.sh
**Purpose:** Quality control for your BAM files before analysis

**What it checks:**
- Mapping statistics
- Duplicate rates
- Coverage per chromosome
- Insert size distributions
- GC bias
- Alignment quality

**Resources:**
- Time: 12 hours
- CPUs: 8
- Memory: 32GB

**Input:** Files in `/home/users/u103499/Project/grp1/*.final.bam`

**Output:** 
- Summary report: `/home/users/u103499/Project/bam_qc/BAM_QC_Report.txt`
- CSV for R: `/home/users/u103499/Project/bam_qc/qc_metrics.csv`

### 2. bam_processing.sh
**Purpose:** Process both tumor and normal BAMs for each patient

**Processing steps:**
1. Adds read groups (RGID, RGSM, etc.)
2. Filters low quality reads (MAPQ < 10)
3. Removes supplementary, secondary, and unmapped reads
4. Removes duplicates with Picard
5. Indexes final BAMs

**Resources:**
- Time: 48 hours
- CPUs: 16
- Memory: 64GB
- Array: 0-3 (processes all 4 pairs)

**Input:** 
- Tumor: `T_24RV18.final.bam`, etc.
- Normal: `B_24RV18.final.bam`, etc.

**Output:**
- `24RV18_tumor_final.bam` and `24RV18_normal_final.bam`
- Located in: `/home/users/u103499/Project/processed_bams/final_bams/`

### 3. coverage.sh
**Purpose:** Calculate genome-wide coverage

**Resources:**
- Time: 24 hours
- CPUs: 8
- Memory: 32GB

**Output:**
- Per-sample stats: `coverage/coverage_stats/`
- R-ready data: `coverage/plots_data/all_samples_coverage.txt`
- Summary: `coverage/combined_coverage_summary.txt`

### 4. cnvkit.sh
**Purpose:** Copy number variation calling

**Resources:**
- Time: 48 hours
- CPUs: 16
- Memory: 128GB

**Output:**
- `CNVkit/results/*.cnr` - Copy number ratios
- `CNVkit/results/revised/*.revised.call.cns` - Final CNV calls
- `CNVkit/results/revised/*.scatter.png` - Visualization
- `CNVkit/results/revised/*.seg` - For IGV

**Two phases:**
1. Initial calling
2. Refined with purity normalization

### 5. snp_pileup.sh
**Purpose:** Generate pileups for FACETS purity/ploidy analysis

**Resources:**
- Time: 72 hours (WGS files are large!)
- CPUs: 4
- Memory: 64GB
- Array: 0-3 (one per patient pair)

**Output:**
- `FACETS/pileups/24RV18.sorted.pileup`
- One pileup file per patient for R analysis

## 📊 Monitoring Your Jobs

```bash
# Check all your jobs
squeue -u u103499

# Detailed job info
scontrol show job <JOB_ID>

# View log in real-time
tail -f logs/bam_qc_<JOB_ID>.out

# Check job efficiency after completion
seff <JOB_ID>

# Cancel a job
scancel <JOB_ID>

# Cancel all your jobs
scancel -u u103499
```

## ⏱️ Expected Runtimes for Your 4 Pairs

Based on your file sizes (703MB to 47GB per file):

| Step | Total Time |
|------|-----------|
| BAM QC | 6-12 hours |
| BAM Processing | 24-48 hours |
| Coverage | 12-24 hours |
| CNVkit | 24-48 hours |
| SNP Pileup | 48-72 hours |

**Total pipeline time:** ~5-8 days for complete analysis of grp1

## 💾 Your Meluxina Quotas

Current usage:
- CPU node-hours: 0/430 per month
- Home directory: 2/100 GiB (keep scripts here)
- Project p201093: 314/1000 GiB (use for data)
- Project p200272: 3064/5120 GiB

**Recommendations:**
- Store BAM files in `/project/home/p201093/` if needed
- Monitor with: `myquota`
- Clean intermediate files after each batch

## 🔄 Processing Multiple Batches

You have 40 samples total in batches of 4 pairs:

### Batch 1 (grp1) - Current
Already set up: T_24RV18, T_24RV21, T_24RV22, T_24RV24

### For Next Batches (grp2, grp3, etc.)

```bash
# Create new group directory
mkdir -p /home/users/u103499/Project/grp2

# Upload your next 4 pairs to grp2
# Then update BAM_DIR in bam_qc.sh:
BAM_DIR=/home/users/u103499/Project/grp2

# Resubmit the pipeline
sbatch bam_qc.sh
# ... etc
```

**Alternative:** Create a batch submission script:
```bash
#!/bin/bash
for grp in grp1 grp2 grp3 grp4 grp5 grp6 grp7 grp8 grp9 grp10; do
    # Update BAM_DIR in each script
    sed -i "s|BAM_DIR=.*|BAM_DIR=/home/users/u103499/Project/$grp|" bam_qc.sh
    
    # Submit with dependency
    sbatch bam_qc.sh
    sleep 60  # Wait between submissions
done
```

## ⚠️ Common Issues & Solutions

### Issue: "Module not found"
```bash
# Check available modules
module spider SAMtools
module avail 2>&1 | grep -i samtools

# Try different module names
module load SAMtools/1.18
module load samtools
```

### Issue: "Permission denied"
```bash
# Make scripts executable
chmod +x *.sh

# Check directory permissions
ls -la /home/users/u103499/Project/
```

### Issue: "No space left on device"
```bash
# Check your quota
myquota

# Clean up intermediate files
rm -rf processed_bams/intermediate_files/*

# Remove old logs
rm -f logs/*.out logs/*.err
```

### Issue: "Job failed immediately"
```bash
# Check the error log
cat logs/bam_qc_<JOB_ID>.err

# Common causes:
# - Module not loaded
# - Wrong path
# - Missing reference file
```

### Issue: Reference genome not found
```bash
# Check if reference exists
ls -lh /home/users/u103499/Project/reference/

# Re-download if needed
cd /home/users/u103499/Project/reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip *.gz
```

## 🎯 After Pipeline Completes

### 1. Review QC Report
```bash
cat /home/users/u103499/Project/bam_qc/BAM_QC_Report.txt
```

Check for:
- Mapping rate > 95%
- Duplicate rate < 30%
- Mean coverage appropriate for WGS (typically 30-60x)

### 2. Visualize Coverage (Run Locally in RStudio)
```bash
# Download coverage data
scp u103499@login.lxp.lu:/home/users/u103499/Project/coverage/plots_data/* ./

# Run your coverage.R script with the data
```

### 3. Review CNVkit Results
```bash
# View scatter plots
ls CNVkit/results/revised/*.scatter.png

# Check segment files for IGV
ls CNVkit/results/revised/*.seg
```

### 4. Run FACETS R Analysis
```bash
# Download pileup files
scp u103499@login.lxp.lu:/home/users/u103499/Project/FACETS/pileups/*.sorted.pileup ./

# Run FACETS.R script locally (provided in your original scripts)
```

### 5. Compare Results
- CNVkit vs FACETS CNV calls
- Use FACETS purity estimates to refine CNVkit (rerun with --purity flag)

## 📝 Quick Command Reference

```bash
# Submit job
sbatch script.sh

# Check status
squeue -u u103499

# View output
tail -f logs/job_*.out

# Cancel job
scancel <JOB_ID>

# Check quota
myquota

# Load modules
module load SAMtools

# Check job efficiency
seff <JOB_ID>

# Download results
scp u103499@login.lxp.lu:/path/to/file ./
```

## ✅ Pre-Flight Checklist

Before running the pipeline:

- [ ] All scripts uploaded to `/home/users/u103499/Project/scripts/`
- [ ] Scripts made executable (`chmod +x *.sh`)
- [ ] Reference genome downloaded to `/home/users/u103499/Project/reference/`
- [ ] BAM files in `/home/users/u103499/Project/grp1/`
- [ ] All BAM files have corresponding .bai index files
- [ ] Logs directory created
- [ ] Modules available (check with `module avail`)
- [ ] Sufficient disk quota (check with `myquota`)

## 🆘 Getting Help

**Meluxina Documentation:** https://docs.lxp.lu/

**Support Email:** support@lxp.lu

**Check job logs:** Always check both `.out` and `.err` files in the logs directory

**Test first:** Run BAM QC on just one sample to verify everything works before processing all 4 pairs

## 📚 Next Steps

1. ✅ Run BAM QC first - this is fast and will catch any setup issues
2. ✅ Review QC report before proceeding
3. ✅ Process BAMs - this takes longest
4. ✅ Run coverage and CNVkit in parallel (both depend on processed BAMs)
5. ✅ Generate FACETS pileups (can run while CNVkit is running)
6. ✅ Download results and run R analyses locally

Good luck with your analysis! 🚀

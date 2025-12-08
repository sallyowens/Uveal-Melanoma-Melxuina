# WGS Analysis Pipeline for Meluxina HPC

Whole genome sequencing (WGS) analysis pipeline adapted for the Meluxina supercomputer in Luxembourg. This pipeline processes matched tumor-normal BAM files for CNV analysis and somatic variant calling.

## 📋 Overview

This pipeline processes pre-aligned WGS BAM files (T_*.final.bam for tumor, B_*.final.bam for normal) through:

1. **Quality Control** - Comprehensive BAM file quality assessment
2. **BAM Processing** - Add read groups, filter low-quality reads, remove duplicates
3. **Coverage Analysis** - Genome-wide coverage statistics and metrics
4. **CNVkit** - Copy number variation detection and calling
5. **FACETS** - Tumor purity, ploidy estimation, and CNV analysis

## 🖥️ System Requirements

**Meluxina HPC Cluster**
- Account with project allocation (e.g., p201093)
- Access to CPU partition
- Sufficient quota for large BAM files (703MB - 47GB per file)

**Software (via modules)**
- samtools/1.20-foss-2023a
- picard/3.2
- Java/21.0.5
- Python/3.12.3-GCCcore-13.3.0
- R/4.4.2-gfbf-2024a

## 📁 Directory Structure

```
/home/users/u103499/Project/
├── grp1/                          # Original BAM files (T_ and B_ prefixed)
│   ├── T_24RV18.final.bam
│   ├── T_24RV18.final.bam.bai
│   ├── B_24RV18.final.bam
│   ├── B_24RV18.final.bam.bai
│   └── ... (additional sample pairs)
├── scripts/                       # SLURM job scripts
│   ├── bam_qc.sh
│   ├── bam_processing.sh
│   ├── coverage.sh
│   ├── cnvkit.sh
│   └── snp_pileup.sh
├── logs/                          # SLURM output/error logs
├── reference/                     # Reference genome (GRCh38)
│   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
│   └── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
├── bam_qc/                       # QC results
├── processed_bams/               # Processed BAM files
│   ├── final_bams/
│   │   ├── 24RV18_tumor_final.bam
│   │   ├── 24RV18_normal_final.bam
│   │   └── ...
│   └── metrics/
├── coverage/                     # Coverage analysis results
├── CNVkit/                       # CNVkit CNV results
└── FACETS/                       # FACETS purity/ploidy results
```

## 🚀 Quick Start

### 1. Initial Setup

```bash
# Navigate to project directory
cd /home/users/u103499/Project

# Create necessary directories
mkdir -p scripts logs reference

# Upload your BAM files to grp1/
# Upload scripts to scripts/
# Make scripts executable
cd scripts
chmod +x *.sh
```

### 2. Download Reference Genome (One-time)

```bash
cd /home/users/u103499/Project/reference

# Download GRCh38 reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Download index
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
```

### 3. Submit Jobs

```bash
cd /home/users/u103499/Project/scripts

# Step 1: Quality Control (12 hours)
sbatch bam_qc.sh

# Step 2: Process BAMs (48 hours) - after QC completes
sbatch bam_processing.sh

# Step 3: Coverage Analysis (24 hours) - after processing
sbatch coverage.sh

# Step 4: CNVkit (48 hours) - after processing
sbatch cnvkit.sh

# Step 5: FACETS SNP Pileup (72 hours) - after processing
sbatch snp_pileup.sh
```

### 4. Monitor Jobs

```bash
# Check job status
squeue -u u103499

# View output log (replace JOB_ID)
tail -f logs/bam_qc_JOB_ID.out

# Check job history
sacct -u u103499 --format=JobID,JobName,State,ExitCode,Elapsed
```

## 🔧 Scripts Overview

### bam_qc.sh
**Purpose:** Comprehensive quality control for BAM files

**Resources:**
- Time: 12 hours
- CPUs: 8
- Memory: 32GB

**Metrics Collected:**
- Mapping statistics (flagstat)
- Per-base statistics
- Coverage per chromosome
- Insert size distribution
- Alignment summary
- GC bias

**Outputs:**
- `bam_qc/BAM_QC_Report.txt` - Summary report
- `bam_qc/qc_metrics.csv` - CSV for R analysis
- `bam_qc/stats/` - Individual sample stats
- `bam_qc/metrics/` - Picard metrics
- `bam_qc/plots/` - QC plots (PDF)

---

### bam_processing.sh
**Purpose:** Add read groups, filter reads, remove duplicates

**Resources:**
- Time: 48 hours
- CPUs: 16
- Memory: 64GB
- Array: 0-3 (processes 4 pairs simultaneously)

**Processing Steps:**
1. Add read groups (RGID, RGSM, RGPL, etc.)
2. Filter supplementary alignments (flag 0x800)
3. Filter secondary alignments (flag 0x100)
4. Filter unmapped reads (flag 0x4)
5. Filter low MAPQ reads (<10)
6. Remove duplicates with Picard MarkDuplicates
7. Index final BAMs

**Outputs:**
- `processed_bams/final_bams/{patient_id}_tumor_final.bam`
- `processed_bams/final_bams/{patient_id}_normal_final.bam`
- `processed_bams/metrics/{patient_id}_{tumor|normal}_dup_metrics.txt`

---

### coverage.sh
**Purpose:** Calculate genome-wide coverage statistics

**Resources:**
- Time: 24 hours
- CPUs: 8
- Memory: 32GB

**Outputs:**
- `coverage/coverage_stats/` - Per-sample coverage
- `coverage/combined_coverage_summary.txt` - Overall summary
- `coverage/plots_data/all_samples_coverage.txt` - R-ready data
- `coverage/plots_data/sample_metadata.txt` - Sample classifications

---

### cnvkit.sh
**Purpose:** Copy number variation detection

**Resources:**
- Time: 48 hours
- CPUs: 16
- Memory: 128GB

**Two-Phase Analysis:**
1. **Initial calling:** Standard CNV detection
2. **Refined calling:** Purity normalization and segmentation refinement

**Outputs:**
- `CNVkit/results/*.cnr` - Copy number ratios
- `CNVkit/results/revised/*.revised.call.cns` - Final CNV calls
- `CNVkit/results/revised/*.scatter.png` - Visualization plots
- `CNVkit/results/revised/*.seg` - SEG format for IGV
- `CNVkit/results/revised/*.bed` - BED format
- `CNVkit/results/revised/*.vcf` - VCF format
- `CNVkit/results/revised/all_samples_heatmap.pdf` - Multi-sample heatmap

---

### snp_pileup.sh
**Purpose:** Generate SNP pileups for FACETS purity/ploidy analysis

**Resources:**
- Time: 72 hours
- CPUs: 4
- Memory: 64GB
- Array: 0-3 (processes 4 pairs)

**Processing:**
1. Downloads dbSNP common variants (first run only)
2. Filters and sorts VCF file
3. Generates SNP pileup for each tumor-normal pair
4. Sorts output for FACETS compatibility

**Outputs:**
- `FACETS/pileups/{patient_id}.sorted.pileup` - Ready for R analysis
- `FACETS/VCF/out_sorted.vcf` - Processed SNP database

**Next Step:** Run FACETS R script locally with pileup files

## ⚙️ Key Differences from Standard Linux

### 1. Shebang Line
All scripts MUST use:
```bash
#!/bin/bash -l
```
The `-l` flag loads the login environment, which initializes the module system.

### 2. SLURM Headers
Every script requires:
```bash
#SBATCH --partition=cpu
#SBATCH --account=p201093
#SBATCH --qos=default
```

### 3. Module Loading
Modules only work on compute nodes, not login nodes:
```bash
module purge
module load samtools/1.20-foss-2023a
module load picard/3.2
module load Java/21.0.5
```

### 4. Picard Syntax
Picard is run via Java:
```bash
# Correct
java -jar $PICARD_JAR CollectInsertSizeMetrics ...

# Incorrect (will fail)
picard CollectInsertSizeMetrics ...
```

### 5. Job Submission
Use `sbatch` instead of `nohup`:
```bash
# Meluxina
sbatch script.sh

# Old way (doesn't work on HPC)
nohup ./script.sh > run.log 2>&1 &
```

## 📊 Expected Runtimes

For 4 sample pairs (8 BAM files):

| Step | Time per Pair | Total Time |
|------|--------------|------------|
| BAM QC | 1-2 hours | 6-12 hours |
| BAM Processing | 6-12 hours | 24-48 hours |
| Coverage | 3-6 hours | 12-24 hours |
| CNVkit | 6-12 hours | 24-48 hours |
| SNP Pileup | 12-24 hours | 48-72 hours |

**Total pipeline:** ~5-8 days for 4 pairs

## 💾 Resource Usage

**Typical per job:**
- BAM QC: ~32GB RAM, 8 CPUs
- BAM Processing: ~64GB RAM, 16 CPUs per pair
- Coverage: ~32GB RAM, 8 CPUs
- CNVkit: ~128GB RAM, 16 CPUs
- SNP Pileup: ~64GB RAM, 4 CPUs per pair

**Disk Space (per 4 pairs):**
- Original BAMs: ~100-200 GB
- Processed BAMs: ~100-200 GB
- QC outputs: ~5-10 GB
- CNVkit results: ~10-20 GB
- FACETS pileups: ~20-40 GB
- **Total: ~250-500 GB**

## 🔍 Troubleshooting

### Module Command Not Found
**On login nodes:** This is normal - modules only work on compute nodes

**In SLURM jobs:** Check shebang line is `#!/bin/bash -l`

### Job Fails Immediately
```bash
# Check error log
cat logs/script_name_JOBID.err

# Check job details
sacct -j JOBID --format=JobID,JobName,State,ExitCode,Elapsed
```

### Out of Memory
Increase `#SBATCH --mem=` value in script header

### Job Timeout
Increase `#SBATCH --time=` value in script header

### Picard Command Not Found
Make sure you're using:
```bash
java -jar $PICARD_JAR CommandName ...
```
Not:
```bash
picard CommandName ...
```

## 📧 Support

**Meluxina Support:** servicedesk@lxp.lu  
**Documentation:** https://docs.lxp.lu/

## 🔗 Useful Commands

```bash
# Check job status
squeue -u u103499

# Check job history
sacct -u u103499 | tail -20

# Job details
scontrol show job JOBID

# Cancel job
scancel JOBID

# Cancel all jobs
scancel -u u103499

# Check quota
myquota

# Interactive session (for testing)
srun --partition=cpu --time=01:00:00 --account=p201093 --qos=default --pty /bin/bash -l

# Watch job queue
watch -n 10 'squeue -u u103499'
```

## 📝 Citation

If you use this pipeline, please cite:

- **Meluxina:** LuxProvide's Meluxina Supercomputer
- **SAMtools:** Li H, et al. (2009) Bioinformatics
- **Picard:** Broad Institute
- **CNVkit:** Talevich E, et al. (2016) PLOS Computational Biology
- **FACETS:** Shen R & Seshan VE (2016) Nucleic Acids Research

## 📄 License

This pipeline is based on original scripts by S. Owens (2025) and adapted for Meluxina HPC.

## ✨ Acknowledgments

- Meluxina HPC support team for module system guidance
- Original pipeline development: S. Owens
- Meluxina adaptation: December 2025

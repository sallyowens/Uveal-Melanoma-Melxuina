# Downloading Analysis Results from MeluXina

This guide explains how to download your WGS analysis results from MeluXina to your local computer.

## Prerequisites

- SSH access to MeluXina configured
- External hard drive mounted (for large files)
- `rsync` installed (comes with macOS/Linux)

---

## Step 1: Verify SSH Configuration

Check your SSH config file has MeluXina configured:
```bash
cat ~/.ssh/config
```

You should see something like:
```
Host meluxina
Hostname login.lxp.lu
User u103499
Port 8822
IdentityFile /Users/yourusername/.ssh/id_ed25519_mlux
IdentitiesOnly yes
ForwardAgent no
```

---

## Step 2: Create Local Directory Structure

On your **local computer**, create a directory for results:
```bash
# For Mac/Linux
mkdir -p /Volumes/External_Drive/UVM_WGS_Analysis/grp1_results

# For Windows (Git Bash)
mkdir -p /e/UVM_WGS_Analysis/grp1_results
```

---

## Step 3: Download Results

Run these commands **on your local computer** (not logged into MeluXina):

### Download BAM QC Results (~2 MB)
```bash
rsync -avz --progress meluxina:Project/bam_qc /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

### Download CNVkit Results (~5-10 GB)
```bash
rsync -avz --progress meluxina:Project/cnvkit /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

### Download FACETS Results (~500 MB - 2 GB)
```bash
rsync -avz --progress meluxina:Project/FACETS /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

### Download Coverage Results (Optional - Large! ~50-100 GB)
```bash
rsync -avz --progress meluxina:Project/coverage /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

**Note**: You'll be prompted for your SSH key passphrase for each download.

---

## Step 4: Verify Downloads

Check that files downloaded correctly:
```bash
# Navigate to results directory
cd /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/

# List downloaded directories
ls -lh

# Check key results files
ls -lh FACETS/reports/
ls -lh cnvkit/revised/
ls -lh bam_qc/multiqc_report/
```

---

## Expected Download Sizes

| Directory | Approximate Size | Contents |
|-----------|-----------------|----------|
| `bam_qc/` | ~2 MB | QC metrics, MultiQC report |
| `cnvkit/` | ~5-10 GB | CNV calls, plots (PDFs) |
| `FACETS/` | ~500 MB - 2 GB | Purity/ploidy, segments, plots |
| `coverage/` | ~50-100 GB | Depth files (optional) |

**Total (without coverage)**: ~5-12 GB per batch  
**Total (with coverage)**: ~55-112 GB per batch

---

## Troubleshooting

### Permission Denied Error
```bash
# Make sure your SSH key has correct permissions
chmod 600 ~/.ssh/id_ed25519_mlux
```

### Connection Timeout
```bash
# Check you're not connected to VPN or firewall blocking port 8822
# Try alternative command:
rsync -avz --progress -e "ssh -p 8822" u103499@login.lxp.lu:Project/bam_qc /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

### Slow Transfer Speed
```bash
# Use compression (already included with -z flag)
# Or limit bandwidth if on shared network:
rsync -avz --progress --bwlimit=10000 meluxina:Project/bam_qc /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
# (limits to 10 MB/s)
```

### Resume Interrupted Transfer
```bash
# rsync automatically resumes - just run the same command again
rsync -avz --progress meluxina:Project/cnvkit /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

---

## Alternative: GUI File Transfer (FileZilla)

If you prefer a graphical interface:

### 1. Download FileZilla
https://filezilla-project.org/download.php?type=client

### 2. Configure Connection
- **Protocol**: SFTP
- **Host**: login.lxp.lu
- **Port**: 8822
- **Logon Type**: Key file
- **User**: u103499
- **Key file**: Browse to `~/.ssh/id_ed25519_mlux`

### 3. Navigate and Download
- **Left pane**: Your computer (`/Volumes/Expansion/UVM_WGS_Analysis/grp1_results/`)
- **Right pane**: MeluXina (`/home/users/u103499/Project/`)
- Right-click folders → Download

---

## After Downloading

### Verify Results
1. Open MultiQC report: `bam_qc/multiqc_report/bam_qc_report.html`
2. Check FACETS summary: `FACETS/reports/facets_summary_all_samples.csv`
3. View CNV plots: `cnvkit/revised/*.pdf`
4. View FACETS plots: `FACETS/plots/*.png`

### Clean Up MeluXina
Once you've verified downloads, free space on MeluXina for the next batch:
```bash
# On MeluXina (after confirming downloads)
~/Project/scripts/reset_for_next_batch.sh
```

See [Processing Next Batch](processing_next_batch.md) for details.

---

## Key Result Files

### FACETS Summary
```
FACETS/reports/facets_summary_all_samples.csv
```
Contains tumor purity, ploidy, and segment counts for all samples.

### CNVkit Segments
```
cnvkit/revised/T_[SAMPLE]_revised.call.cns
```
Copy number segments with log2 ratios and integer copy numbers.

### FACETS Segments
```
FACETS/results/segments_[SAMPLE].csv
```
Allele-specific copy number segments with purity/ploidy information.

### Visualizations
```
cnvkit/revised/*.scatter.pdf         # Genome-wide CNV plots
cnvkit/revised/*.diagram.pdf         # Chromosome diagrams
FACETS/plots/facets_*.png           # FACETS genome plots
FACETS/plots/spider_*.png           # Purity-ploidy solutions
```

---

## Batch Management

For multiple batches, use this naming convention:
```
UVM_WGS_Analysis/
├── grp1_results/
│   ├── bam_qc/
│   ├── cnvkit/
│   └── FACETS/
├── grp2_results/
│   ├── bam_qc/
│   ├── cnvkit/
│   └── FACETS/
└── grp3_results/
    ├── bam_qc/
    ├── cnvkit/
    └── FACETS/
```

---

## Quick Reference
```bash
# Download essential results only (~5-12 GB)
rsync -avz --progress meluxina:Project/{bam_qc,cnvkit,FACETS} /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/

# Download everything including coverage (~55-112 GB)
rsync -avz --progress meluxina:Project/{bam_qc,cnvkit,FACETS,coverage} /Volumes/Expansion/UVM_WGS_Analysis/grp1_results/
```

---

## Notes

- **Always verify downloads** before deleting files on MeluXina
- **Coverage files are optional** - they're very large and mainly for detailed QC
- **Compress old results** on your external drive to save space
- **Keep a log** of which batches you've processed and downloaded

---

**Last updated**: December 2025

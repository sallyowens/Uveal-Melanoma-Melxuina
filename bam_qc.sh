#!/bin/bash
#SBATCH --job-name=bam_qc
#SBATCH --account=p201093
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --output=/home/users/u103499/Project/logs/bam_qc_%j.out
#SBATCH --error=/home/users/u103499/Project/logs/bam_qc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G

set -e

# === PATH CONFIGURATION ===
PROJECT_DIR=${HOME}/Project
BAM_DIR=${PROJECT_DIR}/grp1
QC_DIR=${PROJECT_DIR}/bam_qc
REFERENCE=${PROJECT_DIR}/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Initialize conda
source ${PROJECT_DIR}/software/miniconda3/etc/profile.d/conda.sh
conda activate wgs_analysis

# Create output directories
mkdir -p ${QC_DIR}/{logs,stats,multiqc_report}

echo "=== BAM Quality Check ==="
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${QC_DIR}"
date
echo ""

cd ${BAM_DIR}

# Function to check BAM file
check_bam_file() {
    local bamfile=$1
    local sample=$(basename ${bamfile} .final.bam)
    
    echo "=========================================="
    echo "Checking: ${sample}"
    echo "=========================================="
    
    # 1. Check file integrity
    echo "1. Checking BAM integrity..."
    if samtools quickcheck ${bamfile}; then
        echo "✓ BAM file is valid"
    else
        echo "✗ ERROR: BAM file is corrupted!"
        return 1
    fi
    
    # 2. Check if indexed
    echo "2. Checking BAM index..."
    if [ -f "${bamfile}.bai" ]; then
        echo "✓ Index file exists"
    else
        echo "⚠ Index missing - creating..."
        samtools index ${bamfile}
    fi
    
    # 3. Generate basic stats
    echo "3. Generating flagstat..."
    samtools flagstat ${bamfile} > ${QC_DIR}/stats/${sample}_flagstat.txt
    
    echo "4. Generating samtools stats..."
    samtools stats ${bamfile} > ${QC_DIR}/stats/${sample}_stats.txt
    
    echo "5. Generating idxstats..."
    samtools idxstats ${bamfile} > ${QC_DIR}/stats/${sample}_idxstats.txt
    
    # 4. Check read groups
    echo "6. Checking read groups..."
    samtools view -H ${bamfile} | grep "^@RG" > ${QC_DIR}/stats/${sample}_readgroups.txt 2>&1 || \
        echo "⚠ WARNING: No read groups found" > ${QC_DIR}/stats/${sample}_readgroups.txt
    
    # 5. Calculate coverage statistics
    echo "7. Calculating coverage (per chromosome)..."
    samtools coverage ${bamfile} > ${QC_DIR}/stats/${sample}_coverage.txt
    
    # 6. Alignment metrics with Picard
    echo "8. Running Picard CollectAlignmentSummaryMetrics..."
    picard CollectAlignmentSummaryMetrics \
        I=${bamfile} \
        O=${QC_DIR}/stats/${sample}_alignment_metrics.txt \
        R=${REFERENCE} \
        VALIDATION_STRINGENCY=LENIENT
    
    # 7. Check for duplicates
    echo "9. Checking duplicate rate..."
    picard MarkDuplicates \
        I=${bamfile} \
        O=/dev/null \
        M=${QC_DIR}/stats/${sample}_duplication_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        ASSUME_SORTED=true
    
    echo "✓ QC completed for ${sample}"
    echo ""
}

# Process all BAM files
for bamfile in *.final.bam; do
    if [ -f "$bamfile" ]; then
        check_bam_file "$bamfile" || echo "⚠ Failed QC for $bamfile"
    fi
done

# Generate summary report
echo "Generating summary report..."

{
    echo "=================================="
    echo "BAM Quality Check Summary Report"
    echo "=================================="
    echo "Date: $(date)"
    echo "BAM Location: ${BAM_DIR}"
    echo ""
    echo "Files analyzed:"
    ls -lh ${BAM_DIR}/*.final.bam
    echo ""
    echo "=== SAMPLE OVERVIEW ==="
    echo ""
    for bamfile in ${BAM_DIR}/*.final.bam; do
        sample=$(basename ${bamfile} .final.bam)
        size=$(du -h ${bamfile} | cut -f1)
        echo "${sample}: ${size}"
    done
    echo ""
    echo "=== ALIGNMENT STATISTICS ==="
    echo ""
    
    for flagstat in ${QC_DIR}/stats/*_flagstat.txt; do
        sample=$(basename ${flagstat} _flagstat.txt)
        echo "--- ${sample} ---"
        cat ${flagstat}
        echo ""
    done
    
    echo "=== COVERAGE SUMMARY ==="
    echo ""
    
    for coverage in ${QC_DIR}/stats/*_coverage.txt; do
        sample=$(basename ${coverage} _coverage.txt)
        echo "--- ${sample} ---"
        awk 'NR==1 {print; next}
             NR>1 {
                 total_bases+=$3; 
                 covered_bases+=$4; 
                 total_depth+=($4*$7)
             } 
             END {
                 if(total_bases>0) {
                     printf "Overall genome coverage: %.2f%%\n", (covered_bases/total_bases)*100;
                     printf "Mean depth of covered regions: %.2fx\n", total_depth/covered_bases;
                     printf "Total bases in reference: %d\n", total_bases;
                     printf "Total covered bases: %d\n\n", covered_bases
                 }
             }' ${coverage}
    done
    
    echo "=== DUPLICATION RATES ==="
    echo ""
    
    for dup in ${QC_DIR}/stats/*_duplication_metrics.txt; do
        sample=$(basename ${dup} _duplication_metrics.txt)
        echo "--- ${sample} ---"
        grep -A 1 "^LIBRARY" ${dup} | tail -1
        echo ""
    done
    
} > ${QC_DIR}/BAM_QC_Summary_Report.txt

# Run MultiQC
echo "Running MultiQC to aggregate results..."
multiqc ${QC_DIR}/stats -o ${QC_DIR}/multiqc_report -n bam_qc_report --force

echo ""
echo "=== BAM QC Complete! ==="
echo "Results in: ${QC_DIR}"
echo "Summary report: ${QC_DIR}/BAM_QC_Summary_Report.txt"
echo "MultiQC HTML report: ${QC_DIR}/multiqc_report/bam_qc_report.html"
echo ""
date

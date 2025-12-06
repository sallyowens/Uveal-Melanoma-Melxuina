#!/bin/bash
#SBATCH --job-name=bam_qc
#SBATCH --output=logs/bam_qc_%j.out
#SBATCH --error=logs/bam_qc_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu
#SBATCH --qos=default

# BAM Quality Check Script for Meluxina
# Checks BAM files before CNV and variant calling analysis

set -e  # Exit on error

echo "Starting BAM quality check..."
date

# Load required modules (adjust versions as available on Meluxina)
module purge
module load SAMtools
module load Picard
module load Python  # for multiqc

# Set directories
BAM_DIR=/home/users/u103499/Project/grp1
QC_DIR=/home/users/u103499/Project/bam_qc
STATS_DIR=${QC_DIR}/stats
METRICS_DIR=${QC_DIR}/metrics
PLOTS_DIR=${QC_DIR}/plots
REFERENCE=/home/users/u103499/Project/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Create output directories
mkdir -p ${QC_DIR}
mkdir -p ${STATS_DIR}
mkdir -p ${METRICS_DIR}
mkdir -p ${PLOTS_DIR}
mkdir -p logs

echo "Output directory: ${QC_DIR}"
echo "Processing BAM files from: ${BAM_DIR}"

# Check if BAM directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "ERROR: BAM directory does not exist: $BAM_DIR"
    exit 1
fi

# Count BAM files
n_bams=$(ls ${BAM_DIR}/*.bam 2>/dev/null | wc -l)
echo "Found ${n_bams} BAM files to process"

if [ $n_bams -eq 0 ]; then
    echo "ERROR: No BAM files found in ${BAM_DIR}"
    exit 1
fi

# Process each BAM file
for bamfile in ${BAM_DIR}/*.bam; do
    if [ -f "$bamfile" ]; then
        sample=$(basename "$bamfile" .bam)
        echo ""
        echo "========================================="
        echo "Processing: $sample"
        echo "========================================="
        
        # Check if BAM is indexed, if not create index
        if [ ! -f "${bamfile}.bai" ]; then
            echo "Creating BAM index..."
            samtools index -@ 8 "$bamfile"
        fi
        
        # 1. Basic BAM statistics
        echo "Generating basic statistics..."
        samtools flagstat -@ 8 "$bamfile" > "${STATS_DIR}/${sample}_flagstat.txt"
        samtools stats -@ 8 "$bamfile" > "${STATS_DIR}/${sample}_stats.txt"
        
        # 2. Coverage statistics per chromosome
        echo "Calculating coverage statistics..."
        samtools coverage "$bamfile" > "${STATS_DIR}/${sample}_coverage.txt"
        
        # 3. Insert size metrics (for paired-end data)
        echo "Collecting insert size metrics..."
        picard CollectInsertSizeMetrics \
            I="$bamfile" \
            O="${METRICS_DIR}/${sample}_insert_size_metrics.txt" \
            H="${PLOTS_DIR}/${sample}_insert_size_histogram.pdf" \
            VALIDATION_STRINGENCY=LENIENT
        
        # 4. Alignment summary metrics
        echo "Collecting alignment summary metrics..."
        picard CollectAlignmentSummaryMetrics \
            I="$bamfile" \
            O="${METRICS_DIR}/${sample}_alignment_metrics.txt" \
            R=/project/home/p201093/alignment/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
            VALIDATION_STRINGENCY=LENIENT
        
        # 5. GC bias metrics
        echo "Collecting GC bias metrics..."
        picard CollectGcBiasMetrics \
            I="$bamfile" \
            O="${METRICS_DIR}/${sample}_gc_bias_metrics.txt" \
            CHART="${PLOTS_DIR}/${sample}_gc_bias_metrics.pdf" \
            S="${METRICS_DIR}/${sample}_gc_bias_summary.txt" \
            R=/project/home/p201093/alignment/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
            VALIDATION_STRINGENCY=LENIENT
        
        # 6. Extract key metrics for summary
        echo "Extracting key QC metrics..."
        
        # Total reads
        total_reads=$(grep "^SN" "${STATS_DIR}/${sample}_stats.txt" | grep "raw total sequences:" | awk '{print $NF}')
        
        # Mapped reads
        mapped_reads=$(grep "^SN" "${STATS_DIR}/${sample}_stats.txt" | grep "reads mapped:" | awk '{print $NF}')
        
        # Mapping percentage
        mapping_pct=$(grep "mapped (" "${STATS_DIR}/${sample}_flagstat.txt" | head -1 | sed 's/.*(\(.*\):.*/\1/')
        
        # Duplicate rate
        duplicates=$(grep "duplicates" "${STATS_DIR}/${sample}_flagstat.txt" | awk '{print $1}')
        dup_rate=$(awk "BEGIN {printf \"%.2f\", ($duplicates/$total_reads)*100}")
        
        # Mean coverage
        mean_cov=$(awk 'NR>1 {sum+=$7*$3; bases+=$3} END {if(bases>0) print sum/bases; else print 0}' "${STATS_DIR}/${sample}_coverage.txt")
        
        # Properly paired
        properly_paired_pct=$(grep "properly paired (" "${STATS_DIR}/${sample}_flagstat.txt" | sed 's/.*(\(.*\):.*/\1/')
        
        # Mean insert size
        mean_insert=$(grep "^MEDIAN_INSERT_SIZE" "${METRICS_DIR}/${sample}_insert_size_metrics.txt" | tail -1 | awk '{print $1}')
        
        echo "Sample: $sample" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "  Total reads: $total_reads" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "  Mapped reads: $mapped_reads (${mapping_pct})" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "  Duplicate rate: ${dup_rate}%" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "  Mean coverage: ${mean_cov}x" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "  Properly paired: ${properly_paired_pct}" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "  Mean insert size: ${mean_insert}" >> "${QC_DIR}/qc_summary_temp.txt"
        echo "" >> "${QC_DIR}/qc_summary_temp.txt"
        
        # Create CSV entry for R analysis
        echo "${sample},${total_reads},${mapped_reads},${mapping_pct},${dup_rate},${mean_cov},${properly_paired_pct},${mean_insert}" >> "${QC_DIR}/qc_metrics.csv"
        
        echo "Completed QC for: $sample"
    fi
done

# Create CSV header
sed -i '1i\sample,total_reads,mapped_reads,mapping_pct,duplicate_rate,mean_coverage,properly_paired_pct,mean_insert_size' "${QC_DIR}/qc_metrics.csv"

# Generate combined summary report
{
    echo "=========================================="
    echo "BAM Quality Control Summary Report"
    echo "=========================================="
    echo "Analysis completed: $(date)"
    echo "Total BAM files analyzed: ${n_bams}"
    echo ""
    echo "Output directories:"
    echo "  - Statistics: ${STATS_DIR}"
    echo "  - Metrics: ${METRICS_DIR}"
    echo "  - Plots: ${PLOTS_DIR}"
    echo ""
    echo "=========================================="
    echo "Individual Sample Metrics:"
    echo "=========================================="
    cat "${QC_DIR}/qc_summary_temp.txt"
    echo ""
    echo "=========================================="
    echo "Quality Flags to Review:"
    echo "=========================================="
    
    # Check for potential issues
    awk -F',' 'NR>1 {
        if ($4 > 30) print "WARNING: High duplicate rate for " $1 ": " $4 "%"
        if ($3 < 95) print "WARNING: Low mapping rate for " $1 ": " $3
        if ($6 < 70) print "WARNING: Low properly paired rate for " $1 ": " $6
        if ($5 < 20) print "WARNING: Low coverage for " $1 ": " $5 "x"
    }' "${QC_DIR}/qc_metrics.csv"
    
    echo ""
    echo "CSV summary saved: ${QC_DIR}/qc_metrics.csv"
    echo "Use this file for R visualization and detailed analysis"
    
} > "${QC_DIR}/BAM_QC_Report.txt"

# Clean up temp file
rm -f "${QC_DIR}/qc_summary_temp.txt"

echo ""
echo "========================================="
echo "BAM Quality Check Completed Successfully!"
echo "========================================="
echo "Summary report: ${QC_DIR}/BAM_QC_Report.txt"
echo "CSV metrics: ${QC_DIR}/qc_metrics.csv"
echo ""

date

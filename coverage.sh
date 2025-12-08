#!/bin/bash -l
#SBATCH --job-name=coverage
#SBATCH --output=logs/coverage_%j.out
#SBATCH --error=logs/coverage_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu
#SBATCH --account=p201093
#SBATCH --qos=default

# Coverage Assessment for WGS - Meluxina HPC
# Author: S. Owens (adapted for Meluxina, December 2025)
# Description: Calculate genome-wide coverage statistics

set -e

echo "Starting coverage analysis..."
date

# Load required modules
module purge
module load samtools/1.20-foss-2023a

# Set paths
ALIGNED=/home/users/u103499/Project
FILTERED_BAM=${ALIGNED}/processed_bams/final_bams
COVERAGE_DIR=/home/users/u103499/Project/coverage
DEPTH_DIR=${COVERAGE_DIR}/depth_files
SUMMARY_DIR=${COVERAGE_DIR}/summary

# Create directories
mkdir -p ${COVERAGE_DIR}/coverage_stats
mkdir -p ${COVERAGE_DIR}/plots_data
mkdir -p ${COVERAGE_DIR}/logs
mkdir -p ${DEPTH_DIR}
mkdir -p ${SUMMARY_DIR}
mkdir -p logs

echo "Coverage analysis directories created"
echo "Processing BAM files from: ${FILTERED_BAM}"

# Count BAM files
n_bams=$(ls ${FILTERED_BAM}/*_final.bam 2>/dev/null | wc -l)
echo "Found ${n_bams} BAM files"

if [ $n_bams -eq 0 ]; then
    echo "ERROR: No BAM files found"
    exit 1
fi

# Process each BAM file
for bamfile in ${FILTERED_BAM}/*_final.bam; do
    if [ -f "$bamfile" ]; then
        sample=$(basename "$bamfile" | sed 's/_final.bam//')
        echo ""
        echo "========================================="
        echo "Processing sample: $sample"
        echo "========================================="
        
        # Index if needed
        if [ ! -f "${bamfile}.bai" ]; then
            echo "Indexing $bamfile"
            samtools index -@ ${SLURM_CPUS_PER_TASK} "$bamfile"
        fi
        
        # Generate coverage statistics per chromosome
        echo "Generating coverage statistics..."
        samtools coverage "$bamfile" > "${COVERAGE_DIR}/coverage_stats/${sample}_coverage.txt"
        
        # Generate coverage histogram (faster than full depth)
        echo "Generating coverage histogram..."
        samtools depth -@ ${SLURM_CPUS_PER_TASK} "$bamfile" | \
            awk '{print $3}' | sort -n | uniq -c > "${COVERAGE_DIR}/plots_data/${sample}_depth_histogram.txt"
        
        # Generate coverage summary
        echo "Generating coverage summary..."
        {
            echo "Sample: $sample"
            echo "Date: $(date)"
            echo "BAM file: $bamfile"
            echo ""
            echo "=== Coverage Summary ==="
            samtools coverage "$bamfile" | head -1
            samtools coverage "$bamfile" | tail -n +2 | awk '
            BEGIN {
                total_bases = 0
                covered_bases = 0
                total_depth = 0
            }
            {
                total_bases += $3
                covered_bases += $4
                total_depth += ($4 * $7)
            }
            END {
                if (total_bases > 0) {
                    overall_coverage = (covered_bases / total_bases) * 100
                    mean_depth = total_depth / covered_bases
                    print "Overall genome coverage: " overall_coverage "%"
                    print "Mean depth of covered regions: " mean_depth "x"
                    print "Total bases in reference: " total_bases
                    print "Total covered bases: " covered_bases
                }
            }'
            echo ""
            echo "=== Per-chromosome coverage ==="
            samtools coverage "$bamfile"
        } > "${COVERAGE_DIR}/coverage_stats/${sample}_summary.txt"
        
        echo "Completed: $sample"
    fi
done

# Generate combined summary
echo ""
echo "Generating combined summary..."
{
    echo "WGS Coverage Analysis Summary"
    echo "Generated on: $(date)"
    echo "=============================="
    echo ""
    
    for summary_file in ${COVERAGE_DIR}/coverage_stats/*_summary.txt; do
        if [ -f "$summary_file" ]; then
            cat "$summary_file"
            echo ""
            echo "---"
            echo ""
        fi
    done
} > "${COVERAGE_DIR}/combined_coverage_summary.txt"

# Create R-ready data files
echo "Preparing R-ready data..."

# Combine all coverage data
{
    echo -e "sample\tchromosome\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
    for coverage_file in ${COVERAGE_DIR}/coverage_stats/*_coverage.txt; do
        if [ -f "$coverage_file" ]; then
            sample=$(basename "$coverage_file" | sed 's/_coverage.txt//')
            tail -n +2 "$coverage_file" | awk -v sample="$sample" 'BEGIN{OFS="\t"} {print sample, $0}'
        fi
    done
} > "${COVERAGE_DIR}/plots_data/all_samples_coverage.txt"

# Create sample metadata
echo "Creating sample metadata..."
{
    echo -e "sample\tsample_type\tbam_file"
    for bamfile in ${FILTERED_BAM}/*_final.bam; do
        if [ -f "$bamfile" ]; then
            sample=$(basename "$bamfile" | sed 's/_final.bam//')
            
            # Classify as tumor or normal based on filename
            if [[ "$sample" =~ _tumor$ ]]; then
                sample_type="tumor"
            elif [[ "$sample" =~ _normal$ ]]; then
                sample_type="normal"
            else
                sample_type="unknown"
            fi
            echo -e "${sample}\t${sample_type}\t${bamfile}"
        fi
    done
} > "${COVERAGE_DIR}/plots_data/sample_metadata.txt"

echo ""
echo "========================================="
echo "Coverage analysis completed!"
echo "========================================="
echo "Results in: $COVERAGE_DIR"
echo "  - Individual stats: ${COVERAGE_DIR}/coverage_stats/"
echo "  - Combined summary: ${COVERAGE_DIR}/combined_coverage_summary.txt"
echo "  - R-ready data: ${COVERAGE_DIR}/plots_data/"
echo ""
date

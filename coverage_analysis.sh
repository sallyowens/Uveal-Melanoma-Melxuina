#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --account=p201093
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --output=/home/users/u103499/Project/logs/coverage_%j.out
#SBATCH --error=/home/users/u103499/Project/logs/coverage_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=64G

set -e

# === PATH CONFIGURATION ===
PROJECT_DIR=${HOME}/Project
BAM_DIR=${PROJECT_DIR}/grp1
COVERAGE_DIR=${PROJECT_DIR}/coverage

# Initialize conda
source ${PROJECT_DIR}/software/miniconda3/etc/profile.d/conda.sh
conda activate wgs_analysis

# Create directories
mkdir -p ${COVERAGE_DIR}/{depth_files,coverage_stats,plots_data}

echo "=== Coverage Analysis ==="
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${COVERAGE_DIR}"
date
echo ""

cd ${BAM_DIR}

# Generate coverage for each BAM
for bamfile in *.final.bam; do
    if [ -f "$bamfile" ]; then
        sample=$(basename ${bamfile} .final.bam)
        echo "=========================================="
        echo "Processing: ${sample}"
        echo "=========================================="
        
        # Ensure indexed
        if [ ! -f "${bamfile}.bai" ]; then
            echo "Indexing BAM..."
            samtools index ${bamfile}
        fi
        
        # Coverage statistics per chromosome
        echo "Generating coverage statistics..."
        samtools coverage ${bamfile} > ${COVERAGE_DIR}/coverage_stats/${sample}_coverage.txt
        
        # Depth file (compressed to save space)
        echo "Generating depth file (this will take a while)..."
        samtools depth ${bamfile} | gzip > ${COVERAGE_DIR}/depth_files/${sample}_depth.txt.gz
        
        # Histogram of depth distribution
        echo "Generating depth histogram..."
        zcat ${COVERAGE_DIR}/depth_files/${sample}_depth.txt.gz | \
            awk '{print $3}' | \
            sort -n | uniq -c > ${COVERAGE_DIR}/plots_data/${sample}_depth_histogram.txt
        
        echo "✓ Completed: ${sample}"
        echo ""
    fi
done

# Create combined coverage file for R
echo "Creating combined coverage summary..."
{
    echo -e "sample\tchromosome\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
    for covfile in ${COVERAGE_DIR}/coverage_stats/*_coverage.txt; do
        sample=$(basename ${covfile} _coverage.txt)
        tail -n +2 ${covfile} | awk -v sample="${sample}" 'BEGIN{OFS="\t"} {print sample, $0}'
    done
} > ${COVERAGE_DIR}/plots_data/all_samples_coverage.txt

# Create sample metadata
echo "Creating sample metadata..."
{
    echo -e "sample\tsample_type\tbam_file"
    for bamfile in ${BAM_DIR}/*.final.bam; do
        sample=$(basename ${bamfile} .final.bam)
        if [[ ${sample} =~ ^T_ ]]; then
            sample_type="tumor"
        elif [[ ${sample} =~ ^B_ ]]; then
            sample_type="normal"
        else
            sample_type="unknown"
        fi
        echo -e "${sample}\t${sample_type}\t${bamfile}"
    done
} > ${COVERAGE_DIR}/plots_data/sample_metadata.txt

# Generate summary
{
    echo "Coverage Analysis Summary"
    echo "========================="
    echo "Date: $(date)"
    echo ""
    echo "Samples processed:"
    cat ${COVERAGE_DIR}/plots_data/sample_metadata.txt
    echo ""
    echo "Output files:"
    echo "- Coverage stats: ${COVERAGE_DIR}/coverage_stats/"
    echo "- Depth files: ${COVERAGE_DIR}/depth_files/"
    echo "- Plot data: ${COVERAGE_DIR}/plots_data/"
} > ${COVERAGE_DIR}/coverage_summary.txt

echo ""
echo "=== Coverage Analysis Complete! ==="
echo "Results in: ${COVERAGE_DIR}"
echo "Summary: ${COVERAGE_DIR}/coverage_summary.txt"
echo ""
date

#!/bin/bash -l
#SBATCH --job-name=cnvkit
#SBATCH --output=logs/cnvkit_%j.out
#SBATCH --error=logs/cnvkit_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=cpu
#SBATCH --account=p201093
#SBATCH --qos=default

# CNVkit Analysis for WGS on Meluxina HPC
# Author: S. Owens (adapted for Meluxina, December 2025)
# Description: Copy number variation detection and calling

set -e

echo "Starting CNVkit analysis..."
date

# Load required modules
module purge
module load Python/3.12.3-GCCcore-13.3.0
module load samtools/1.20-foss-2023a
module load R/4.4.2-gfbf-2024a

# Set directories
REFERENCE=/home/users/u103499/Project/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
BEDFILES=/home/users/u103499/Project/CNVkit/bedfiles
BAMFILES=/home/users/u103499/Project/processed_bams/final_bams
CNVKIT=/home/users/u103499/Project/CNVkit/results

# Create directories
mkdir -p ${BEDFILES}
mkdir -p ${CNVKIT}
mkdir -p ${CNVKIT}/revised
mkdir -p logs

echo "CNVkit directories created"
echo "Reference: ${REFERENCE}"
echo "BAM files: ${BAMFILES}"
echo "Output: ${CNVKIT}"

# Create access bedfile (only needs to be done once)
if [ ! -f "${BEDFILES}/access.bed" ]; then
    echo "Creating access BED file..."
    cnvkit.py access ${REFERENCE} -o ${BEDFILES}/access.bed
    echo "Access BED file created"
fi

# Separate tumor and normal BAM files
cd ${BAMFILES}
TUMOR_BAMS=$(ls *_tumor_final.bam 2>/dev/null | tr '\n' ' ')
NORMAL_BAMS=$(ls *_normal_final.bam 2>/dev/null | tr '\n' ' ')

echo "Tumor BAMs found:"
echo "$TUMOR_BAMS"
echo ""
echo "Normal BAMs found:"
echo "$NORMAL_BAMS"
echo ""

if [ -z "$TUMOR_BAMS" ]; then
    echo "ERROR: No tumor BAM files found"
    exit 1
fi

if [ -z "$NORMAL_BAMS" ]; then
    echo "ERROR: No normal BAM files found"
    exit 1
fi

#### PHASE ONE: Initial CNV calling ####
echo "========================================="
echo "PHASE ONE: Initial CNV calling"
echo "========================================="
date

cnvkit.py batch ${BAMFILES}/*_tumor_final.bam \
    --normal ${BAMFILES}/*_normal_final.bam \
    --method wgs \
    --fasta ${REFERENCE} \
    --access ${BEDFILES}/access.bed \
    --output-dir ${CNVKIT} \
    --scatter \
    --diagram \
    -p ${SLURM_CPUS_PER_TASK}

echo "Phase one completed"
date
echo ""

#### PHASE TWO: Refine with purity normalization ####
echo "========================================="
echo "PHASE TWO: Purity normalization"
echo "========================================="
date

# Process each .cnr file
for cnr_file in ${CNVKIT}/*.cnr; do
    if [ -f "$cnr_file" ]; then
        sample=$(basename "$cnr_file" .cnr)
        echo "Processing: $sample"
        
        # Resegment and drop low coverage
        echo "  Resegmenting..."
        cnvkit.py segment "$cnr_file" \
            --drop-low-coverage \
            -o ${CNVKIT}/revised/${sample}.revised.cns \
            -p ${SLURM_CPUS_PER_TASK}
        
        # Call copy numbers
        echo "  Calling copy numbers..."
        cnvkit.py call ${CNVKIT}/revised/${sample}.revised.cns \
            -y \
            -m clonal \
            -o ${CNVKIT}/revised/${sample}.revised.call.cns
        
        # Generate scatter plot
        echo "  Generating scatter plot..."
        cnvkit.py scatter "$cnr_file" \
            -s ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.revised.call.scatter.png
        
        # Generate diagram
        echo "  Generating diagram..."
        cnvkit.py diagram "$cnr_file" \
            -s ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.revised.call.diagram.pdf
        
        # Export formats
        echo "  Exporting results..."
        
        cnvkit.py export bed ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.bed
        
        cnvkit.py export seg ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.seg
        
        cnvkit.py export vcf ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.vcf
        
        echo "  Completed: $sample"
        echo ""
    fi
done

# Generate metrics
echo "Generating summary metrics..."
cnvkit.py metrics ${CNVKIT}/*.cnr -s ${CNVKIT}/*.cns \
    -o ${CNVKIT}/cnvkit_metrics.txt

# Create heatmap
echo "Creating heatmap..."
cnvkit.py heatmap ${CNVKIT}/revised/*.revised.call.cns \
    -o ${CNVKIT}/revised/all_samples_heatmap.pdf

# Generate summary report
{
    echo "CNVkit WGS Analysis Summary"
    echo "==========================="
    echo "Analysis completed: $(date)"
    echo ""
    echo "Parameters:"
    echo "  Method: wgs"
    echo "  Reference: ${REFERENCE}"
    echo ""
    echo "Samples processed:"
    ls -1 ${CNVKIT}/revised/*.revised.call.cns | wc -l
    echo ""
    echo "Output files:"
    echo "  - Results: ${CNVKIT}/"
    echo "  - Refined: ${CNVKIT}/revised/"
    echo ""
} > ${CNVKIT}/cnvkit_summary.txt

echo "========================================="
echo "CNVkit analysis completed!"
echo "========================================="
echo "Summary: ${CNVKIT}/cnvkit_summary.txt"
date

#!/bin/bash
#SBATCH --job-name=cnvkit
#SBATCH --output=logs/cnvkit_%j.out
#SBATCH --error=logs/cnvkit_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=cpu
#SBATCH --qos=default

# CNVkit Analysis for WGS on Meluxina

set -e

echo "Starting CNVkit analysis..."
date

# Load required modules
module purge
module load Python
module load SAMtools
module load R

# Activate conda environment for CNVkit (if using conda)
# Or ensure CNVkit is in your PATH
# conda activate cnvkit_env  # Uncomment if using conda

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

# For WGS, we don't need a target BED file
# CNVkit will use bins across the genome

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

#### PHASE ONE: Initial CNV calling (no normalization) ####
echo "========================================="
echo "PHASE ONE: Initial CNV calling"
echo "========================================="
date

# For WGS, use method 'wgs' and appropriate parameters
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
        
        # Drop low coverage regions and resegment
        echo "  Resegmenting and dropping low coverage..."
        cnvkit.py segment "$cnr_file" \
            --drop-low-coverage \
            -o ${CNVKIT}/revised/${sample}.revised.cns \
            -p ${SLURM_CPUS_PER_TASK}
        
        # Call copy numbers with clonal purity normalization
        # If you have known purity values from FACETS, add --purity <value>
        echo "  Calling copy numbers with purity normalization..."
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
        
        # Export to various formats
        echo "  Exporting results..."
        
        # BED format
        cnvkit.py export bed ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.bed
        
        # SEG format (for IGV)
        cnvkit.py export seg ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.seg
        
        # VCF format
        cnvkit.py export vcf ${CNVKIT}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT}/revised/${sample}.vcf
        
        echo "  Completed: $sample"
        echo ""
    fi
done

# Generate summary metrics
echo "Generating summary metrics..."
cnvkit.py metrics ${CNVKIT}/*.cnr -s ${CNVKIT}/*.cns \
    -o ${CNVKIT}/cnvkit_metrics.txt

# Create heatmap of all samples
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
    echo "  Access file: ${BEDFILES}/access.bed"
    echo ""
    echo "Output files:"
    echo "  Initial results: ${CNVKIT}/"
    echo "  Refined results: ${CNVKIT}/revised/"
    echo ""
    echo "File types generated:"
    echo "  - .cnr: Copy number ratios"
    echo "  - .cns: Segmented copy numbers"
    echo "  - .bed: BED format export"
    echo "  - .seg: SEG format for IGV"
    echo "  - .vcf: VCF format"
    echo "  - .png: Scatter plots"
    echo "  - .pdf: Diagrams and heatmap"
    echo ""
    echo "Samples processed:"
    ls -1 ${CNVKIT}/revised/*.revised.call.cns | wc -l
    echo ""
    echo "Next steps:"
    echo "1. Review scatter plots and diagrams"
    echo "2. Load .seg files into IGV for visualization"
    echo "3. Compare with FACETS results"
    echo "4. If you have FACETS purity estimates, rerun with --purity flag"
} > ${CNVKIT}/cnvkit_summary.txt

echo "========================================="
echo "CNVkit analysis completed!"
echo "========================================="
echo "Summary: ${CNVKIT}/cnvkit_summary.txt"
echo "Results: ${CNVKIT}/revised/"
echo ""
date

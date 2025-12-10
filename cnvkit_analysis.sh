#!/bin/bash
#SBATCH --job-name=cnvkit
#SBATCH --account=p201093
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --output=/home/users/u103499/Project/logs/cnvkit_%j.out
#SBATCH --error=/home/users/u103499/Project/logs/cnvkit_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=128G

set -e

# === PATH CONFIGURATION ===
PROJECT_DIR=${HOME}/Project
BAM_DIR=${PROJECT_DIR}/grp1
REFERENCE=${PROJECT_DIR}/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
CNVKIT_DIR=${PROJECT_DIR}/cnvkit

# Initialize conda
source ${PROJECT_DIR}/software/miniconda3/etc/profile.d/conda.sh
conda activate wgs_analysis

# Create directories
mkdir -p ${CNVKIT_DIR}/{raw,revised,logs}

echo "=== CNVkit Analysis for WGS ==="
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${CNVKIT_DIR}"
date
echo ""

cd ${BAM_DIR}

# Separate tumor and normal files
TUMOR_BAMS=(T_*.final.bam)
NORMAL_BAMS=(B_*.final.bam)

echo "Found ${#TUMOR_BAMS[@]} tumor samples:"
printf '%s\n' "${TUMOR_BAMS[@]}"
echo ""
echo "Found ${#NORMAL_BAMS[@]} normal samples:"
printf '%s\n' "${NORMAL_BAMS[@]}"
echo ""

# Create access file if it doesn't exist
if [ ! -f "${CNVKIT_DIR}/access.bed" ]; then
    echo "Creating access file..."
    cnvkit.py access ${REFERENCE} -o ${CNVKIT_DIR}/access.bed
fi

# Run CNVkit batch mode for WGS
echo "Running CNVkit batch analysis..."
echo "This will take several hours for WGS data..."

cnvkit.py batch ${TUMOR_BAMS[@]} \
    --normal ${NORMAL_BAMS[@]} \
    --method wgs \
    --fasta ${REFERENCE} \
    --access ${CNVKIT_DIR}/access.bed \
    --output-dir ${CNVKIT_DIR}/raw \
    --processes 16 \
    --diagram \
    --scatter

echo ""
echo "CNVkit batch complete. Starting refinement..."
echo ""

# Phase 2: Refine segments
cd ${CNVKIT_DIR}/raw

for cnr in *.cnr; do
    if [ -f "$cnr" ]; then
        sample=$(basename ${cnr} .cnr)
        echo "=========================================="
        echo "Refining: ${sample}"
        echo "=========================================="
        
        # Drop low coverage regions
        echo "Dropping low coverage regions..."
        cnvkit.py segment ${cnr} \
            --drop-low-coverage \
            -o ${CNVKIT_DIR}/revised/${sample}.revised.cns
        
        # Call copy numbers (clonal mode, no purity correction)
        echo "Calling copy numbers..."
        cnvkit.py call ${CNVKIT_DIR}/revised/${sample}.revised.cns \
            -y -m clonal \
            -o ${CNVKIT_DIR}/revised/${sample}.revised.call.cns
        
        # Generate scatter plot
        echo "Generating scatter plot..."
        cnvkit.py scatter ${cnr} \
            -s ${CNVKIT_DIR}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT_DIR}/revised/${sample}.revised.call.scatter.pdf
        
        # Generate chromosome diagrams
        echo "Generating chromosome diagram..."
        cnvkit.py diagram ${cnr} \
            -s ${CNVKIT_DIR}/revised/${sample}.revised.call.cns \
            -o ${CNVKIT_DIR}/revised/${sample}.revised.call.diagram.pdf
        
        echo "✓ Completed: ${sample}"
        echo ""
    fi
done

# Generate summary
{
    echo "CNVkit Analysis Summary"
    echo "======================="
    echo "Date: $(date)"
    echo ""
    echo "Analysis parameters:"
    echo "- Method: WGS"
    echo "- Reference: ${REFERENCE}"
    echo "- Tumor samples: ${#TUMOR_BAMS[@]}"
    echo "- Normal samples: ${#NORMAL_BAMS[@]}"
    echo ""
    echo "Output files:"
    echo "- Raw results: ${CNVKIT_DIR}/raw/"
    echo "- Refined segments: ${CNVKIT_DIR}/revised/*.revised.call.cns"
    echo "- Plots: ${CNVKIT_DIR}/revised/*.pdf"
    echo ""
    echo "Files generated:"
    ls -lh ${CNVKIT_DIR}/revised/*.cns 2>/dev/null || echo "No files yet"
} > ${CNVKIT_DIR}/cnvkit_summary.txt

echo ""
echo "=== CNVkit Analysis Complete! ==="
echo "Results in: ${CNVKIT_DIR}"
echo "Summary: ${CNVKIT_DIR}/cnvkit_summary.txt"
echo ""
date

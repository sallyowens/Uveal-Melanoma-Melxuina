#!/bin/bash -l
#SBATCH --job-name=snp_pileup
#SBATCH --output=logs/snp_pileup_%A_%a.out
#SBATCH --error=logs/snp_pileup_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=cpu
#SBATCH --account=p201093
#SBATCH --qos=default
#SBATCH --array=0-3

# SNP Pileup for FACETS WGS Analysis on Meluxina HPC
# Author: S. Owens (adapted for Meluxina, December 2025)
# Description: Generate SNP pileups for FACETS purity/ploidy analysis

set -e

echo "Starting SNP pileup for FACETS..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
date

# Load required modules
module purge
module load samtools/1.20-foss-2023a

# Set directories
BAMFILES=/home/users/u103499/Project/processed_bams/final_bams
FACETS_DIR=/home/users/u103499/Project/FACETS
VCF_DIR=${FACETS_DIR}/VCF
PILEUP_DIR=${FACETS_DIR}/pileups
LOG_DIR=${FACETS_DIR}/logs

# Create directories
mkdir -p ${VCF_DIR}
mkdir -p ${PILEUP_DIR}
mkdir -p ${LOG_DIR}
mkdir -p logs

echo "FACETS directories created"

# Download and prepare VCF file (only array task 0)
if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    if [ ! -f "${VCF_DIR}/out_sorted.vcf" ]; then
        echo "Downloading and processing VCF file..."
        
        # Download common SNPs
        if [ ! -f "${VCF_DIR}/00-common_all.vcf.gz" ]; then
            wget -P ${VCF_DIR} https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
        fi
        
        # Process and sort VCF
        echo "Processing VCF file..."
        gunzip -c ${VCF_DIR}/00-common_all.vcf.gz | \
        awk 'BEGIN{OFS="\t"} 
             /^#/ {print; next}
             /^[0-9XY]+\t/ && ($6 == "." || $6 >= 30) {
                 if($1 ~ /^(chr)?[0-9XY]+$/) print
             }' > ${VCF_DIR}/temp_filtered.vcf
        
        # Sort
        (grep '^#' ${VCF_DIR}/temp_filtered.vcf; 
         grep -v '^#' ${VCF_DIR}/temp_filtered.vcf | sort -k1,1V -k2,2n) > ${VCF_DIR}/out_sorted.vcf
        
        rm -f ${VCF_DIR}/temp_filtered.vcf
        echo "VCF file ready"
    fi
fi

# Wait for VCF
sleep 60

# Get list of tumor files
cd ${BAMFILES}
TUMOR_FILES=($(ls *_tumor_final.bam 2>/dev/null | sort))

# Get current tumor file
CURRENT_TUMOR=${TUMOR_FILES[${SLURM_ARRAY_TASK_ID}]}

if [ -z "$CURRENT_TUMOR" ]; then
    echo "No tumor file for array task ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

echo "Processing tumor: $CURRENT_TUMOR"

# Extract patient ID
patient_id=$(echo "$CURRENT_TUMOR" | sed 's/_tumor_final.bam//')
echo "Patient ID: $patient_id"

# Find normal file
NORMAL_BAM="${patient_id}_normal_final.bam"

if [ ! -f "$NORMAL_BAM" ]; then
    echo "ERROR: No matching normal: $NORMAL_BAM"
    exit 1
fi

echo "Found pair:"
echo "  Tumor: $CURRENT_TUMOR"
echo "  Normal: $NORMAL_BAM"

# Check BAM integrity
echo "Checking BAM files..."
samtools quickcheck "$CURRENT_TUMOR" || { echo "ERROR: Tumor BAM corrupted"; exit 1; }
samtools quickcheck "$NORMAL_BAM" || { echo "ERROR: Normal BAM corrupted"; exit 1; }

# Index if needed
if [ ! -f "${CURRENT_TUMOR}.bai" ]; then
    echo "Indexing tumor..."
    samtools index -@ ${SLURM_CPUS_PER_TASK} "$CURRENT_TUMOR"
fi

if [ ! -f "${NORMAL_BAM}.bai" ]; then
    echo "Indexing normal..."
    samtools index -@ ${SLURM_CPUS_PER_TASK} "$NORMAL_BAM"
fi

# Run snp-pileup
echo "Running snp-pileup for $patient_id..."
pileup_output="${PILEUP_DIR}/${patient_id}.pileup"
log_output="${LOG_DIR}/${patient_id}_pileup.log"

snp-pileup -q 25 -Q 20 -r 10 \
    "${VCF_DIR}/out_sorted.vcf" \
    "$pileup_output" \
    "$NORMAL_BAM" \
    "$CURRENT_TUMOR" > "$log_output" 2>&1

if [ $? -eq 0 ] && [ -s "$pileup_output" ]; then
    echo "✓ SNP pileup completed"
    echo "  Size: $(du -h $pileup_output | cut -f1)"
    
    # Sort pileup
    echo "Sorting pileup..."
    temp_file="${PILEUP_DIR}/${patient_id}.sorted.tmp"
    {
        head -n 1 "$pileup_output"
        tail -n +2 "$pileup_output" | sort -t$'\t' -k1,1V -k2,2n
    } > "$temp_file"
    
    if [ $? -eq 0 ] && [ -s "$temp_file" ]; then
        mv "$temp_file" "${PILEUP_DIR}/${patient_id}.sorted.pileup"
        rm "$pileup_output"
        echo "✓ Sorted pileup created"
    else
        echo "✗ Error sorting"
        exit 1
    fi
else
    echo "✗ SNP pileup failed"
    exit 1
fi

echo ""
echo "========================================="
echo "Completed: $patient_id"
echo "========================================="
date

#!/bin/bash
#SBATCH --job-name=bam_process
#SBATCH --output=logs/bam_process_%A_%a.out
#SBATCH --error=logs/bam_process_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --array=0-3  # Process 4 pairs at a time

# BAM Processing Script for Meluxina
# Adds read groups, filters, removes duplicates
# For files with T_ (tumor) and B_ (blood/normal) prefixes

set -e

echo "Starting BAM processing on Meluxina..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
date

# Load required modules
module purge
module load SAMtools
module load Picard
module load Java

# Set directories
BAM_DIR=/home/users/u103499/Project/grp1
OUTPUT_DIR=/home/users/u103499/Project/processed_bams

# Create subdirectories
mkdir -p ${OUTPUT_DIR}/intermediate_files
mkdir -p ${OUTPUT_DIR}/final_bams
mkdir -p ${OUTPUT_DIR}/metrics
mkdir -p logs

# Get list of tumor BAM files (T_ prefix)
cd ${BAM_DIR}
TUMOR_BAMS=($(ls T_*.final.bam 2>/dev/null | sort))

# Calculate which tumor file to process based on array task ID
CURRENT_TUMOR=${TUMOR_BAMS[${SLURM_ARRAY_TASK_ID}]}

if [ -z "$CURRENT_TUMOR" ]; then
    echo "No tumor BAM file assigned to array task ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

echo "Processing tumor BAM: $CURRENT_TUMOR"

# Extract patient ID (e.g., T_24RV18.final.bam -> 24RV18)
patient_id=$(echo "$CURRENT_TUMOR" | sed 's/T_\(.*\)\.final\.bam/\1/')
echo "Patient ID: $patient_id"

# Find corresponding normal (blood) file
NORMAL_BAM="B_${patient_id}.final.bam"

if [ ! -f "$NORMAL_BAM" ]; then
    echo "ERROR: No matching normal file found: $NORMAL_BAM"
    exit 1
fi

echo "Found matching pair:"
echo "  Tumor: $CURRENT_TUMOR"
echo "  Normal: $NORMAL_BAM"
echo ""

# Process tumor sample
echo "========================================="
echo "Processing TUMOR: $CURRENT_TUMOR"
echo "========================================="

# Read group information for tumor
TUMOR_RGID="ID:${patient_id}_tumor"
TUMOR_RGLB="LB:WGS_${patient_id}_tumor"
TUMOR_RGSM="SM:${patient_id}_tumor"
TUMOR_RGPL="PL:ILLUMINA"
TUMOR_RGPU="PU:${patient_id}_tumor"

echo "Tumor Read Groups:"
echo "  $TUMOR_RGID"
echo "  $TUMOR_RGLB"
echo "  $TUMOR_RGSM"
echo "  $TUMOR_RGPL"
echo "  $TUMOR_RGPU"
echo ""

# Add read groups to tumor
echo "Adding read groups to tumor BAM..."
picard AddOrReplaceReadGroups \
    I=${BAM_DIR}/${CURRENT_TUMOR} \
    O=${OUTPUT_DIR}/intermediate_files/${patient_id}_tumor_RG.bam \
    RGID=${TUMOR_RGID#ID:} \
    RGLB=${TUMOR_RGLB#LB:} \
    RGPL=${TUMOR_RGPL#PL:} \
    RGPU=${TUMOR_RGPU#PU:} \
    RGSM=${TUMOR_RGSM#SM:} \
    VALIDATION_STRINGENCY=LENIENT

# Filter tumor BAM
echo "Filtering tumor BAM..."
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 2048 -q 10 \
    ${OUTPUT_DIR}/intermediate_files/${patient_id}_tumor_RG.bam | \
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 256 | \
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 4 \
    > ${OUTPUT_DIR}/intermediate_files/${patient_id}_tumor_filtered.bam

# Remove duplicates from tumor
echo "Removing duplicates from tumor BAM..."
picard MarkDuplicates \
    INPUT=${OUTPUT_DIR}/intermediate_files/${patient_id}_tumor_filtered.bam \
    OUTPUT=${OUTPUT_DIR}/final_bams/${patient_id}_tumor_final.bam \
    REMOVE_DUPLICATES=true \
    METRICS_FILE=${OUTPUT_DIR}/metrics/${patient_id}_tumor_dup_metrics.txt \
    VALIDATION_STRINGENCY=LENIENT

# Index tumor BAM
echo "Indexing tumor BAM..."
samtools index -@ ${SLURM_CPUS_PER_TASK} ${OUTPUT_DIR}/final_bams/${patient_id}_tumor_final.bam

echo "Tumor processing completed"
echo ""

# Process normal sample
echo "========================================="
echo "Processing NORMAL: $NORMAL_BAM"
echo "========================================="

# Read group information for normal
NORMAL_RGID="ID:${patient_id}_normal"
NORMAL_RGLB="LB:WGS_${patient_id}_normal"
NORMAL_RGSM="SM:${patient_id}_normal"
NORMAL_RGPL="PL:ILLUMINA"
NORMAL_RGPU="PU:${patient_id}_normal"

echo "Normal Read Groups:"
echo "  $NORMAL_RGID"
echo "  $NORMAL_RGLB"
echo "  $NORMAL_RGSM"
echo "  $NORMAL_RGPL"
echo "  $NORMAL_RGPU"
echo ""

# Add read groups to normal
echo "Adding read groups to normal BAM..."
picard AddOrReplaceReadGroups \
    I=${BAM_DIR}/${NORMAL_BAM} \
    O=${OUTPUT_DIR}/intermediate_files/${patient_id}_normal_RG.bam \
    RGID=${NORMAL_RGID#ID:} \
    RGLB=${NORMAL_RGLB#LB:} \
    RGPL=${NORMAL_RGPL#PL:} \
    RGPU=${NORMAL_RGPU#PU:} \
    RGSM=${NORMAL_RGSM#SM:} \
    VALIDATION_STRINGENCY=LENIENT

# Filter normal BAM
echo "Filtering normal BAM..."
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 2048 -q 10 \
    ${OUTPUT_DIR}/intermediate_files/${patient_id}_normal_RG.bam | \
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 256 | \
samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 4 \
    > ${OUTPUT_DIR}/intermediate_files/${patient_id}_normal_filtered.bam

# Remove duplicates from normal
echo "Removing duplicates from normal BAM..."
picard MarkDuplicates \
    INPUT=${OUTPUT_DIR}/intermediate_files/${patient_id}_normal_filtered.bam \
    OUTPUT=${OUTPUT_DIR}/final_bams/${patient_id}_normal_final.bam \
    REMOVE_DUPLICATES=true \
    METRICS_FILE=${OUTPUT_DIR}/metrics/${patient_id}_normal_dup_metrics.txt \
    VALIDATION_STRINGENCY=LENIENT

# Index normal BAM
echo "Indexing normal BAM..."
samtools index -@ ${SLURM_CPUS_PER_TASK} ${OUTPUT_DIR}/final_bams/${patient_id}_normal_final.bam

echo "Normal processing completed"
echo ""

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f ${OUTPUT_DIR}/intermediate_files/${patient_id}_tumor_RG.bam
rm -f ${OUTPUT_DIR}/intermediate_files/${patient_id}_tumor_filtered.bam
rm -f ${OUTPUT_DIR}/intermediate_files/${patient_id}_normal_RG.bam
rm -f ${OUTPUT_DIR}/intermediate_files/${patient_id}_normal_filtered.bam

echo "========================================="
echo "Processing completed for: $patient_id"
echo "Final BAMs:"
echo "  Tumor: ${OUTPUT_DIR}/final_bams/${patient_id}_tumor_final.bam"
echo "  Normal: ${OUTPUT_DIR}/final_bams/${patient_id}_normal_final.bam"
echo "========================================="
date

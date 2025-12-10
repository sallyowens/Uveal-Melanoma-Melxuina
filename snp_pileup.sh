#!/bin/bash
#SBATCH --job-name=snp_pileup
#SBATCH --account=p201093
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --output=/home/users/u103499/Project/logs/snp_pileup_%j.out
#SBATCH --error=/home/users/u103499/Project/logs/snp_pileup_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=64G

set -e

# === PATH CONFIGURATION ===
PROJECT_DIR=${HOME}/Project
BAM_DIR=${PROJECT_DIR}/grp1
VCF_DIR=${PROJECT_DIR}/FACETS/VCF
PILEUP_DIR=${PROJECT_DIR}/FACETS/pileups
LOG_DIR=${PROJECT_DIR}/FACETS/logs

# Initialize conda
source ${PROJECT_DIR}/software/miniconda3/etc/profile.d/conda.sh
conda activate wgs_analysis

# Create directories
mkdir -p ${VCF_DIR} ${PILEUP_DIR} ${LOG_DIR}

echo "=== SNP Pileup for FACETS WGS ==="
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${PILEUP_DIR}"
date
echo ""

# Download and prepare VCF if not exists
if [ ! -f "${VCF_DIR}/out_sorted.vcf" ]; then
    echo "Downloading common SNPs VCF (this is large, ~2 GB)..."
    
    if [ ! -f "${VCF_DIR}/00-common_all.vcf.gz" ]; then
        wget -P ${VCF_DIR} \
            https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
    fi
    
    echo "Processing VCF file (this will take 10-15 minutes)..."
    gunzip -c ${VCF_DIR}/00-common_all.vcf.gz | \
    awk 'BEGIN{OFS="\t"} 
         /^#/ {print; next}
         /^[0-9XY]+\t/ && ($6 == "." || $6 >= 30) {
             if($1 ~ /^(chr)?[0-9XY]+$/) print
         }' > ${VCF_DIR}/temp_filtered.vcf
    
    echo "Sorting VCF..."
    (grep '^#' ${VCF_DIR}/temp_filtered.vcf; 
     grep -v '^#' ${VCF_DIR}/temp_filtered.vcf | sort -k1,1V -k2,2n) \
    > ${VCF_DIR}/out_sorted.vcf
    
    rm ${VCF_DIR}/temp_filtered.vcf
    
    echo "✓ VCF prepared"
else
    echo "✓ VCF file already exists"
fi

echo "VCF lines: $(wc -l < ${VCF_DIR}/out_sorted.vcf)"
echo ""

cd ${BAM_DIR}

# Match tumor-normal pairs and run pileup
processed=0
failed=0

for tumor_bam in T_*.final.bam; do
    if [ -f "$tumor_bam" ]; then
        # Extract patient ID (T_24RV18 -> 24RV18)
        patient_id=$(echo ${tumor_bam} | sed 's/T_\(.*\)\.final\.bam/\1/')
        normal_bam="B_${patient_id}.final.bam"
        
        echo "=========================================="
        echo "Processing: ${patient_id}"
        echo "=========================================="
        echo "Tumor: ${tumor_bam}"
        echo "Normal: ${normal_bam}"
        
        if [ -f "${normal_bam}" ]; then
            # Ensure BAMs are indexed
            if [ ! -f "${tumor_bam}.bai" ]; then
                echo "Indexing tumor BAM..."
                samtools index ${tumor_bam}
            fi
            
            if [ ! -f "${normal_bam}.bai" ]; then
                echo "Indexing normal BAM..."
                samtools index ${normal_bam}
            fi
            
            # Run snp-pileup
            echo "Running snp-pileup (this will take several hours for WGS)..."
            echo "Start time: $(date)"
            
            timeout 43200 snp-pileup \
                -q 25 -Q 20 -r 10 \
                "${VCF_DIR}/out_sorted.vcf" \
                "${PILEUP_DIR}/${patient_id}.pileup" \
                "${normal_bam}" \
                "${tumor_bam}" \
                > "${LOG_DIR}/${patient_id}_pileup.log" 2>&1
            
            exit_code=$?
            
            if [ $exit_code -eq 0 ]; then
                # Sort pileup
                echo "Sorting pileup file..."
                {
                    head -n 1 "${PILEUP_DIR}/${patient_id}.pileup"
                    tail -n +2 "${PILEUP_DIR}/${patient_id}.pileup" | sort -k1,1V -k2,2n
                } > "${PILEUP_DIR}/${patient_id}.sorted.pileup"
                
                rm "${PILEUP_DIR}/${patient_id}.pileup"
                
                pileup_size=$(du -h "${PILEUP_DIR}/${patient_id}.sorted.pileup" | cut -f1)
                echo "✓ Completed: ${patient_id} (size: ${pileup_size})"
                processed=$((processed + 1))
            elif [ $exit_code -eq 124 ]; then
                echo "✗ Error: Timed out after 12 hours"
                failed=$((failed + 1))
            else
                echo "✗ Error: snp-pileup failed (exit code: ${exit_code})"
                echo "  Check log: ${LOG_DIR}/${patient_id}_pileup.log"
                failed=$((failed + 1))
            fi
        else
            echo "✗ Normal file not found: ${normal_bam}"
            failed=$((failed + 1))
        fi
        echo "End time: $(date)"
        echo ""
    fi
done

# Generate summary
{
    echo "SNP Pileup Summary"
    echo "=================="
    echo "Date: $(date)"
    echo ""
    echo "Results:"
    echo "- Successfully processed: ${processed}"
    echo "- Failed: ${failed}"
    echo ""
    echo "Pileup files:"
    ls -lh ${PILEUP_DIR}/*.sorted.pileup 2>/dev/null || echo "No pileup files generated"
    echo ""
    echo "Next step: Run FACETS R analysis"
} > ${PILEUP_DIR}/pileup_summary.txt

echo ""
echo "=== SNP Pileup Complete! ==="
echo "Successfully processed: ${processed}"
echo "Failed: ${failed}"
echo "Results in: ${PILEUP_DIR}"
echo "Summary: ${PILEUP_DIR}/pileup_summary.txt"
echo ""
date

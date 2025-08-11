#!/bin/bash
#SBATCH --job-name=rmats_siz1
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00

# Activate environment
mamba activate rmats

# Define paths for rMATS analysis
GTF_FILE=/global/scratch/users/enricocalvane/minRNAseq/as/Arabidopsis_thaliana.TAIR10.61.gtf
ALIGNED_DIR=/global/scratch/users/enricocalvane/minRNAseq/raw_reads/aligned_bams
OUTPUT_BASE=/global/scratch/users/enricocalvane/minRNAseq/as/rmats_output
TMP_DIR=/global/scratch/users/enricocalvane/minRNAseq/as/rmats_tmp

# Create output directories
mkdir -p $OUTPUT_BASE $TMP_DIR

# Create BAM list files for each experimental group
echo "Creating BAM list files..."

# Wild type samples (2 replicates)
echo "${ALIGNED_DIR}/WT_1.unique.sorted.bam,${ALIGNED_DIR}/WT_2.unique.sorted.bam" > $OUTPUT_BASE/wt_bams.txt

# siz1 mutant samples (2 replicates)
echo "${ALIGNED_DIR}/siz1_2_1.unique.sorted.bam,${ALIGNED_DIR}/siz1_2_2.unique.sorted.bam" > $OUTPUT_BASE/siz1_bams.txt

# SIZ1 overexpression line 1 (SIZ1_1 + SIZ1_2)
echo "${ALIGNED_DIR}/SIZ1_1.unique.sorted.bam,${ALIGNED_DIR}/SIZ1_2.unique.sorted.bam" > $OUTPUT_BASE/siz1_oe1_bams.txt

# SIZ1 overexpression line 2 (SIZ1_3 + SIZ1_4)
echo "${ALIGNED_DIR}/SIZ1_3.unique.sorted.bam,${ALIGNED_DIR}/SIZ1_4.unique.sorted.bam" > $OUTPUT_BASE/siz1_oe2_bams.txt

# Verify BAM files exist
echo "Verifying BAM files exist..."
for bam_list in $OUTPUT_BASE/*_bams.txt; do
    echo "Checking files in $(basename $bam_list):"
    while IFS=',' read -ra BAMS; do
        for bam in "${BAMS[@]}"; do
            if [[ -f "$bam" ]]; then
                echo "  ✓ Found: $(basename $bam)"
            else
                echo "  ✗ Missing: $(basename $bam)"
                exit 1
            fi
        done
    done < "$bam_list"
done

# Function to run rMATS analysis
run_rmats() {
    local comparison=$1
    local b1_file=$2
    local b2_file=$3
    local output_dir=${OUTPUT_BASE}/${comparison}
    local tmp_dir=${TMP_DIR}/${comparison}
    
    echo "=========================================="
    echo "Starting rMATS analysis for comparison: $comparison at $(date)"
    echo "Control group (b1): $(cat $b1_file)"
    echo "Treatment group (b2): $(cat $b2_file)"
    echo "=========================================="
    
    mkdir -p $output_dir $tmp_dir
    
    python $RMATS_PATH \
        --b1 $b1_file \
        --b2 $b2_file \
        --gtf $GTF_FILE \
        --od $output_dir \
        --tmp $tmp_dir \
        -t paired \
        --readLength 150 \
        --nthread $SLURM_NTASKS \
        --tstat $SLURM_NTASKS \
        --libType fr-unstranded \
        --task both \
        --novelSS \
        --individual-counts
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Completed rMATS analysis for $comparison at $(date)"
        echo "Output directory: $output_dir"
    else
        echo "✗ Failed rMATS analysis for $comparison at $(date)"
        exit 1
    fi
}

# Run all three comparisons
echo "=========================================="
echo "Starting rMATS comparative analysis at $(date)"
echo "Total comparisons to perform: 3"
echo "=========================================="

# Comparison 1: WT vs siz1 mutant
echo "Running comparison 1/3: WT vs siz1 mutant..."
run_rmats "wt_vs_siz1" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/siz1_bams.txt

# Comparison 2: WT vs SIZ1_OE1
echo "Running comparison 2/3: WT vs SIZ1_OE1..."
run_rmats "wt_vs_siz1_oe1" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/siz1_oe1_bams.txt

# Comparison 3: WT vs SIZ1_OE2
echo "Running comparison 3/3: WT vs SIZ1_OE2..."
run_rmats "wt_vs_siz1_oe2" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/siz1_oe2_bams.txt

echo "=========================================="
echo "All rMATS analyses completed successfully at $(date)"
echo "=========================================="

# Summary of output directories
echo "Output summary:"
echo "- WT vs siz1 mutant: $OUTPUT_BASE/wt_vs_siz1/"
echo "- WT vs SIZ1_OE1: $OUTPUT_BASE/wt_vs_siz1_oe1/"
echo "- WT vs SIZ1_OE2: $OUTPUT_BASE/wt_vs_siz1_oe2/"
echo ""
echo "Key files in each directory:"
echo "- SE.MATS.JC.txt (Skipped Exons)"
echo "- A5SS.MATS.JC.txt (Alternative 5' Splice Sites)"
echo "- A3SS.MATS.JC.txt (Alternative 3' Splice Sites)"
echo "- MXE.MATS.JC.txt (Mutually Exclusive Exons)"
echo "- RI.MATS.JC.txt (Retained Introns)"
echo ""
echo "PSI interpretation:"
echo "- Negative ΔPSI: More splicing variant in mutant/OE vs WT"
echo "- Positive ΔPSI: More splicing variant in WT vs mutant/OE"

#!/bin/bash
#SBATCH --job-name=rmats_gbpl3
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00

#Activate environment
mamba activate rmats

# Define paths for rMATS analysis
GTF_FILE=/global/scratch/users/enricocalvane/heejae_as/reference/TAIR10_fixed.gtf
ALIGNED_DIR=/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams
OUTPUT_BASE=/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/rmats_output
TMP_DIR=/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/rmats_tmp

# Create output directories
mkdir -p $OUTPUT_BASE $TMP_DIR

# Create BAM list files
# For wild type samples
echo "${ALIGNED_DIR}/SRR18516933_wild_type.unique.sorted.bam,${ALIGNED_DIR}/SRR18516934_wild_type.unique.sorted.bam,${ALIGNED_DIR}/SRR18516935_wild_type.unique.sorted.bam" > $OUTPUT_BASE/wt_bams.txt

# For gbpl3 samples
echo "${ALIGNED_DIR}/SRR18516930_gbpl3_3.unique.sorted.bam,${ALIGNED_DIR}/SRR18516931_gbpl3_2.unique.sorted.bam,${ALIGNED_DIR}/SRR18516932_gbpl3_1.unique.sorted.bam" > $OUTPUT_BASE/gbpl3_bams.txt

# Function to run rMATS analysis
run_rmats() {
    local comparison=$1
    local b1_file=$2
    local b2_file=$3
    local output_dir=${OUTPUT_BASE}/${comparison}
    local tmp_dir=${TMP_DIR}/${comparison}

    echo "Starting rMATS analysis for comparison: $comparison at $(date)"
    
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

    echo "Completed rMATS analysis for $comparison at $(date)"
}

# Run the comparison
echo "Starting rMATS comparison at $(date)"

echo "Running gbpl3 vs WT comparison..."
run_rmats "gbpl3_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/gbpl3_bams.txt

echo "rMATS analysis completed at $(date)"

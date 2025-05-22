#!/bin/bash
#SBATCH --job-name=rmats_pnet4
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00

# Activate conda environment and set up library paths
source ~/.bashrc
conda activate rmats_env
export LD_LIBRARY_PATH=/global/home/users/enricocalvane/.conda/envs/rmats_env/lib:$LD_LIBRARY_PATH

# Define paths for rMATS analysis
RMATS_PATH=/global/scratch/users/enricocalvane/rmats_turbo/rmats-turbo/rmats.py
GTF_FILE=/global/scratch/users/enricocalvane/heejae_as/reference/TAIR10_fixed.gtf
ALIGNED_DIR=/global/scratch/users/enricocalvane/heejae_as/pnet4_data/pnet4_RNAseq/aligned_bams
OUTPUT_BASE=/global/scratch/users/enricocalvane/heejae_as/pnet4_data/rmats_output
TMP_DIR=/global/scratch/users/enricocalvane/heejae_as/pnet4_data/rmats_tmp

# Create output directories
mkdir -p $OUTPUT_BASE $TMP_DIR

# Create BAM list files - using the exact filenames from your directory
# For wild type samples (Col_RT replicates)
echo "${ALIGNED_DIR}/Col_RT1.unique.sorted.bam,${ALIGNED_DIR}/Col_RT2.unique.sorted.bam" > $OUTPUT_BASE/wt_bams.txt

# For pnet4 samples
echo "${ALIGNED_DIR}/pnet4_RT1.unique.sorted.bam,${ALIGNED_DIR}/pnet4_RT2.unique.sorted.bam" > $OUTPUT_BASE/pnet4_bams.txt

# Function to run rMATS analysis (unchanged from original)
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
echo "Running pnet4 vs WT comparison..."
run_rmats "pnet4_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/pnet4_bams.txt
echo "rMATS analysis completed at $(date)"

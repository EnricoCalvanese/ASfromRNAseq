#!/bin/bash
#SBATCH --job-name=rmats_pnet1
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
ALIGNED_DIR=/global/scratch/users/enricocalvane/heejae_as/pnet1_data/aligned_bams
OUTPUT_BASE=/global/scratch/users/enricocalvane/heejae_as/pnet1_data/rmats_output
TMP_DIR=/global/scratch/users/enricocalvane/heejae_as/pnet1_data/rmats_tmp

# Create output directories
mkdir -p $OUTPUT_BASE $TMP_DIR

# Create BAM list files - carefully matching your file naming pattern
# For wild type samples (Col-0 replicates)
echo "${ALIGNED_DIR}/Col-0-1_FRRB192015921-1a.unique.sorted.bam,${ALIGNED_DIR}/Col-0-2_FRRB192015922-1a.unique.sorted.bam,${ALIGNED_DIR}/Col-0-3_FRRB192015923-1a.unique.sorted.bam" > $OUTPUT_BASE/wt_bams.txt

# For pnet1 samples (t3 replicates)
echo "${ALIGNED_DIR}/t3-1_FRRB192015924-1a.unique.sorted.bam,${ALIGNED_DIR}/t3-2_FRRB192015925-1a.unique.sorted.bam,${ALIGNED_DIR}/t3-3_FRRB192015926-1a.unique.sorted.bam" > $OUTPUT_BASE/pnet1_bams.txt

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
echo "Running pnet1 vs WT comparison..."
run_rmats "pnet1_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/pnet1_bams.txt
echo "rMATS analysis completed at $(date)"

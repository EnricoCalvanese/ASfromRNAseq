#!/bin/bash
#SBATCH --job-name=rmats_analysis
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Define paths
GTF_FILE=/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf
OUTPUT_BASE=/global/scratch/users/enricocalvane/heejae_as/crwn_data/rmats_output

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_BASE

# Function to run rMATS for a comparison
run_rmats() {
    local comparison=$1
    local b1_file=$2
    local b2_file=$3
    local output_dir=${OUTPUT_BASE}/${comparison}

    echo "Starting rMATS analysis for comparison: $comparison at $(date)"
    
    # Create the output directory for this comparison
    mkdir -p $output_dir
    
    rmats.py \
        --b1 $b1_file \
        --b2 $b2_file \
        --gtf $GTF_FILE \
        --od $output_dir \
        -t single \
        --readLength 51 \
        --nthread $SLURM_NTASKS \
        --tstat $SLURM_NTASKS \
        --cstat 0.0001 \
        --libType fr-unstranded \
        --novelSS \
        --statoff

    echo "Completed rMATS analysis for $comparison at $(date)"
}

# First, create the BAM list files
echo "Creating BAM list files..."

# Wild type samples
echo "${ALIGNED_DIR}/SRR7657889_02_WT.sorted.bam,${ALIGNED_DIR}/SRR7657890_01_WT.sorted.bam,${ALIGNED_DIR}/SRR7657898_03_WT.sorted.bam" > $OUTPUT_BASE/wt_bams.txt

# CRWN1 samples
echo "${ALIGNED_DIR}/SRR7657893_04_crwn1.sorted.bam,${ALIGNED_DIR}/SRR7657894_05_crwn1.sorted.bam,${ALIGNED_DIR}/SRR7657895_06_crwn1.sorted.bam" > $OUTPUT_BASE/crwn1_bams.txt

# CRWN2 samples
echo "${ALIGNED_DIR}/SRR7657897_09_crwn2.sorted.bam,${ALIGNED_DIR}/SRR7657900_07_crwn2.sorted.bam,${ALIGNED_DIR}/SRR7657901_08_crwn2.sorted.bam" > $OUTPUT_BASE/crwn2_bams.txt

# CRWN4 samples
echo "${ALIGNED_DIR}/SRR7657891_12_crwn4.sorted.bam,${ALIGNED_DIR}/SRR7657892_11_crwn4.sorted.bam,${ALIGNED_DIR}/SRR7657896_10_crwn4.sorted.bam" > $OUTPUT_BASE/crwn4_bams.txt

# CRWN1/CRWN2 double mutant samples
echo "${ALIGNED_DIR}/SRR7657887_13_crwn1_crwn2.sorted.bam,${ALIGNED_DIR}/SRR7657888_15_crwn1_crwn2.sorted.bam,${ALIGNED_DIR}/SRR7657902_14_crwn1_crwn2.sorted.bam" > $OUTPUT_BASE/crwn1_crwn2_bams.txt

# CRWN1/CRWN4 double mutant samples
echo "${ALIGNED_DIR}/SRR7657885_18_crwn1_crwn4.sorted.bam,${ALIGNED_DIR}/SRR7657886_17_crwn1_crwn4.sorted.bam,${ALIGNED_DIR}/SRR7657899_16_crwn1_crwn4.sorted.bam" > $OUTPUT_BASE/crwn1_crwn4_bams.txt

# Run all comparisons
echo "Starting all rMATS comparisons at $(date)"

# CRWN1 vs WT
run_rmats "crwn1_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn1_bams.txt

# CRWN2 vs WT
run_rmats "crwn2_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn2_bams.txt

# CRWN4 vs WT
run_rmats "crwn4_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn4_bams.txt

# CRWN1/CRWN2 vs WT
run_rmats "crwn1_crwn2_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn1_crwn2_bams.txt

# CRWN1/CRWN4 vs WT
run_rmats "crwn1_crwn4_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn1_crwn4_bams.txt

echo "All rMATS comparisons completed at $(date)"

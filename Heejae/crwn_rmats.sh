#!/bin/bash
#SBATCH --job-name=rmats_unique
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Define our key directories and files
# The GTF file provides the gene model that rMATS uses to identify possible splicing events
GTF_FILE=/global/scratch/users/enricocalvane/heejae_as/reference/TAIR10_fixed.gtf

# The aligned_bams directory contains our uniquely mapped reads
ALIGNED_DIR=/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams

# This is where rMATS will store its analysis results
OUTPUT_BASE=/global/scratch/users/enricocalvane/heejae_as/crwn_data/rmats_output

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_BASE
mkdir -p $OUTPUT_BASE/tmp

# Create BAM list files for each condition, using the unique mapped reads
echo "Creating BAM list files for uniquely mapped reads..."

# Wild type samples (our control group)
echo "${ALIGNED_DIR}/SRR7657889_02_WT.unique.sorted.bam,${ALIGNED_DIR}/SRR7657890_01_WT.unique.sorted.bam,${ALIGNED_DIR}/SRR7657898_03_WT.unique.sorted.bam" > $OUTPUT_BASE/wt_bams.txt

# CRWN1 mutant samples
echo "${ALIGNED_DIR}/SRR7657893_04_crwn1.unique.sorted.bam,${ALIGNED_DIR}/SRR7657894_05_crwn1.unique.sorted.bam,${ALIGNED_DIR}/SRR7657895_06_crwn1.unique.sorted.bam" > $OUTPUT_BASE/crwn1_bams.txt

# CRWN2 mutant samples
echo "${ALIGNED_DIR}/SRR7657897_09_crwn2.unique.sorted.bam,${ALIGNED_DIR}/SRR7657900_07_crwn2.unique.sorted.bam,${ALIGNED_DIR}/SRR7657901_08_crwn2.unique.sorted.bam" > $OUTPUT_BASE/crwn2_bams.txt

# CRWN4 mutant samples
echo "${ALIGNED_DIR}/SRR7657891_12_crwn4.unique.sorted.bam,${ALIGNED_DIR}/SRR7657892_11_crwn4.unique.sorted.bam,${ALIGNED_DIR}/SRR7657896_10_crwn4.unique.sorted.bam" > $OUTPUT_BASE/crwn4_bams.txt

# CRWN1/CRWN2 double mutant samples
echo "${ALIGNED_DIR}/SRR7657887_13_crwn1_crwn2.unique.sorted.bam,${ALIGNED_DIR}/SRR7657888_15_crwn1_crwn2.unique.sorted.bam,${ALIGNED_DIR}/SRR7657902_14_crwn1_crwn2.unique.sorted.bam" > $OUTPUT_BASE/crwn1_crwn2_bams.txt

# CRWN1/CRWN4 double mutant samples
echo "${ALIGNED_DIR}/SRR7657885_18_crwn1_crwn4.unique.sorted.bam,${ALIGNED_DIR}/SRR7657886_17_crwn1_crwn4.unique.sorted.bam,${ALIGNED_DIR}/SRR7657899_16_crwn1_crwn4.unique.sorted.bam" > $OUTPUT_BASE/crwn1_crwn4_bams.txt

# Function to run rMATS analysis for each comparison
# This function handles all the parameter settings and execution for one comparison
run_rmats() {
    local comparison=$1      # Name of the comparison (e.g., "crwn1_vs_wt")
    local b1_file=$2        # Path to file containing BAM files for condition 1 (WT)
    local b2_file=$3        # Path to file containing BAM files for condition 2 (mutant)
    local output_dir=${OUTPUT_BASE}/${comparison}

    echo "Starting rMATS analysis for comparison: $comparison at $(date)"
    
    mkdir -p $output_dir
    
    # Run rMATS with parameters optimized for our data
    rmats.py \
        --b1 $b1_file \
        --b2 $b2_file \
        --gtf $GTF_FILE \
        --od $output_dir \
        -t single \
        --readLength 51 \
        --nthread $SLURM_NTASKS \
        --tstat $SLURM_NTASKS \
        --cstat 0.05 \
        --libType fr-unstranded \
        --tmp $OUTPUT_BASE/tmp

    echo "Completed rMATS analysis for $comparison at $(date)"
}

# Run all comparisons sequentially
echo "Starting all rMATS comparisons at $(date)"

echo "Running CRWN1 vs WT comparison..."
run_rmats "crwn1_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn1_bams.txt

echo "Running CRWN2 vs WT comparison..."
run_rmats "crwn2_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn2_bams.txt

echo "Running CRWN4 vs WT comparison..."
run_rmats "crwn4_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn4_bams.txt

echo "Running CRWN1/CRWN2 vs WT comparison..."
run_rmats "crwn1_crwn2_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn1_crwn2_bams.txt

echo "Running CRWN1/CRWN4 vs WT comparison..."
run_rmats "crwn1_crwn4_vs_wt" \
    $OUTPUT_BASE/wt_bams.txt \
    $OUTPUT_BASE/crwn1_crwn4_bams.txt

echo "All rMATS comparisons completed at $(date)"

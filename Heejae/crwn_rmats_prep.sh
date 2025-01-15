#!/bin/bash
#SBATCH --job-name=rmats_prep
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Define directories
ALIGNED_DIR=/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams
GTF_FILE=/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf
OUTPUT_BASE=/global/scratch/users/enricocalvane/heejae_as/crwn_data/rmats_output

# Create output directory
mkdir -p $OUTPUT_BASE

# Create separate directories for each comparison
mkdir -p $OUTPUT_BASE/{crwn1_vs_wt,crwn2_vs_wt,crwn4_vs_wt,crwn1_crwn2_vs_wt,crwn1_crwn4_vs_wt}

# Create BAM list files for each condition
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

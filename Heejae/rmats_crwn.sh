#!/bin/bash
#SBATCH --job-name=rMATS_CRWN
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Define our key directories and files
# The GTF file provides the gene model information for identifying splice events
GTF_PATH=/global/scratch/users/enricocalvane/heejae_as/reference/TAIR10_fixed.gtf

# Set up the output directory where rMATS will store its results
OUTPUT_FOLDER=/global/scratch/users/enricocalvane/heejae_as/crwn_data/rmats_output

# Create text files listing the BAM files for each condition
# First, let's make a directory to store these lists if it doesn't exist
mkdir -p $OUTPUT_FOLDER

# Create the WT (control) sample list
echo "/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams/SRR7657889_02_WT.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams/SRR7657890_01_WT.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams/SRR7657898_03_WT.unique.sorted.bam" > $OUTPUT_FOLDER/wt_samples.txt

# Create each mutant's sample list
echo "/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams/SRR7657893_04_crwn1.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams/SRR7657894_05_crwn1.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams/SRR7657895_06_crwn1.unique.sorted.bam" > $OUTPUT_FOLDER/crwn1_samples.txt

# Create a temporary directory for rMATS processing
mkdir -p $OUTPUT_FOLDER/tmp

cd /global/scratch/users/enricocalvane/software/Conda/envs/colabfold/rMATS

# Run rMATS for CRWN1 vs WT comparison
# Note that we're using single-end reads (unlike the example which used paired-end)
python rmats.py \
    --b1 $OUTPUT_FOLDER/wt_samples.txt \
    --b2 $OUTPUT_FOLDER/crwn1_samples.txt \
    --gtf $GTF_PATH \
    -t single \
    --readLength 51 \
    --nthread $SLURM_NTASKS \
    --od $OUTPUT_FOLDER/crwn1_vs_wt \
    --tmp $OUTPUT_FOLDER/tmp

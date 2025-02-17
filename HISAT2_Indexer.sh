#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=04:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Define input files
FASTA_FILE=/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
GTF_FILE=/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf
INDEX_PREFIX=/global/scratch/users/enricocalvane/hisat2_index/tair10

# Create output directory if it doesn't exist
mkdir -p $(dirname $INDEX_PREFIX)

# Load required module
module load bio/hisat2/2.2.1-gcc-11.4.0

# Extract splice sites and exons from GTF file
hisat2_extract_splice_sites.py $GTF_FILE > ${INDEX_PREFIX}_splice_sites.txt
hisat2_extract_exons.py $GTF_FILE > ${INDEX_PREFIX}_exons.txt

# Build the index using both splice sites and exons
hisat2-build \
    -p $SLURM_NTASKS \
    --ss ${INDEX_PREFIX}_splice_sites.txt \
    --exon ${INDEX_PREFIX}_exons.txt \
    $FASTA_FILE \
    $INDEX_PREFIX

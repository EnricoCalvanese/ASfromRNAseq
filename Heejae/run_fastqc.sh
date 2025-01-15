#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Load required modules
module load bio/fastqc/0.12.1-gcc-11.4.0
module load parallel

# Create output directory for FastQC results
mkdir -p fastqc_results

# Create a list of all fastq.gz files
find . -name "*.fastq.gz" > fastq_files.txt

# Run FastQC in parallel using all available cores
cat fastq_files.txt | parallel -j $SLURM_NTASKS \
    "fastqc {} -o fastqc_results"

# Create a summary report of all FastQC results
multiqc fastqc_results -o fastqc_results/multiqc

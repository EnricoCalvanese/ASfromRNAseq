#!/bin/bash
#SBATCH --job-name=HISAT2_indexer
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#directories
FASTA_FILE=/global/scratch/users/enricocalvane/RNAseq_imb2/Athaliana_447_TAIR10.fa
GTF_FILE=/global/scratch/users/enricocalvane/RNAseq_imb2/Araport11_GTF_genes_transposons.current.gtf
#run hisat2 build
$ hisat2-build -p $SLURM_NTASKS $FASTA_FILE $GTF_FILE

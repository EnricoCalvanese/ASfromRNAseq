#!/bin/bash
#SBATCH --job-name=HISAT
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#directories
INDEX_FOLDER=/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-index
FASTQ_FOLDER=/global/scratch/users/enricocalvane/RNAseq_imb2/CleanData
OUTPUT_FOLDER=/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment

cd $FASTQ_FOLDER
module load hisat2/2.1.0

#run HISAT2 alignment
hisat2 -x $INDEX_FOLDER/Araport11_GTF_genes_transposons.current.gtf -1 Col-0-2_FRRB190012050-1a_1.clean.fq.gz,Col-0-3_FRRB190012051-1a_1.clean.fq.gz,imb2-1_FRRB190012052-1a_1.clean.fq.gz,imb2-2_FRRB190012053-1a_1.clean.fq.gz,imb2-3_FRRB190012054-1a_1.clean.fq.gz -2 Col-0-2_FRRB190012050-1a_2.clean.fq.gz,Col-0-3_FRRB190012051-1a_2.clean.fq.gz,imb2-1_FRRB190012052-1a_2.clean.fq.gz,imb2-2_FRRB190012053-1a_2.clean.fq.gz,imb2-3_FRRB190012054-1a_2.clean.fq.gz -S $OUTPUT_FOLDER/Col-0-1_FRRB190012049-1a_hisat.sam -p $SLURM_NTASKS

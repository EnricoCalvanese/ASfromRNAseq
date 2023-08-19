#!/bin/bash
#SBATCH --job-name=Convert-Sort-to-BAM
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#directories
OUTPUT_FOLDER=/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment

cd $OUTPUT_FOLDER
module load samtools/1.8

#convert from SAM to BAM
samtools view -bS Col-0-1_FRRB190012049_hisat.sam > Col-0-1_FRRB190012049_hisat.bam
samtools view -bS Col-0-2_FRRB190012050_hisat.sam > Col-0-2_FRRB190012050_hisat.bam
samtools view -bS Col-0-3_FRRB190012051_hisat.sam > Col-0-3_FRRB190012051_hisat.bam
samtools view -bS imb2-1_FRRB190012052_hisat.sam > imb2-1_FRRB190012052_hisat.bam
samtools view -bS imb2-2_FRRB190012053_hisat.sam > imb2-2_FRRB190012053_hisat.bam
samtools view -bS imb2-3_FRRB190012054_hisat.sam > imb2-3_FRRB190012054_hisat.bam

#sort BAM files
samtools sort Col-0-1_FRRB190012049_hisat.bam -o Col-0-1_FRRB190012049_hisat_sorted.bam
samtools sort Col-0-2_FRRB190012050_hisat.bam -o Col-0-2_FRRB190012050_hisat_sorted.bam
samtools sort Col-0-3_FRRB190012051_hisat.bam -o Col-0-3_FRRB190012051_hisat_sorted.bam
samtools sort imb2-1_FRRB190012052_hisat.bam -o imb2-1_FRRB190012052_hisat_sorted.bam
samtools sort imb2-2_FRRB190012053_hisat.bam -o imb2-2_FRRB190012053_hisat_sorted.bam
samtools sort imb2-3_FRRB190012054_hisat.bam -o imb2-3_FRRB190012054_hisat_sorted.bam


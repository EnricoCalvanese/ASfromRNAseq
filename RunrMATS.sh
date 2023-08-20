#!/bin/bash
#SBATCH --job-name=rMATS
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#directories
COL_0_PATH=/global/scratch/users/enricocalvane/RNAseq_imb2/rMATS/Col-0.txt
##Col-0.txt is /global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment/Col-0-1_FRRB190012049_hisat.bam,/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment/Col-0-2_FRRB190012050_hisat.bam,/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment/Col-0-3_FRRB190012051_hisat.bam
imb2_PATH=/global/scratch/users/enricocalvane/RNAseq_imb2/rMATS/imb2.txt
##imb2 is /global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment/imb2-1_FRRB190012052_hisat.bam,/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment/imb2-2_FRRB190012053_hisat.bam,/global/scratch/users/enricocalvane/RNAseq_imb2/HISAT-alignment/imb2-3_FRRB190012054_hisat.bam
OUTPUT_FOLDER=/global/scratch/users/enricocalvane/RNAseq_imb2/rMATS
GTF_PATH=/global/scratch/users/enricocalvane/RNAseq_imb2/Araport11_GTF_genes_transposons.current.gtf

cd /global/scratch/users/enricocalvane/software/Conda/envs/colabfold/rMATS

#Run rMATS
python rmats.py --b1 $COL_0_PATH --b2 $imb2_PATH --gtf GTF_PATH -t paired --readLength 150 --nthread $SLURM_NTASKS --od $OUTPUT_FOLDER --tmp $OUTPUT_FOLDER/tmp

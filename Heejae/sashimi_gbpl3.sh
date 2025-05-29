#!/bin/bash
#SBATCH --job-name=sashimi_gbpl3
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00

# Activate conda environment
source ~/.bashrc
conda activate rmats_env

# Define paths
OUTDIR=/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/sashimi_plot_output
GFF3=/global/scratch/users/enricocalvane/heejae_as/reference/Arabidopsis_thaliana.TAIR10.61.gff3
COORDINATE="4:-:11860409:11866475:$GFF3"

# BAM files (WT = b1, gbpl3 = b2)
B1_BAMS="/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams/SRR18516933_wild_type.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams/SRR18516934_wild_type.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams/SRR18516935_wild_type.unique.sorted.bam"
B2_BAMS="/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams/SRR18516930_gbpl3_3.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams/SRR18516931_gbpl3_2.unique.sorted.bam,/global/scratch/users/enricocalvane/heejae_as/gbpl3_data/aligned_bams/SRR18516932_gbpl3_1.unique.sorted.bam"

# Create output directory
mkdir -p $OUTDIR

# Create grouping file
GROUPING_FILE=${OUTDIR}/grouping.gf
cat <<EOF > $GROUPING_FILE
WT: 1-3
gbpl3: 4-6
EOF

# Run rmats2sashimiplot
rmats2sashimiplot \
    -o $OUTDIR \
    --b1 $B1_BAMS \
    --b2 $B2_BAMS \
    -c $COORDINATE \
    --l1 WT \
    --l2 gbpl3 \
    --exon_s 1 \
    --intron_s 1 \
    --group-info $GROUPING_FILE \
    --font-size 10 \
    --fig-height 7 \
    --fig-width 10

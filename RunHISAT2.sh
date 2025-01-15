#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=08:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Load required modules
module load bio/hisat2/2.2.1-gcc-11.4.0
module load bio/samtools/1.17-gcc-11.4.0

# Define directories and index
INDEX_PREFIX=/global/scratch/users/enricocalvane/heejae_as/hisat2_index/tair10
FASTQ_DIR=/global/scratch/users/enricocalvane/heejae_as/crwn_data
OUTPUT_DIR=/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams
mkdir -p $OUTPUT_DIR

# Create a log directory for mapping statistics
mkdir -p ${OUTPUT_DIR}/logs

# Loop through all fastq files and perform alignment
for FASTQ in ${FASTQ_DIR}/*.fastq.gz
do
    # Extract sample name from full path
    SAMPLE=$(basename $FASTQ .fastq.gz)
    
    echo "Processing sample: $SAMPLE"
    
    # Perform alignment with HISAT2
    # The --dta flag is crucial for downstream splice site analysis
    # -p uses all available cores for parallel processing
    # --no-mixed and --no-discordant ensure consistent handling of reads
    hisat2 \
        -x $INDEX_PREFIX \
        -U $FASTQ \
        --dta \
        -p $SLURM_NTASKS \
        --no-mixed \
        --no-discordant \
        --new-summary \
        --summary-file ${OUTPUT_DIR}/logs/${SAMPLE}_alignment_summary.txt \
        2> ${OUTPUT_DIR}/logs/${SAMPLE}_hisat2.log \
        | samtools view -bS - \
        | samtools sort -@ $SLURM_NTASKS -o ${OUTPUT_DIR}/${SAMPLE}.sorted.bam -

    # Index the sorted BAM file
    samtools index ${OUTPUT_DIR}/${SAMPLE}.sorted.bam
    
    # Generate alignment statistics
    samtools flagstat ${OUTPUT_DIR}/${SAMPLE}.sorted.bam > ${OUTPUT_DIR}/logs/${SAMPLE}_flagstat.txt
done

# Create a summary of alignment rates
echo "Sample,Total Reads,Aligned Reads,Alignment Rate" > ${OUTPUT_DIR}/logs/alignment_summary.csv
for LOG in ${OUTPUT_DIR}/logs/*_alignment_summary.txt
do
    SAMPLE=$(basename $LOG _alignment_summary.txt)
    TOTAL=$(grep "Total reads" $LOG | cut -f2)
    ALIGNED=$(grep "Aligned exactly 1 time" $LOG | cut -f2)
    RATE=$(grep "Overall alignment rate" $LOG | cut -f2)
    echo "$SAMPLE,$TOTAL,$ALIGNED,$RATE" >> ${OUTPUT_DIR}/logs/alignment_summary.csv
done

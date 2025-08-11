#!/bin/bash
#SBATCH --job-name=hisat2_align_paired
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=HISAT_alignment.log

# Load required modules
module load bio/hisat2/2.2.1-gcc-11.4.0
module load bio/samtools/1.17-gcc-11.4.0

# Define directories and index
INDEX_PREFIX=/global/scratch/users/enricocalvane/minRNAseq/as/hisat2_index/tair10
FASTQ_DIR=/global/scratch/users/enricocalvane/minRNAseq/raw_reads
OUTPUT_DIR=/global/scratch/users/enricocalvane/minRNAseq/raw_reads/aligned_bams
mkdir -p $OUTPUT_DIR

# Create a log directory for mapping statistics
mkdir -p ${OUTPUT_DIR}/logs

# Process paired-end reads
# We look for read1 files and then find their matching read2 files
for READ1 in ${FASTQ_DIR}/*_1.fq
do
    # Get the base sample name by removing '_1.fq'
    BASE_SAMPLE=$(basename $READ1 _1.fq)
    
    # Construct the read2 filename
    READ2=${FASTQ_DIR}/${BASE_SAMPLE}_2.fq
    
    echo "Processing sample: $BASE_SAMPLE"
    echo "Read 1: $READ1"
    echo "Read 2: $READ2"
    
    # Perform paired-end alignment with HISAT2
    # Key differences for paired-end data:
    # -1 specifies the first read file
    # -2 specifies the second read file
    # --no-mixed and --no-discordant ensure proper paired-end alignment
    hisat2 \
        -x $INDEX_PREFIX \
        -1 $READ1 \
        -2 $READ2 \
        --dta \
        -p $SLURM_NTASKS \
        --no-mixed \
        --no-discordant \
        --new-summary \
        --summary-file ${OUTPUT_DIR}/logs/${BASE_SAMPLE}_alignment_summary.txt \
        2> ${OUTPUT_DIR}/logs/${BASE_SAMPLE}_hisat2.log \
        | samtools view -bS - \
        | samtools sort -@ $SLURM_NTASKS -o ${OUTPUT_DIR}/${BASE_SAMPLE}.sorted.bam -

    # Index the sorted BAM file
    samtools index ${OUTPUT_DIR}/${BASE_SAMPLE}.sorted.bam
    
    # Generate alignment statistics
    samtools flagstat ${OUTPUT_DIR}/${BASE_SAMPLE}.sorted.bam > ${OUTPUT_DIR}/logs/${BASE_SAMPLE}_flagstat.txt
done

# Create a summary of alignment rates
echo "Sample,Total Reads,Aligned Reads,Alignment Rate" > ${OUTPUT_DIR}/logs/alignment_summary.csv
for LOG in ${OUTPUT_DIR}/logs/*_alignment_summary.txt
do
    SAMPLE=$(basename $LOG _alignment_summary.txt)
    TOTAL=$(grep "Total pairs" $LOG | cut -f2)
    ALIGNED=$(grep "Aligned concordantly exactly 1 time" $LOG | cut -f2)
    RATE=$(grep "Overall alignment rate" $LOG | cut -f2)
    echo "$SAMPLE,$TOTAL,$ALIGNED,$RATE" >> ${OUTPUT_DIR}/logs/alignment_summary.csv
done

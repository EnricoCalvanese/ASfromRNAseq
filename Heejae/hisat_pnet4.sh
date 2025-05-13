
#!/bin/bash
#SBATCH --job-name=hisat2_pnet4
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=08:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Load required modules for sequence alignment and BAM file processing
module load bio/hisat2/2.2.1-gcc-11.4.0
module load bio/samtools/1.17-gcc-11.4.0

# Define our working directories and reference index
# The index prefix points to the HISAT2 index we built for the Arabidopsis genome
INDEX_PREFIX=/global/scratch/users/enricocalvane/heejae_as/hisat2_index/tair10

# Set up input and output directories specific to the pnet4 dataset
FASTQ_DIR=/global/scratch/users/enricocalvane/heejae_as/pnet4_data/pnet4_RNAseq
OUTPUT_DIR=/global/scratch/users/enricocalvane/heejae_as/pnet4_data/pnet4_RNAseq/aligned_bams

# Create output directories for BAM files and alignment logs
mkdir -p $OUTPUT_DIR
mkdir -p ${OUTPUT_DIR}/logs

# Process each pair of paired-end reads
# We search for read1 files and construct the matching read2 filename
for READ1 in ${FASTQ_DIR}/*_1.fq.gz
do
    # Extract the base sample name by removing the '_1.fq.gz' suffix
    BASE_SAMPLE=$(basename $READ1 _1.fq.gz)
    
    # Construct the path to the matching read2 file
    READ2=${FASTQ_DIR}/${BASE_SAMPLE}_2.fq.gz
    
    # Log which sample we're processing to help with monitoring
    echo "Processing sample: $BASE_SAMPLE"
    echo "Read 1: $READ1"
    echo "Read 2: $READ2"
    
    # Perform paired-end alignment with HISAT2
    # This command aligns both read files simultaneously, maintaining pair information
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

    # Create an index for the sorted BAM file
    # This index is required for many downstream tools including rMATS
    samtools index ${OUTPUT_DIR}/${BASE_SAMPLE}.sorted.bam
    
    # Generate detailed alignment statistics
    # This helps us verify the quality of our alignments
    samtools flagstat ${OUTPUT_DIR}/${BASE_SAMPLE}.sorted.bam > ${OUTPUT_DIR}/logs/${BASE_SAMPLE}_flagstat.txt
done

# Create a comprehensive summary of alignment rates for all samples
# This helps us quickly verify that all samples aligned with similar efficiency
echo "Sample,Total Reads,Aligned Reads,Alignment Rate" > ${OUTPUT_DIR}/logs/alignment_summary.csv
for LOG in ${OUTPUT_DIR}/logs/*_alignment_summary.txt
do
    SAMPLE=$(basename $LOG _alignment_summary.txt)
    TOTAL=$(grep "Total pairs" $LOG | cut -f2)
    ALIGNED=$(grep "Aligned concordantly exactly 1 time" $LOG | cut -f2)
    RATE=$(grep "Overall alignment rate" $LOG | cut -f2)
    echo "$SAMPLE,$TOTAL,$ALIGNED,$RATE" >> ${OUTPUT_DIR}/logs/alignment_summary.csv
done

#!/bin/bash
#SBATCH --job-name=hisat2_align_paired
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=03:00:00
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

# Function to filter BAM files for uniquely mapped reads
filter_bam() {
    # Get input BAM file path
    input_bam=$1
    
    # Create output filename by inserting 'unique' before '.sorted.bam'
    # For example: sample.sorted.bam becomes sample.unique.sorted.bam
    output_bam=${input_bam/.sorted.bam/.unique.sorted.bam}
    
    # Create a temporary filename to avoid overwriting our input while processing
    temp_bam=${input_bam/.sorted.bam/.temp.sorted.bam}
    
    echo "Processing ${input_bam##*/}"
    echo "Started at $(date)"
    
    # Filter for uniquely mapped reads (MAPQ > 10)
    # The -bq 10 flag keeps only reads with mapping quality above 10
    # This effectively removes ambiguously mapped reads
    samtools view -bq 10 $input_bam | \
    samtools sort -@ $SLURM_NTASKS -o $temp_bam
    
    # Create index for the filtered BAM
    samtools index $temp_bam
    
    # Move the temporary file to its final name
    mv $temp_bam $output_bam
    mv ${temp_bam}.bai ${output_bam}.bai
    
    # Calculate and display mapping statistics
    echo "Generating mapping statistics for ${output_bam##*/}"
    samtools flagstat $output_bam > ${output_bam}.flagstat
    
    echo "Completed processing ${input_bam##*/} at $(date)"
}

# Process all sorted BAM files that don't already have 'unique' in their name
echo "Starting BAM filtering process at $(date)"
for bam in ${OUTPUT_DIR}/*.sorted.bam; do
    # Skip if this is already a unique BAM file
    if [[ $bam != *"unique"* ]]; then
        filter_bam $bam
    fi
done

echo "All BAM filtering completed at $(date)"

# Create a summary of the results
echo "Creating summary report..."
echo "File,Total Reads,Uniquely Mapped Reads,Mapping Rate" > ${OUTPUT_DIR}/unique_mapping_summary.csv
for stats in ${OUTPUT_DIR}/*.unique.sorted.bam.flagstat; do
    bam_file=$(basename ${stats%.flagstat})
    total=$(grep "total" $stats | awk '{print $1}')
    mapped=$(grep "mapped" $stats | head -n 1 | awk '{print $1}')
    rate=$(awk "BEGIN {printf \"%.2f\", ($mapped/$total)*100}")
    echo "$bam_file,$total,$mapped,$rate%" >> ${OUTPUT_DIR}/unique_mapping_summary.csv
done

echo "Process complete. Check unique_mapping_summary.csv for mapping statistics."

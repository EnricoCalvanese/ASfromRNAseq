#!/bin/bash
#SBATCH --job-name=unique_map
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Load required modules
module load bio/samtools/1.17-gcc-11.4.0

# Define the directory containing our BAM files
ALIGNED_DIR=/global/scratch/users/enricocalvane/heejae_as/crwn_data/aligned_bams

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
for bam in ${ALIGNED_DIR}/*.sorted.bam; do
    # Skip if this is already a unique BAM file
    if [[ $bam != *"unique"* ]]; then
        filter_bam $bam
    fi
done

echo "All BAM filtering completed at $(date)"

# Create a summary of the results
echo "Creating summary report..."
echo "File,Total Reads,Uniquely Mapped Reads,Mapping Rate" > ${ALIGNED_DIR}/unique_mapping_summary.csv
for stats in ${ALIGNED_DIR}/*.unique.sorted.bam.flagstat; do
    bam_file=$(basename ${stats%.flagstat})
    total=$(grep "total" $stats | awk '{print $1}')
    mapped=$(grep "mapped" $stats | head -n 1 | awk '{print $1}')
    rate=$(awk "BEGIN {printf \"%.2f\", ($mapped/$total)*100}")
    echo "$bam_file,$total,$mapped,$rate%" >> ${ALIGNED_DIR}/unique_mapping_summary.csv
done

echo "Process complete. Check unique_mapping_summary.csv for mapping statistics."

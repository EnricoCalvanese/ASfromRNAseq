#!/bin/bash
#SBATCH --job-name=bootstrapping
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
# Set the input MSA file
input_msa="Curated.FL-ECLIPSE.msa.fasta"
output_tree="FL-ECLIPSE.tree"

# Ensure input file is in Unix format
dos2unix "$input_msa"

# Generate bootstrap replicates
seqboot << EOF
$input_msa
R
1000
Y
12345
Y
EOF

# Check if bootstrapped file was created
if [ ! -f "bootstrapped.fasta" ]; then
    echo "Error: Bootstrap file not created"
    exit 1
fi

# Split the bootstrapped file into individual files
csplit -z -b "%03d.fasta" -f replicate bootstrapped.fasta '/^>/' '{*}'

# Build trees for each replicate
for i in replicate*.fasta; do
    FastTree -nt "$i" > "$i.newick"
done

# Combine trees
cat replicate*.newick > all_trees.tre

# Generate consensus tree
consense << EOF
all_trees.tre
Y
EOF

# Rename the output to the desired output name
if [ -f outfile ]; then
    mv outfile "$output_tree"
else
    echo "Error: Consensus tree file not created"
    exit 1
fi

# Optional: Clean up intermediate files
rm replicate*.fasta replicate*.newick all_trees.tre bootstrapped.fasta

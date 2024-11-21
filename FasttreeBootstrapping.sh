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

# Set the input MSA file
input_msa="Curated.FL-ECLIPSE.msa.fasta"
output_tree="FL-ECLIPSE.tree"

# Explicitly copy input file to infile
cp "$input_msa" infile

# Run seqboot with explicit input
seqboot << EOF
D
J
R
1000
W
N
C
N
S
I
Y
12345
Y
EOF

# Check if bootstrapped file was created
if [ ! -f "outfile" ]; then
    echo "Error: Bootstrap file not created"
    exit 1
fi

# Rename outfile to bootstrapped.fasta
mv outfile bootstrapped.fasta

# Split the bootstrapped file
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

# Rename the output
if [ -f outfile ]; then
    mv outfile "$output_tree"
else
    echo "Error: Consensus tree file not created"
    exit 1
fi

# Optional: Clean up intermediate files
rm infile replicate*.fasta replicate*.newick all_trees.tre bootstrapped.fasta

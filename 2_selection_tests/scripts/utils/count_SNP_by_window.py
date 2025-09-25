#!/usr/bin/env python
# coding=utf-8
from Bio import SeqIO

####################################################################################
# Your input
window_size = 10
alignment_file = "example.fa"
output_file = f"snp_count_by_{window_size}bp_windows.txt"

####################################################################################
# Read the alignment file using SeqIO and call the other functions
def count_substitutions_from_file(alignment, window_size):
    aligned_sequences = []
    for record in SeqIO.parse(alignment, "fasta"):  # Adapt if your file is not in FASTA format, i.e: "clustal"
        aligned_sequences.append(record.seq)

    return count_snp(aligned_sequences, window_size)


def count_snp(sequences, window_size):
    num_snp = []
    alignment_length = len(sequences[0])  # Assuming all sequences have the same length

    # Split the alignment by window
    for start in range(0, alignment_length, window_size):
        end = min(start + window_size, alignment_length)
        window_sequences = [seq[start:end] for seq in sequences]

        # Count SNP for each window
        snp = 0
        for pos in range(len(window_sequences[0])):
            bases = set(window_sequences[i][pos].upper() for i in range(len(window_sequences)))  # upper is for eventual repeats
            # will be between 1 (conserved) and 4 (ATCG). Maybe more if indel, N,...
            if len(bases) > 1:
                snp += 1
        num_snp.append(snp)

    return num_snp, alignment_length


# Write output
def write_output(output, counts):
    with open(output, "w") as file:
        file.write("Window\tSNP\n")
        for i, mutations in enumerate(counts):
            file.write(f"{i+1}\t{mutations}\n")


####################################################################################
# Run everything
snp_counts, total_length = count_substitutions_from_file(alignment_file, window_size)
write_output(output_file, snp_counts)

print(f"Done! Found a total of {sum(snp_counts)} SNP in the alignment of {total_length} bases")

####################################################################################

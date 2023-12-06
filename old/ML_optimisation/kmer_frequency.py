from Bio import SeqIO
from Bio.Seq import Seq
import pandas
from collections import Counter
import time

start_time = time.time()
# Datas
path = "/Users/alaverre/Documents/Detecting_positive_selection/results/positive_selection/"
fasta_file = path + "human/CEBPA/Model/test_pos.fasta"
fasta_neg_file = path + "human/CEBPA/Model/test_neg.fasta"
output_file = path + "human/CEBPA/Model/test_kmer_frequency.csv"
output_file_hdf = path + "human/CEBPA/Model/test_kmer_frequency.hdf5"
kmer_file = path + "kmer.fa"  # File containing all possible k-mers
kmer_length = 10

def count_kmers(seqs, k, possible_kmers):
    kmer_freq = {}
    seq_ids = []
    observed_kmers = set()
    for seq_record in seqs:
        seq = str(seq_record.seq)
        seq_len = len(seq)
        seq_kmer_freq = Counter()  # Initialize Counter
        seq_ids.append(seq_record.id)
        for i in range(seq_len - k + 1):
            kmer = seq[i:i + k]
            # Get the reverse complement if needed
            if kmer not in possible_kmers:
                kmer = str(Seq(kmer).reverse_complement())

            seq_kmer_freq.update([kmer])
            observed_kmers.update([kmer])

        kmer_freq[seq_record.id] = dict(seq_kmer_freq)  # Store frequency dictionary for current sequence

    # Create DataFrame from result dictionary
    df = pandas.DataFrame(kmer_freq, index=list(observed_kmers), columns=seq_ids).fillna(0).T

    print("Adding missing k-mers...")
    missing_kmers = list(possible_kmers - observed_kmers)
    missing_df = pandas.DataFrame(0, index=df.index, columns=missing_kmers)
    df = pandas.concat([df, missing_df], axis=1)

    return df


print("Reading all non-redundant k-mers...")
all_kmers = [str(record.seq) for record in SeqIO.parse(kmer_file, "fasta")]

print("Reading sequences from FASTA file...")
sequences = SeqIO.parse(fasta_file, "fasta")
sequences_neg = SeqIO.parse(fasta_neg_file, "fasta")

print("Counting k-mers and get the result...")
result_pos = count_kmers(sequences, kmer_length, set(all_kmers))
result_pos['match'] = 1

result_neg = count_kmers(sequences_neg, kmer_length, set(all_kmers))
result_neg['match'] = 0
result = pandas.concat([result_pos, result_neg])

print("Writing output...")
result.to_hdf(output_file_hdf, key='data', mode='w', complevel=5)

"""
new_df = pd.read_hdf(output_file_hdf, key="data")
new_hdf = time.time()
print("Time spent HDF5:", new_hdf-end_hdf, "seconds")
"""
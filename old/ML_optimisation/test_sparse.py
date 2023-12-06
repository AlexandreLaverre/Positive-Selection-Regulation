from Bio import SeqIO
from Bio.Seq import Seq
import pandas
from collections import Counter
import time
from scipy.sparse import csr_matrix

start_time = time.time()
# Datas
path = "/Users/alaverre/Documents/Detecting_positive_selection/results/positive_selection/"
fasta_file = path + "human/CEBPA/sequences/test.fasta"
output_file = path + "human/CEBPA/sequences/test_kmer_frequency.csv"
output_file_hdf = path + "human/CEBPA/sequences/test_kmer_frequency.hdf5"
kmer_file = path + "kmer.fa"  # File containing all possible k-mers
kmer_length = 10


def count_kmers(seqs, k, possible_kmers):
    kmer_freq = {}
    seq_ids = []
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

        kmer_freq[seq_record.id] = dict(seq_kmer_freq)  # Store frequency dictionary for current sequence

    return kmer_freq, seq_ids


print("Reading all non-redundant k-mers...")
all_kmers = [str(record.seq) for record in SeqIO.parse(kmer_file, "fasta")]

print("Reading sequences from FASTA file...")
sequences = SeqIO.parse(fasta_file, "fasta")

print("Counting k-mers and get the result...")
result, IDs = count_kmers(sequences, kmer_length, set(all_kmers))

end_calcul = time.time()
print("Time spent calcul:", end_calcul-start_time, "seconds")

# Iterate over the result dictionary and collect the keys from each dictionary
observed_kmers = set()
for seq in result.values():
    observed_kmers.update(seq.keys())

# Create an empty sparse matrix with the appropriate shape
num_rows = len(IDs)
num_cols = len(observed_kmers)
data = []
row_indices = []
col_indices = []
# Iterate over each counter and populate the sparse matrix
for row, (seq_id, freq_dict) in enumerate(result.items()):
    for kmer, count in freq_dict.items():
        col = list(observed_kmers).index(kmer)  # Get the column index for the kmer
        data.append(count)
        row_indices.append(row)
        col_indices.append(col)

# Create the sparse matrix using CSR format
sparse_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(num_rows, num_cols))

end_sparse = time.time()
print("Time spent df:", end_sparse-end_calcul, "seconds")

# Create a DataFrame from the sparse matrix
df = pandas.DataFrame.sparse.from_spmatrix(sparse_matrix, columns=list(observed_kmers))
print(df)

end_df = time.time()
print("Time spent df:", end_df-end_sparse, "seconds")

print("Adding missing k-mers...")
missing_kmers = list(set(all_kmers) - observed_kmers)
missing_df = pandas.DataFrame(0, index=df.index, columns=missing_kmers)
df = pandas.concat([df, missing_df], axis=1)

end = time.time()
print("Time spent adding missing:", end-end_df, "seconds")
print("Total time spent :", end-start_time, "seconds")

"""
print("Writing output...")

df.to_hdf(output_file_hdf, key='data', mode='w', complevel=5)
end_hdf = time.time()
print("Time spent HDF5:", end_hdf-end, "seconds")


new_df = pd.read_hdf(output_file_hdf, key="data")
new_hdf = time.time()
print("Time spent HDF5:", new_hdf-end_hdf, "seconds")
"""
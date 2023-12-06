from Bio.Seq import Seq
import os
import pandas
import numpy


# Get number of substitutions per sequence
def get_sub_number(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Focal and ancestral sequences don't have the same length!")
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq_ref, seq_alt))


# Get ID of substitutions in a sequence
def get_sub_ids(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Focal and ancestral sequences don't have the same length!")
    sub_ids = []
    loc = 0
    for pos_ref, pos_alt in zip(seq_ref, seq_alt):
        if pos_ref != pos_alt:
            id = "pos" + str(loc) + ":" + pos_ref + "-" + pos_alt
            sub_ids.append(id)
        loc += 1

    return sub_ids


# Get SVM score for each kmer of a given model
def get_svm_dict(path_model):
    svm_dict = {}
    with open(path_model, 'r') as model:
        for i in model.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            kmer = i[0]
            svm_score = float(i[1])
            rev_kmer = str(Seq(kmer).reverse_complement())  # add the reverse complement kmer

            svm_dict[kmer] = svm_score
            svm_dict[rev_kmer] = svm_score


# Calculate SVM from sliding windows
def calculate_svm(seq, svm_dict, kmer_len):
    svm = 0
    for pos in range(len(seq) - kmer_len + 1):   # sliding window of kmer length
        kmer = seq[pos:pos+kmer_len]
        svm += svm_dict[kmer]   # sum of SVM for each kmer
    return round(svm, 7)


# Calculate delta SVM from sliding windows
def calculate_delta(seq_ref, seq_alt, svm_dict, kmer_len):
    delta_svm = 0
    for pos in range(len(seq_ref) - kmer_len + 1):   # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos+kmer_len]
        kmer_alt = seq_alt[pos:pos+kmer_len]
        if kmer_ref != kmer_alt:
            delta_svm += svm_dict[kmer_alt] - svm_dict[kmer_ref]  # sum of delta between sequences for each kmer
    return round(delta_svm, 7)


# Get the Substitution Matrix and the normalised one for each chromosome
def get_sub_matrix(path_matrix):
    submats = {}
    submats_norm = {}
    for file in os.listdir(path_matrix):
        if file.endswith('.txt'):
            chrom = file.strip('.txt')

            chrom_Table = pandas.read_table(path_matrix + file, sep=' ')
            chrom_Table.index = ['A', 'C', 'G', 'T']  # change row values
            numpy.fill_diagonal(chrom_Table.values, 0)  # assign 0 to diagonal
            chrom_SubMat = chrom_Table.to_dict('index')  # get dict by matrix rows
            submats[chrom] = chrom_SubMat
            submats_norm[chrom] = chrom_Table.div(chrom_Table.sum(axis=1), axis=0).to_dict(
                'index')  # get normalized dict by matrix rows

    if len(submats) == 0:
        raise ValueError("Substitution matrix not found!")

    return submats, submats_norm

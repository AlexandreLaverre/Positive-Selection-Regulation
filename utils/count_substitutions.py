#!/usr/bin/env python
# coding=utf-8
from Bio import SeqIO
import itertools
import argparse

####################################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("sp", help="Species name: human dog")
parser.add_argument("sample", help="Study name: Wilson Schmidt...")
parser.add_argument("TF", help="Transcription Factor: CEBPA CTCF ...")
parser.add_argument("peakType", help="NarrowPeaks or BroadPeaks")
args = parser.parse_args()

path = f"/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/positive_selection/" \
       f"{args.peakType}/{args.sp}/{args.sample}/{args.TF}/"

Ref_Path = f"{path}/sequences/filtered_ancestral_sequences.fa"
Target_Path = f"{path}/sequences/filtered_focal_sequences_upper.fa"
Output = open(f"{path}/sequences/focal_substitutions_stats.txt", "w")

nuc = ['A', 'T', 'C', 'G']
combinations = itertools.permutations(nuc, 2)

########################################################################################################################
# Get substitutions of each type per sequence
def get_sub_number(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Reference and target sequences don't have the same length!")

    subs = {"".join(sub): 0 for sub in combinations}

    for base1, base2 in zip(seq_ref, seq_alt):
        if base1 != base2:
            subs[f"{base1}{base2}"] += 1

    subs["Weak2Weak"] = subs["AT"] + subs["TA"]
    subs["Weak2Strong"] = subs["AC"] + subs["AG"] + subs["TC"] + subs["TG"]
    subs["Strong2Strong"] = subs["GC"] + subs["CG"]
    subs["Strong2Weak"] = subs["CA"] + subs["CT"] + subs["GA"] + subs["GT"]

    return subs

########################################################################################################################
# Get Reference sequences
ReferenceSeqs = SeqIO.to_dict(SeqIO.parse(open(Ref_Path), "fasta"))
if len(ReferenceSeqs) == 0:
    raise ValueError("Reference file is empty!")

# Get Focal sequences
TargetSeqs = SeqIO.to_dict(SeqIO.parse(open(Target_Path), "fasta"))
SeqIDs = TargetSeqs.keys()
if len(TargetSeqs) == 0:
    raise ValueError("Focal sequence file is empty!")

########################################################################################################################

permutation = ["".join(sub) for sub in combinations]
Output.write('\t'.join(["ID", "Length", "GC_Content", "NbSub"]) + '\t' + '\t'.join(permutation) + '\t' +
             '\t'.join(["Weak2Weak", "Weak2Strong", "Strong2Strong", "Strong2Weak"]) + "\n")

for ID in SeqIDs:
    ref_seq = str(ReferenceSeqs[ID].seq)
    target_seq = str(TargetSeqs[ID].seq)
    GC_content = str(float((target_seq.count('G') + target_seq.count('C'))) / len(target_seq) * 100)

    substitutions = get_sub_number(ref_seq, target_seq)
    substitution_number = list(substitutions.values())
    Nsub = str(sum(substitution_number))

    Output.write('\t'.join([ID, str(len(ref_seq)), GC_content, Nsub, '\t'.join(map(str, substitution_number))]) + "\n")


Output.close()

########################################################################################################################

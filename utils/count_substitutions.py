#!/usr/bin/env python
# coding=utf-8
import sys
from pathlib import Path
from Bio import SeqIO

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"

Ref_Path = path + sys.argv[1]
Target_Path = path + sys.argv[2]
Output = open(str(Path(Ref_Path).parent) + "/substitutions.txt", "w")


# Get substitutions of each type per sequence
def get_sub_number(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Reference and target sequences don't have the same length!")

    substitutions = {"Weak2Weak": 0, "Weak2Strong": 0, "Strong2Weak": 0, "Strong2Strong": 0}
    for base1, base2 in zip(seq_ref, seq_alt):
        if base1 != base2:
            base1_type = True if base1 in "AT" else False  # Weak == True, Strong == False
            base2_type = True if base2 in "AT" else False

            if base1_type:
                sub_type = "Weak2Weak" if base2_type else "Weak2Strong"
            else:
                sub_type = "Strong2Weak" if base2_type else "Strong2Strong"

            substitutions[sub_type] += 1

    return substitutions


# Get Reference sequences
ReferenceSeqs = SeqIO.to_dict(SeqIO.parse(open(Ref_Path), "fasta"))
if len(ReferenceSeqs) == 0:
    raise ValueError("Reference file is empty!")

# Get Focal sequences
TargetSeqs = SeqIO.to_dict(SeqIO.parse(open(Target_Path), "fasta"))
SeqIDs = TargetSeqs.keys()
if len(TargetSeqs) == 0:
    raise ValueError("Focal sequence file is empty!")


Output.write('\t'.join(["ID", "Length", "NbSub", "Weak2Weak", "Weak2Strong", "Strong2Weak", "Strong2Strong"]) + "\n")
for ID in SeqIDs:
    ref_seq = str(ReferenceSeqs[ID].seq)
    target_seq = str(TargetSeqs[ID].seq)

    substitutions = get_sub_number(ref_seq, target_seq)
    substitution_number = list(substitutions.values())
    Nsub = str(sum(substitution_number))

    Output.write('\t'.join([ID, str(len(ref_seq)), Nsub, '\t'.join(map(str, substitution_number))]) + "\n")


Output.close()

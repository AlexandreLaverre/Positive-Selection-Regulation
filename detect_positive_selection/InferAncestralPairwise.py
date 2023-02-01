#!/usr/bin/env python
# coding=utf-8

from Bio.Seq import Seq
from Bio import AlignIO
import os
import sys
from collections import defaultdict

sp = sys.argv[1]    # i.e: human or mouse
sample = sys.argv[2]    # i.e: CEBPA or HNF4A
BED_file = sys.argv[3]  # i.e: path + "Tools/JialinTool/data/" + sp + "/" + sp + "_ChIP-Seq/" + prefix + sample + ".bed2"
method = sys.argv[4]    # i.e: parsimony

####################################################################################
# Define variables and path
sister_sp = ["gorGor5", "panTro6"] if sp == "human" else ["gorGor5", "panTro6"]
prefix = "hsap_" if sp == "human" else "b6-"
genome = "hg38" if sp == "human" else "mm10"

path = "/Users/alaverre/Documents/Detecting_positive_selection/"
pathAlign = path + "data/genome_alignments/"
pathResults = path + "results/" + sp + "/" + sample + "/Alignments/"

PositiveSeq = path + "results/" + sp + "/" + sample + "/Model/posSet.fa"
Tree_file = pathAlign + sp + '_tree.nk'
Tree = open(Tree_file).read()

os.makedirs(pathResults + "MFA", exist_ok=True)
os.makedirs(pathResults + "running", exist_ok=True)
os.makedirs(pathResults + "focal_sequences", exist_ok=True)
os.makedirs(pathResults + "ancestral_sequences", exist_ok=True)
error = open(pathResults + "error_log.txt", "w+")

####################################################################################
# Extract regions from the genomes alignment
for sis in sister_sp:
    pathMAFs = pathAlign + sample + "_" + genome + "To" + sis + "/"
    GenomeAlignment = pathAlign + genome + "." + sis + ".maf"
    if not os.path.exists(pathMAFs):
        os.system("mafsInRegion " + BED_file + " -outDir " + pathMAFs + " " + GenomeAlignment)

# Extract sequences and estimate ancestral state by region
with open(BED_file, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        ID = i[3]
        PositiveSeqID = i[0] + "_" + str(int(i[1])+1) + "_" + i[2] + "_*"
        print(ID)

        PathRunning = pathResults + '/running/' + ID + "/"
        os.makedirs(PathRunning, exist_ok=True)

        # Extract reference sequence from positive Set
        print('seqkit grep -r -p "' + PositiveSeqID + '" ' + PositiveSeq + " > " + PathRunning + genome + ".fa")
        os.system('seqkit grep -r -p "' + PositiveSeqID + '" ' + PositiveSeq + " > " + PathRunning + genome + ".fa")
        os.system("sed -i '' -e 's/>.*/>" + genome + "/g' " + PathRunning + genome + ".fa")

        # Extract aligned sequences for each species
        for sis in sister_sp:
            MAF = pathAlign + sample + "_" + genome + "To" + sis + "/" + ID + ".maf"
            seq = []
            for multiple_alignment in AlignIO.parse(MAF, "maf"):
                for seqrec in multiple_alignment:
                    if seqrec.id.startswith(sis):
                        seqrec.seq = seqrec.seq.replace("-", "")
                        seq.append(seqrec.seq)

            out = PathRunning + sis + ".fa"
            output = open(out, "w+")
            output.write(">" + sis + '\n')
            output.write(str(Seq('').join(seq)) + '\n')
            output.close()

        # Make alignment with PECAN
        Align_pecan = PathRunning + "pecan.mfa"
        if not os.path.exists(Align_pecan):
            os.chdir(PathRunning)
            # Check if a sequence is found in each species
            if os.path.exists(PathRunning + "panTro6.fa") & os.path.exists(PathRunning + "gorGor5.fa"):
                os.system("java bp.pecan.Pecan -E '" + Tree +
                          "' -F " + '.fa '.join(sister_sp) + ".fa " + genome + ".fa -G pecan.mfa")
            else:
                error.write(ID + '\tMissing_Alignment\n')

        # Infer ancestral sequences with TreeTime
        Complete_alignment = pathResults + "MFA/" + ID + ".mfa"
        Focal_seq = pathResults + "focal_sequences/" + ID
        Ancestral_seq = pathResults + "ancestral_sequences/" + ID
        os.chdir(pathResults)

        if os.path.exists(Align_pecan):
            os.system("treetime ancestral --keep-overhangs --verbose 0 --aln " + Align_pecan
                      + " --tree " + Tree_file + " --method-anc " + method + " --outdir " + PathRunning)

            if os.path.exists(PathRunning + "/ancestral_sequences.fasta"):
                os.system("mv " + PathRunning + "/ancestral_sequences.fasta " + Complete_alignment)

                # Extract ancestral and focal sequences from alignment
                os.system("seqkit grep -p " + genome + " -p NODE_0000001 " + Complete_alignment + " > " + PathRunning + "anc_foc.mfa")

                # Remove gaps
                os.system("trimal -in " + PathRunning + "anc_foc.mfa -out " + PathRunning + "anc_foc_nogap.mfa" + " -nogaps -keepheader")

                # Extract fasta
                os.system("seqkit grep -p " + genome + " " + PathRunning + "anc_foc.mfa > " + Focal_seq + ".fa")
                os.system("seqkit grep -p " + genome + " " + PathRunning + "anc_foc_nogap.mfa > " + Focal_seq + "_nogap.fa")
                os.system("sed -i '' -e 's/" + genome + "/" + ID + "/g' " + Focal_seq + ".fa")
                os.system("sed -i '' -e 's/" + genome + "/" + ID + "/g' " + Focal_seq + "_nogap.fa")

                os.system("seqkit grep -p NODE_0000001 " + PathRunning + "anc_foc.mfa > " + Ancestral_seq + ".fa")
                os.system("seqkit grep -p NODE_0000001 " + PathRunning + "anc_foc_nogap.mfa > " + Ancestral_seq + "_nogap.fa")
                os.system("sed -i '' -e 's/NODE_0000001/" + ID + "/g' " + Ancestral_seq + ".fa")
                os.system("sed -i '' -e 's/NODE_0000001/" + ID + "/g' " + Ancestral_seq + "_nogap.fa")

                #os.system("rm -r " + PathRunning)
            else:
                error.write(ID + '\tMissing_Ancestral\n')

error.close()

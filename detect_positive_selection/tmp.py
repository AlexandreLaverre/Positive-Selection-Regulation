#!/usr/bin/env python
# coding=utf-8

from Bio.Seq import Seq
from Bio import AlignIO
import os
import sys
from collections import defaultdict

sp = "human"
sample = "CEBPA"

####################################################################################
# Define variables and path
sister_sp = ["gorGor5", "panTro6"] if sp == "human" else ["gorGor5", "panTro6"]
prefix = "hsap_" if sp == "human" else "b6-"
genome = "hg38" if sp == "human" else "mm10"

path = "/Users/alaverre/Documents/Detecting_positive_selection/"
pathAlign = path + "data/genome_alignments/"
pathResults = path + "results/" + sp + "/" + sample + "/Alignments-test/"

BED_file = path + "Tools/JialinTool/data/" + sp + "/" + sp + "_ChIP-Seq/" + prefix + sample + "_based0.bed"

PositiveSeq = path + "results/" + sp + "/" + sample + "/Model/posSet.fa"
Tree_file = pathAlign + sp + '_tree.nk'
Tree = open(Tree_file).read()

os.makedirs(pathResults + "MFA", exist_ok=True)
os.makedirs(pathResults + "running", exist_ok=True)
os.makedirs(pathResults + "focal_sequences", exist_ok=True)
os.makedirs(pathResults + "ancestral_sequences", exist_ok=True)
error = open(pathResults + "error_log.txt", "w+")

####################################################################################

with open(BED_file, 'r') as f1:
    for i in f1.readlines()[1:10]:
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

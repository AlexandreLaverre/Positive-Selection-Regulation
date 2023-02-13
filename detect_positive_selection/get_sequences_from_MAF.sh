#!/bin/bash

sp=$1    		# i.e: human or mouse
sample=$2    	# i.e: CEBPA or HNF4A
BED_file=$3	# i.e: 

path="/Users/alaverre/Documents/Detecting_positive_selection/"
pathGenomeAlign="${path}/data/genome_alignments/${sp}/Dog_triplet_ancestor.maf"
pathResults="${path}/results/positive_selection/${sp}/${sample}/"
pathAlign="${pathResults}/Alignments/"

mkdir -p "${pathAlign}"
mkdir -p "${pathResults}/focal_sequences/"
mkdir -p "${pathResults}/ancestral_sequences/"

if [ ${sp} = "dog" ]; then
    sp_name="Canis_lupus_familiaris"
    anc_name="Anc09"
fi

source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
conda activate phyml

#########################################################################

# Retrieve alignments of regions from whole genome alignment
mafsInRegion ${BED_file} -outDir ${pathAlign}/ ${pathGenomeAlign}

for align in `ls ${pathAlign}/*maf`
do
	ID="$(basename "${align}" .maf)"
	
	# Convert to FASTA
	msa_view -i MAF -o FASTA ${align} > ${pathAlign}/${ID}.mfa
	sed -i '' 's/ //g' ${pathAlign}/${ID}.mfa
	
	# Keep only focal and ancestral sequences
	seqkit grep -p ${sp_name} -p ${anc_name} ${pathAlign}/${ID}.mfa > ${pathAlign}/${ID}_anc_foc.mfa
	
	# Remove GAP
	trimal -nogaps -keepheader -in ${pathAlign}/${ID}_anc_foc.mfa  ${pathAlign}/${ID}_anc_foc_nogap.mfa
	
	# Get focal and ancestral sequences separatly
	seqkit grep -p ${sp_name} ${pathAlign}/${ID}_anc_foc_nogap.mfa > ${pathResults}/focal_sequences/${ID}.fa
	seqkit grep -p ${anc_name} ${pathAlign}/${ID}_anc_foc_nogap.mfa > ${pathResults}/ancestral_sequences/${ID}.fa
	
done
	
#########################################################################
	
	
	

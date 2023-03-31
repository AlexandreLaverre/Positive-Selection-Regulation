#!/bin/bash

sp=$1    		# i.e: human or mouse
sample=$2    	# i.e: CEBPA or HNF4A
TF=$3
BED_file=$4
cluster=$5

if [ ${cluster} = "local" ]; then
	export path=/Users/alaverre/Documents/Detecting_positive_selection
	source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
else
	export path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/
	source /users/alaverre/Tools/mambaforge/etc/profile.d/conda.sh
fi

conda activate MAF

pathGenomeAlign="${path}/data/genome_alignments/${sp}/triplet_ancestor.maf"
pathResults="${path}/results/positive_selection/${sp}/${sample}/${TF}/Alignments/"
pathAlign="${pathResults}/MAFs/"

mkdir -p "${pathAlign}"
mkdir -p "${pathResults}/focal_sequences/"
mkdir -p "${pathResults}/ancestral_sequences/"

if [ ${sp} = "dog" ]; then
    sp_name="Canis_lupus_familiaris"
    anc_name="Anc09"
fi

if [ ${sp} = "human" ]; then
    sp_name="Homo_sapiens"
    anc_name="fullTreeAnc105"
fi

if [ ${sp} = "mouse" ]; then
    sp_name="Mus_musculus"
    anc_name="fullTreeAnc35"
fi

if [ ${sp} = "rat" ]; then
    sp_name="Rattus_norvegicus"
    anc_name="fullTreeAnc38"
fi

if [ ${sp} = "macaca" ]; then
    sp_name="Macaca_mulatta"
    anc_name="fullTreeAnc91"
fi

if [ ${sp} = "chicken" ]; then
    sp_name="Gallus_gallus"
    anc_name=""
fi

#########################################################################

# Retrieve alignments of regions from whole genome alignment
mafsInRegion ${BED_file} -outDir ${pathAlign}/ ${pathGenomeAlign}

for ID in `cat "${BED_file}" | cut -f 4`
do
	#ID=$(echo ${line} | cut -f 4 -d ' ')
	echo "${ID}"
	align="${pathAlign}"/"${ID}".maf
  echo "$align"

	# Convert to FASTA
	msa_view --missing-as-indels -i MAF -o FASTA "${align}" > "${pathAlign}"/"${ID}".mfa
	sed -i 's/ //g' "${pathAlign}"/"${ID}".mfa

  # Check if an alignment exists with the other species
	if [ -s "${pathAlign}"/"${ID}".mfa ]; then
	
	    # Keep only focal and ancestral sequences
	    seqkit grep -p ${sp_name} -p "${anc_name}" "${pathAlign}"/"${ID}".mfa > "${pathAlign}"/"${ID}"_anc_foc.mfa

	    # Remove GAP
	    trimal -nogaps -keepheader -in "${pathAlign}"/"${ID}"_anc_foc.mfa -out "${pathAlign}"/"${ID}"_anc_foc_nogap.mfa
	    rm ${pathAlign}/"${ID}"_anc_foc.mfa "${pathAlign}"/"${ID}".mfa

	    # Check if alignment still exists after removing GAPs
	    if [ -f ${pathAlign}/"${ID}"_anc_foc_nogap.mfa ]; then
		    # Get focal and ancestral sequences separatly
		    seqkit grep -p ${sp_name} ${pathAlign}/"${ID}"_anc_foc_nogap.mfa > ${pathResults}/focal_sequences/"${ID}"_nogap.fa
		    seqkit grep -p ${anc_name} ${pathAlign}/"${ID}"_anc_foc_nogap.mfa > ${pathResults}/ancestral_sequences/"${ID}"_nogap.fa

		    # Change species names to ID
		    sed -i "s/${sp_name}/"${ID}"/g" ${pathResults}/focal_sequences/"${ID}"_nogap.fa
		    sed -i "s/${anc_name}/"${ID}"/g" ${pathResults}/ancestral_sequences/"${ID}"_nogap.fa
	    else
	    	echo "${ID}" >> ${pathResults}/missing.txt
	    fi
	    
  	else
    		echo "${ID}" >> ${pathResults}/missing.txt
  	fi
  
done
	
#########################################################################
	
	
	

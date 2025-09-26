#!/bin/bash

sp=$1    		  # e.g: human
pathResults=$2
BED_file=$3
AncNode="${4:-ancestral}"

#########################################################################
# Define paths
export path=${path:-"$(pwd)/../../"}

pathGenomeAlign="${path}/data/genome_alignments/${sp}/triplet_${AncNode}.maf.gz"
pathAlign="${pathResults}/MAFs/"

if [ ! -f "${pathGenomeAlign}" ]; then
  echo "Genome alignment file not found!"
  exit
fi

mkdir -p "${pathAlign}"
echo "Results will be stored in ${pathResults}"
mkdir -p "${pathResults}/focal_sequences/"
mkdir -p "${pathResults}/sister_sequences/"
mkdir -p "${pathResults}/ancestral_sequences/"
# Define focal and ancestral names

source ${path}/scripts/config/params.sh "$sp"
sister_name="$(echo "$close_species" | cut -d "," -f 2)"


if [[ $AncNode != 'ancestral' ]]; then
  anc_name="fullTree${AncNode}"
fi

echo "Focal species: $sp_name"
echo "Ancestral species: $anc_name"
#########################################################################
# Retrieve alignments of regions from whole genome alignment
mafsInRegion "$BED_file" -outDir "${pathAlign}"/ "$pathGenomeAlign"

# Retrieve ancestral and focal sequences
for ID in `cat "$BED_file" | cut -f 4`
do
	align="${pathAlign}"/"${ID}".maf

	# Convert to FASTA
	msa_view --missing-as-indels -i MAF -o FASTA "$align" > "${pathAlign}"/"${ID}".mfa
	sed -i 's/ //g' "${pathAlign}"/"${ID}".mfa

  # Check if an alignment exists with the other species
	if [ -s "${pathAlign}"/"${ID}".mfa ]; then
	
	    # Keep only focal and ancestral sequences
	    seqkit grep -p "$sp_name" -p "$anc_name" "${pathAlign}"/"${ID}".mfa > "${pathAlign}"/"${ID}"_anc_foc.mfa

	    # Remove GAP
	    trimal -nogaps -keepheader -in "${pathAlign}"/"${ID}"_anc_foc.mfa -out "${pathAlign}"/"${ID}"_anc_foc_nogap.mfa
	    rm "${pathAlign}"/"${ID}".mfa "${pathAlign}"/"${ID}"_anc_foc.mfa

	    # Check if alignment still exists after removing GAPs
	    if [ -f "${pathAlign}"/"${ID}"_anc_foc_nogap.mfa ]; then
		    # Get focal and ancestral sequences separately
		    seqkit grep -p "$sp_name" "${pathAlign}"/"${ID}"_anc_foc_nogap.mfa > "${pathResults}"/focal_sequences/"${ID}"_nogap.fa
        seqkit grep -p "$anc_name" "${pathAlign}"/"${ID}"_anc_foc_nogap.mfa > "${pathResults}"/ancestral_sequences/"${ID}"_nogap.fa
        #seqkit grep -p "$sister_name" "${pathAlign}"/"${ID}"_anc_foc_nogap.mfa > "${pathResults}"/sister_sequences/"${ID}"_nogap.fa

		    # Change species names to ID
		    sed -i "s/${sp_name}/${ID}/g" "${pathResults}"/focal_sequences/"${ID}"_nogap.fa
		    sed -i "s/${anc_name}/${ID}/g" "${pathResults}"/ancestral_sequences/"${ID}"_nogap.fa
		    #sed -i "s/${sister_name}/${ID}/g" "${pathResults}"/sister_sequences/"${ID}"_nogap.fa
	    else
	    	echo "${ID}" >> "${pathResults}"/missing.txt
	    fi
	    
  	else
    		echo "${ID}" >> "${pathResults}"/missing.txt
  	fi
  
done

#########################################################################
	
	
	

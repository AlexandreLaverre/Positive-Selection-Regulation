#!/bin/bash

export sp=$1

####################################################################################

export path=/Users/alaverre/Documents/Detecting_positive_selection/
export pathResults=${path}/results/substitution_matrix/${sp}/
export pathAlignment=${path}/data/genome_alignments/${sp}/
export SpeciesTree=${path}/data/species_trees/${sp}_tree.nk

if [ ${sp} = "dog" ]; then
	export Alignment=${pathAlignment}/Dog_triplets.maf
	chroms=({1..38} "MT" "X")
fi

if [ ${sp} = "human" ]; then
    export Alignment=${pathAlignment}/hg38.gorGor5.maf
    chroms=(chr{1..22} "M" "X" "Y")
fi

source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
conda activate phyml

#########################################################################
mkdir -p ${pathAlignment}/per_chrom/PHYLIPs/
mkdir -p ${pathResults}/phyML_logs/
	
# Get alignment per chromosome
if [ -e ${pathAlignment}/per_chrom/MAFs/ ]; then
	echo "MAF split already done!"
else
	mkdir -p ${pathAlignment}/per_chrom/MAFs/scaffolds/	
	mafSplit -byTarget -useFullSequenceName _.bed ${pathAlignment}/per_chrom/MAFs/scaffolds/ ${Alignment}
fi

# Run parameters estimation with phyML for each chromosome
for chr in ${chroms[@]}
do
	echo "########### ${chr} ###########"
	
	if [ -e ${pathResults}/${chr}.txt ]; then
		echo "Substitution matrix already done!"
	else

		# Change MAF to PHYLIP format
		phylip_file="${pathAlignment}/per_chrom/PHYLIPs/${chr}.phylip"
		
		if [ -e ${phylip_file} ]; then
			echo "MAF to PHYLIP already done!"
		else
			mv ${pathAlignment}/per_chrom/MAFs/scaffolds/${chr}.maf ${pathAlignment}/per_chrom/MAFs/${chr}.maf # just for cleaning (running only on non-scaffolds)
			msa_view -i MAF -o PHYLIP ${pathAlignment}/per_chrom/MAFs/${chr}.maf > ${phylip_file}
		fi

		# Run phyML
		phyml --datatype nt --sequential --quiet --model GTR -o r --input ${phylip_file} ${SpeciesTree}
		
		# Make clean substitutions matrix
		grep -A5 "rate matrix" ${phylip_file}_phyml_stats.txt | tail -n5 > ${pathResults}/${chr}.txt
		sed -i '' "s/\[A---------C---------G---------T------\]/A C G T/g" ${pathResults}/${chr}.txt
		sed -i '' "s/  //g" ${pathResults}/${chr}.txt
		sed -i '' "s/-/ -/g" ${pathResults}/${chr}.txt
	fi
	
done

# Move phyML output files
rm ${pathAlignment}/per_chrom/PHYLIPs/*phyml_tree.txt 
mv ${pathAlignment}/per_chrom/PHYLIPs/*phyml_stats.txt ${pathResults}/phyML_logs/

echo "Done! =)"

#########################################################################


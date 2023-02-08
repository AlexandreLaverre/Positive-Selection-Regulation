#!/bin/bash

export sp=$1

####################################################################################

export path=/Users/alaverre/Documents/Detecting_positive_selection/
export pathResults=${path}/results/substitution_matrix/${sp}/
export pathAlignment=${path}/data/genome_alignments/${sp}/
export SpeciesTree=${path}/data/species_trees/${sp}_tree.nk

if [ ${sp} = "dog" ]; then
	export Alignment=${pathAlignment}/Dog_triplets.maf
	export chroms=({1..38} "MT" "X")
fi

if [ ${sp} = "human" ]; then
    export Alignment=${pathAlignment}/hg38.gorGor5.maf
    export chroms=(chr{1..22} "M" "X" "Y")
fi

conda activate phyml

#########################################################################

mkdir -p ${pathAlignment}/per_chrom/MAFs/scaffolds/
mkdir -p ${pathAlignment}/per_chrom/PHYLIPs/
mkdir -p ${pathResults}/phyML_logs/

# Get alignment per chromosome
if [ -e ${pathSites}/${alnprefix}_${chr}_4d-sites.ss ]; then
	echo "MAF split already done!"
else
	mafSplit -byTarget -useFullSequenceName _.bed ${pathResults}/per_chrom/MAFs/scaffolds/ ${Alignment}
fi

# Run parameters estimation with phyML for each chromosome
for chr in ${chroms[@]}
do
	echo "########### ${chr} ###########"
	phylip_file="${pathResults}/per_chrom/PHYLIPs/${chr}.phylip"

	# Change MAF to PHYLIP format
	mv ${pathResults}/per_chrom/MAFs/scaffolds/${chr}.maf ${pathResults}/per_chrom/MAFs/${chr}.maf # just for cleaning (running only on non-scaffolds)
	msa_view -i MAF -o PHYLIP ${pathResults}/per_chrom/MAFs/${chr}.maf > ${phylip_file}

	# Run phyML
	phyml --datatype nt --sequential --quiet --model GTR -o r --input ${phylip_file} ${SpeciesTree}
	
	# Make clean substitutions matrix
	grep -A5 "rate matrix" ${phylip_file}_phyml_stats.txt | tail -n5 > ${pathResults}/${chr}.txt
	sed -i'' -e "s/  //g" ${pathResults}/${chr}.txt
	sed -i'' -e "s/---/ /g" ${pathResults}/${chr}.txt
	sed -i'' -e "s/-/ -/g" ${pathResults}/${chrchrom}.txt
	
done

# Move phyML output files
rm ${pathAlignment}/per_chrom/PHYLIPs/*phyml_tree.txt 
mv ${pathAlignment}/per_chrom/PHYLIPs/*phyml_stats.txt ${pathResults}/phyML_logs

echo "Done! =)"

#########################################################################



#!/bin/bash

export sp=$1
export cluster=$2

####################################################################################

if [ ${cluster} = "local" ]; then
	export path=/Users/alaverre/Documents/Detecting_positive_selection
	source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
else
	export path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/
	source /users/alaverre/Tools/mambaforge/etc/profile.d/conda.sh
fi

export pathResults=${path}/results/substitution_matrix/${sp}/
export pathAlignment=${path}/data/genome_alignments/${sp}/
export pathGFF=${path}/data/genome_sequences/${sp}/

export SpeciesTree=${path}/data/species_trees/${sp}_tree.nk
export Alignment=${pathAlignment}/triplet_ancestor.maf.gz

####################################################################################
# Define names, chromosomes and files
source ${path}/scripts/config/params.sh ${sp} ${cluster}
conda activate MAF

#########################################################################
mkdir -p ${pathAlignment}/per_chrom/PHYLIPs/
mkdir -p ${pathResults}/phyML_logs/
	
# Get alignment per chromosome
if [ -e ${pathAlignment}/per_chrom/MAFs/ ]; then
	echo "MAF split already done!"
else
  echo "Splitting whole genome alignment by chromosome"
	mkdir -p ${pathAlignment}/per_chrom/MAFs/scaffolds/	
	mafSplit -byTarget -useFullSequenceName _.bed ${pathAlignment}/per_chrom/MAFs/scaffolds/ ${Alignment}
fi

# Run parameters estimation with phyML for each chromosome
for chr in "${chroms[@]}"
do
	echo "########### ${chr} ###########"
	
	if [ -e ${pathResults}/${chr}.txt ]; then
		echo "Substitution matrix already done!"
	else

		# Mask exons in MAF
		if [ -e ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf ]; then
			echo "MAF already masked for exons!"
		else
			# Get a GFF with sorted exons per chromosome
			if [ ! -e ${pathGFF}/exons.uniq.sorted.${spID}.gff ]; then
				mkdir -p ${pathGFF}/GFF_per_chrom
				zcat < ${pathGFF}/${spID}${GFF_suffix}  | grep -w "exon" > ${pathGFF}/exons_${spID}.gff
				cut -f 1-8 ${pathGFF}/exons_${spID}.gff | sort -u > ${pathGFF}/exons.uniq.${spID}.gff 
				sort -t $'\t' -k1,1V -k4,4h -k5,5h ${pathGFF}/exons.uniq.${spID}.gff > ${pathGFF}/exons.uniq.sorted.${spID}.gff
				rm ${pathGFF}/exons_${spID}.gff ${pathGFF}/exons.uniq.${spID}.gff
				
			fi

	      		echo "Masking for exons..."
	      		grep -w "^${chr}" ${pathGFF}/exons.uniq.sorted.${spID}.gff > ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${spID}.gff
	      		
	      		if [ ${sp} = "human" ] || [ ${sp} = "mouse" ] || [ ${sp} = "rat" ]; then
							# use UCSC chromosomes names
							sed -i'' 's/^/chr/g' ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${spID}.gff 
							mv ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${spID}.gff  ${pathGFF}/GFF_per_chrom/chr${chr}.exons.uniq.sorted_${spID}.gff 
							chr="chr${chr}"
						fi

	      		if [ ! -s ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${spID}.gff ]; then
	        			echo "GFF file for ${chr} is empty! Is-it normal?"
	        			cp ${pathAlignment}/per_chrom/MAFs/scaffolds/${chr}.maf ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf
	        	else
	        			maf_parse -o MAF --features ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${spID}.gff --mask-features ${close_species} ${pathAlignment}/per_chrom/MAFs/scaffolds/${chr}.maf > ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf
	        	fi
		fi
		
		# Change MAF to PHYLIP format
		phylip_file="${pathAlignment}/per_chrom/PHYLIPs/${chr}.exons_masked.phylip"
		
		if [ -e ${phylip_file} ]; then
			echo "MAF to PHYLIP already done!"
		else
			msa_view -i MAF -o PHYLIP ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf > ${phylip_file}
		fi

		# Run phyML
		phyml --datatype nt --sequential --quiet --model GTR -o r --input ${phylip_file} ${SpeciesTree}
		
		# Make clean substitutions matrix
		grep -A5 "rate matrix" ${phylip_file}_phyml_stats.txt | tail -n5 > ${pathResults}/${chr}.txt
		sed -i "s/\[A---------C---------G---------T------\]/A C G T/g" ${pathResults}/${chr}.txt
		sed -i "s/  //g" ${pathResults}/${chr}.txt
		sed -i "s/-/ -/g" ${pathResults}/${chr}.txt
	fi
	
done

# Move phyML output files
rm ${pathAlignment}/per_chrom/PHYLIPs/*phyml_tree.txt 
mv ${pathAlignment}/per_chrom/PHYLIPs/*phyml_stats.txt ${pathResults}/phyML_logs/

echo "Done! =)"

#########################################################################


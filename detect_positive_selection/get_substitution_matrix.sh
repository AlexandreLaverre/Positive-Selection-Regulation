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
export SpeciesTree=${path}/data/species_trees/${sp}_tree.nk
export pathGFF=${path}/data/genome_sequences/${sp}/

export Alignment=${pathAlignment}/triplet_ancestor.maf

if [ ${sp} = "dog" ]; then
	chroms=({1..38} "MT" "X")
	prefix="Canis_lupus_familiaris.CanFam3.1.104"
	species="Canis_lupus_familiaris,Canis_lupus_lupus,Lycaon_pictus"
fi

if [ ${sp} = "human" ]; then
	chroms=(chr{1..22} "chrM" "chrX" "chrY")
	prefix="Homo_sapiens.GRCh38.104"
	species="Homo_sapiens,Pan_troglodytes,Gorilla_gorilla"
fi

if [ ${sp} = "mouse" ]; then
	chroms=(chr{1..19} "chrM" "chrX" "chrY")
	prefix="Mus_musculus.GRCm39.104"
	species="Mus_musculus,Mus_spretus,Mus_caroli"
fi

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
			if [ ! -e ${pathGFF}/exons.uniq.sorted.${prefix}.gff ]; then
			  mkdir -p ${pathGFF}/GFF_per_chrom
				zcat < ${pathGFF}/${prefix}.gff3.gz | grep -w "exon" > ${pathGFF}/exons_${prefix}.gff
				cut -f 1-8 ${pathGFF}/exons_${prefix}.gff | sort -u > ${pathGFF}/exons.uniq.${prefix}.gff 
				sort -k1,1V -k4,4h -k5,5rh -k3,3r ${pathGFF}/exons.uniq.${prefix}.gff > ${pathGFF}/exons.uniq.sorted.${prefix}.gff
				rm ${pathGFF}/exons_${prefix}.gff ${pathGFF}/exons.uniq.${prefix}.gff

				if [ ${sp} != "dog" ] ; then
					sed -i 's/^/chr/g' ${pathGFF}/exons.uniq.sorted.${prefix}.gff
					sed -i 's/chrMT/chrM/g' ${pathGFF}/exons.uniq.sorted.${prefix}.gff
				fi
			fi

      		echo "Masking for exons..."
      		grep -w "^${chr}" ${pathGFF}/exons.uniq.sorted.${prefix}.gff > ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff

      		if [ ! -s ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff ]; then
        			echo "Weird! GFF file for ${chr} is empty!"
        		fi

		 maf_parse -o MAF --features ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff --mask-features ${species} ${pathAlignment}/per_chrom/MAFs/scaffolds/${chr}.maf > ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf

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


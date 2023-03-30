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
export Alignment=${pathAlignment}/triplet_ancestor.maf

####################################################################################

if [ ${sp} = "dog" ]; then
	chroms=({1..38} "MT" "X")
	prefix="Canis_lupus_familiaris.CanFam3.1.104"
	species="Canis_lupus_familiaris,Canis_lupus_lupus,Lycaon_pictus"
fi

if [ ${sp} = "human" ]; then
	chroms=({1..22} "X" "Y")
	prefix="Homo_sapiens.GRCh38.104"
	species="Homo_sapiens,Pan_troglodytes,Gorilla_gorilla"
fi

if [ ${sp} = "mouse" ]; then
	chroms=({1..19} "X" "Y")
	prefix="Mus_musculus.GRCm39.104"
	species="Mus_musculus,Mus_spretus,Mus_caroli"
fi

if [ ${sp} = "macaca" ]; then
	chroms=(CM002977.3 CM002980.3 CM002984.2 CM002983.2 CM002981.2 CM002982.3 CM002991.3 CM002985.3 CM002987.3 CM002992.3 CM002989.3 CM002979.2 CM002978.2 CM002988.2 CM002986.2 CM002994.2 CM002990.2 CM002995.2 CM002996.3 CM002993.2 CM002997.3 CM003438.1) # from NCBI (1 to 20, X and Y)
	prefix="sup2kb_GCA_000772875.3_Mmul_8.0.1_genomic"
	species="Macaca_mulatta,Macaca_fascicularis,Macaca_nemestrina"
fi

if [ ${sp} = "chicken" ]; then
	chroms=({1..33} "Z" "W" "MT")
	prefix="Gallus_gallus.GRCg6a.104"
	species="Gallus_gallus,Coturnix_japonica,Meleagris_gallopavo"
fi


if [ ${sp} = "rat" ]; then
	chroms=({1..20} "X" "Y")
	prefix="Rattus_norvegicus.Rnor_6.0.104"
	species="Rattus_norvegicus,Mus_musculus,Acomys_cahirinus"
fi

if [ ${sp} = "cat" ]; then
	chroms=(CM0013{78..96}.2) # from NCBI (A1 to F2 and X)
	prefix="GCA_000181335.3_Felis_catus_8.0_genomic"
	species="Felis_catus,Felis_nigripes,Puma_concolor"
fi

if [ ${sp} = "cattle" ]; then
	chroms=(CM000177.6 CM000178.7 CM000179.6 CM000180.6 CM000181.7 CM000182.7 CM000183.7 CM000184.6 CM000185.7 CM000186.6 CM000187.6 CM000188.7 CM000189.7 CM000190.6 CM000191.7 CM000192.5 CM000193.6 CM000194.7 CM000195.6 CM000196.6 CM000197.6 CM000198.6 CM000199.8 CM000200.7 CM000201.6 CM000202.7 CM000203.6 CM000204.6 CM000205.6 CM000206.6 CM001061.2) # from NCBI (1 to 29, X and Y)
	prefix="GCA_000003205.7_Btau_5.0.1_genomic"
	species="Bos_taurus,Bos_indicus,Bos_mutus"
fi

if [ ${sp} = "rabbit" ]; then
	chroms=(CM000{790..811}.1) # from NCBI (chr1 to 21 and X)
	prefix="GCA_000003625.3_OryCun2.0_genomic"
	species="Oryctolagus_cuniculus,Lepus_americanus,Ochotona_princeps"
fi

if [ ${sp} = "macaca" ] || [ ${sp} = "cat" ] || [ ${sp} = "cattle" ] || [ ${sp} = "rabbit" ]; then
	# Use GFF from NCBI
	suffix=".gff.gz"
else
	suffix=".gff3.gz"
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
				zcat < ${pathGFF}/${prefix}${suffix}  | grep -w "exon" > ${pathGFF}/exons_${prefix}.gff
				cut -f 1-8 ${pathGFF}/exons_${prefix}.gff | sort -u > ${pathGFF}/exons.uniq.${prefix}.gff 
				sort -t $'\t' -k1,1V -k4,4h -k5,5h ${pathGFF}/exons.uniq.${prefix}.gff > ${pathGFF}/exons.uniq.sorted.${prefix}.gff
				rm ${pathGFF}/exons_${prefix}.gff ${pathGFF}/exons.uniq.${prefix}.gff
				
			fi

	      		echo "Masking for exons..."
	      		grep -w "^${chr}" ${pathGFF}/exons.uniq.sorted.${prefix}.gff > ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff
	      		
	      		if [ ${sp} = "human" ] || [ ${sp} = "mouse" ] || [ ${sp} = "rat" ]; then
				# use UCSC chromosomes names
				sed -i'' 's/^/chr/g' ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff 
				mv ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff  ${pathGFF}/GFF_per_chrom/chr${chr}.exons.uniq.sorted_${prefix}.gff 
				chr="chr${chr}"
			fi

	      		if [ ! -s ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff ]; then
	        			echo "GFF file for ${chr} is empty! Is-it normal?"
	        			cp ${pathAlignment}/per_chrom/MAFs/scaffolds/${chr}.maf ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf
	        		else
	        			maf_parse -o MAF --features ${pathGFF}/GFF_per_chrom/${chr}.exons.uniq.sorted_${prefix}.gff --mask-features ${species} ${pathAlignment}/per_chrom/MAFs/scaffolds/${chr}.maf > ${pathAlignment}/per_chrom/MAFs/${chr}.exons_masked.maf

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


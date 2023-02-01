#!/bin/bash

export sp=$1		# i.e: dog
export sample=$2	# i.e: FOXA1 CEBPA All
export container=$3	# i.e: docker or conda (only conda working currently)
export resume=$4	# i.e: -resume or nothing

export path=/Users/alaverre/Documents/Detecting_positive_selection
export pathResults=${path}/results/${sp}/peaks_calling
export pathData=${path}/data

export release=104

source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
conda activate nextflow

#########################################################################

if [ ${sp} = "human" ]; then
	export spID="Homo_sapiens.GRCh38"
	export genomesize=2701262066 # from 50: https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/conf/igenomes.config
fi

if [ ${sp} == "dog" ]; then
	export spID="Canis_lupus_familiaris.CanFam3.1"
	export genomesize=2237684358
fi

export sampleID=${pathData}/ChIP-seq/${sp}/${sample}_samples_input.csv
export genome=${pathData}/genome_sequence/${sp}/${spID}.dna_sm.primary_assembly.fa
export GTF=${pathData}/genome_sequence/${sp}/${spID}.${release}.gtf
export blacklist=${pathData}/genome_sequence/${sp}/black_list.txt

if [ -f "${pathResults}/indexes/${spID}.dna_sm.primary_assembly.rev.2.bt2" ]; then
    echo "Indexes already done!"
    export index="--bwa_index ${pathResults}/indexes/ --bowtie2_index ${pathResults}/indexes/"
else
    export index="--save_reference"
fi


#########################################################################
mkdir -p ${pathResults}

nextflow run nf-core/chipseq --input ${sampleID} --outdir ${pathResults}/${sample} --fasta ${genome} --gtf ${GTF} --blacklist ${blacklist} --aligner bowtie2 --macs_gsize ${genomesize} -profile ${container} -with-conda true ${index} --max_memory '32.GB' --max_cpus 8 --skip_igv ${resume}

#########################################################################


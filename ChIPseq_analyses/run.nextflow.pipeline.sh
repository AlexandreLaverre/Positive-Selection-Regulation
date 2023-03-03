#!/bin/bash

export sp=$1		# i.e: dog
export sample=$2	# i.e: FOXA1 CEBPA All
export container=$3	# i.e: docker or conda (only conda working currently)
export threads=$4	# i.e: number of threads to use
export cluster=$5	# i.e: local or cluster
export resume=$6	# i.e: -resume or nothing

if [ ${cluster} = "local" ]; then
	export path=/Users/alaverre/Documents/Detecting_positive_selection
	export pathConda="/Users/alaverre/miniconda3/etc/profile.d/conda.sh"
else
	export path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/
	export pathConda="/users/alaverre/Tools/mambaforge/etc/profile.d/conda.sh"
fi

export pathResults=${path}/results/peaks_calling/${sp}/
export pathData=${path}/data
export pathScripts=${path}/scripts/ChIPseq_analyses

mkdir -p ${pathResults}

export release=104

#########################################################################

if [ ${sp} = "human" ]; then
	export spID="Homo_sapiens.GRCh38"
	export genomesize=2701262066 # from 50: https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/conf/igenomes.config
	export blacklist="--blacklist  ${pathData}/genome_sequences/${sp}/blacklist.txt"
fi

if [ ${sp} == "mouse" ]; then
	export spID="Mus_musculus.GRCm39"
	export genomesize=2307679482
	export blacklist="--blacklist  ${pathData}/genome_sequences/${sp}/blacklist.txt"
fi

if [ ${sp} == "dog" ]; then
	export spID="Canis_lupus_familiaris.CanFam3.1"
	export genomesize=2237684358
	export blacklist=""
fi

export sampleID=${pathData}/ChIP-seq/${sp}/${sample}_samples_input.csv
export genome=${pathData}/genome_sequences/${sp}/${spID}.dna_sm.primary_assembly.fa.gz
export GTF=${pathData}/genome_sequences/${sp}/${spID}.${release}.gtf.gz

if [ -f "${pathResults}/indexes/${spID}.dna_sm.toplevel.rev.2.bt2" ]; then
    echo "Indexes already done!"
    export index="--bwa_index ${pathResults}/indexes/ --bowtie2_index ${pathResults}/indexes/"
else
    export index="--save_reference"
fi

#########################################################################
echo "#!/bin/bash" > ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}

if [ ${cluster} = "cluster" ]; then
	echo "#SBATCH --job-name="ChIP-seq_peaks_calling_${sp}_${sample}" >>  ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --output=${pathScripts}/std_output_peaks_calling_${sp}_${sample}}.txt" >>  ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --error=${pathScripts}/std_error_peaks_calling_${sp}_${sample}.txt" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --partition=cpu" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --mem=30G" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --time=48:00:00" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
fi

echo "source ${pathConda}" > ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
echo "conda activate nextflow" > ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}

#########################################################################

echo "nextflow run nf-core/chipseq --input ${sampleID} --outdir ${pathResults}/${sample} --fasta ${genome} --gtf ${GTF} ${blacklist} --aligner bowtie2 --macs_gsize ${genomesize} -profile ${container} -with-conda true ${index} --max_memory '32.GB' --max_cpus ${threads} ${resume}" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}

#########################################################################


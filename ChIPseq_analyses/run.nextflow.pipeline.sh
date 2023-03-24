#!/bin/bash

export sp=$1		# i.e: dog
export sample=$2	# i.e: Wilson Schmidt
export source=$3	# i.e: Ensembl (for dog and cat) and NCBI (for all others)
export container=$4	# i.e: docker singularity conda
export threads=$5	# i.e: number of threads to use
export cluster=$6	# i.e: local or cluster
export resume=$7	# i.e: -resume or nothing

if [ ${cluster} = "local" ]; then
	export path=/Users/alaverre/Documents/Detecting_positive_selection
	export pathConda="/Users/alaverre/miniconda3/etc/profile.d/conda.sh"
else
	export path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/
	export pathConda="/users/alaverre/Tools/mambaforge/etc/profile.d/conda.sh"
fi

export pathResults=${path}/results/peaks_calling/${sp}/
export pathData=${path}/data
export pathScripts=${path}/scripts/ChIPseq_analyses/logs

mkdir -p ${pathResults}

#########################################################################
# Define parameters according to species and data source

if [ ${sp} = "human" ]; then
	export spID="Homo_sapiens.GRCh38"
	export genomesize=2701262066 # from 50: https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/conf/igenomes.config
	export blacklist="--blacklist  ${pathData}/genome_sequences/${sp}/blacklist.txt"
fi

if [ ${sp} == "mouse" ]; then
	export spID="Mus_musculus.GRCm38"
	export genomesize=2307679482
	export blacklist="--blacklist  ${pathData}/genome_sequences/${sp}/blacklist.txt"
fi

if [ ${sp} == "dog" ]; then
	export spID="Canis_lupus_familiaris.CanFam3.1"
	export genomesize=2237684358
	export blacklist=""
fi

if [ ${sp} == "cat" ]; then
	export spID="Felis_catus.Felis_catus_9.0"
	export genomesize=235000000
	export blacklist=""
fi

if [ ${sp} == "rat" ]; then
	export spID="GCF_000001895.5_Rnor_6.0"
	export genomesize=2375372135
	export blacklist=""
fi

if [ ${sp} == "macaca" ]; then
	export spID="sup2kb_GCA_000772875.3_Mmul_8.0.1"
	export genomesize=2498932238
	export blacklist=""
fi

if [ ${sp} == "cattle" ]; then
	export spID="GCF_000003205.7_Btau_5.0.1"
	export genomesize=2370644326
	export blacklist=""
fi

if [ ${sp} == "pig" ]; then
	export spID="GCF_000003025.5_Sscrofa10.2"
	export genomesize=2105185708
	export blacklist=""
fi

if [ ${sp} == "chicken" ]; then
	export spID="Gallus_gallus.GRCg6a"
	export genomesize=974987959
	export blacklist=""
fi

if [ ${source} == "Ensembl" ]; then
	export genome_suffix=".dna_sm.toplevel.fa.gz"
	export GTF_suffix=".102.gtf.gz"
	export GFF_suffix=".102.gff3.gz"
else
	export genome_suffix="_genomic.fna.gz"
	export GTF_suffix="_genomic.gtf.gz"
	export GFF_suffix="_genomic.gff.gz"
fi

#########################################################################
# Define input files

export sampleID=${pathData}/ChIP-seq/${sp}/${sample}_samples_input.csv
export genome=${pathData}/genome_sequences/${sp}/${spID}${genome_suffix}
export GTF=${pathData}/genome_sequences/${sp}/${spID}${GTF_suffix}
export GFF=${pathData}/genome_sequences/${sp}/${spID}${GFF_suffix}

if [ -f "${GTF}" ]; then
    echo "Using GTF for annotations."
    export annotations="--gtf ${GTF}"
else 
    echo "Using GFF for annotations."
    export annotations="--gff ${GFF}"
fi


if [ 0 -lt $(ls ${pathResults}/indexes/${spID}*.rev.2* 2>/dev/null | wc -w) ]; then
    echo "Indexes already done!"
    export index="--bwa_index ${pathResults}/indexes/ --bowtie2_index ${pathResults}/indexes/"
else
    echo "Need to create indexes!"
    export index="--save_reference"
fi

#########################################################################
echo "#!/bin/bash" > ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}

if [ ${cluster} = "cluster" ]; then
	echo "#SBATCH --job-name=ChIP_calling_${sp}_${sample}" >>  ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --output=${pathScripts}/std_output_peaks_calling_${sp}_${sample}.txt" >>  ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --error=${pathScripts}/std_error_peaks_calling_${sp}_${sample}.txt" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --partition=cpu" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --mem=20G" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
	echo "#SBATCH --time=8:00:00" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
fi

echo "source ${pathConda}" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
echo "conda activate nextflow" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}

echo "nextflow run nf-core/chipseq --input ${sampleID} --outdir ${pathResults}/${sample} --fasta ${genome} ${annotations} ${blacklist} --aligner bowtie2 --macs_gsize ${genomesize} -profile ${container} -with-conda true ${index} --max_memory '20.GB' --max_cpus ${threads} ${resume}" >> ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}

#########################################################################

if [ ${cluster} = "cluster" ]; then
	sbatch ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
else
	bash ${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
fi

#########################################################################

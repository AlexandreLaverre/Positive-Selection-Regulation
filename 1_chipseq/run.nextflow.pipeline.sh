#!/bin/bash

export sp=$1				          # i.e: dog human mouse ...
export sample=$2			        # i.e: Wilson Schmidt Rensch ...
export peaksType=$3           # i.e: Narrow or Broad
export threads=$4			        # i.e: number of threads to use
export cluster=$5			        # i.e: local or cluster
export container=$6			      # i.e: conda or singularity (SLURM works only with singularity)
export resume=${7:-"false"}	  # i.e: resume or false
export skip=${8:-"false"}		  # i.e: skip or false

########################################################################################################################
# Ensure Conda environment exists
if ! conda env list | grep -q '^nextflow'; then
    echo "Creating nextflow Conda environment..."
    conda env create -f env_nextflow.yaml
fi

########################################################################################################################

export path=${path:-"$(pwd)/../../"}
export pathConda="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
export pathResults=${path}/results/peaks_calling/${peaksType}Peaks/${sp}/
export pathData=${path}/data
export pathScripts=${path}/scripts/1_chipseq/logs

mkdir -p ${pathResults}

# Define parameters according to species
source ${path}/scripts/config/params.sh ${sp}

########################################################################################################################
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


if [ 0 -lt "$(ls ${pathResults}/indexes/${spID}*.rev.2* 2>/dev/null | wc -w)" ]; then
    echo "Indexes already done!"
    export index="--bwa_index ${pathResults}/indexes/ --bowtie2_index ${pathResults}/indexes/"
else
    echo "Need to create indexes!"
    export index="--save_reference"
fi

[[ "${skip}" != "false" ]] && skip_flags="--skip_spp --skip_multiqc"
[[ "${resume}" != "false" ]] && resume_flag="-resume"
[[ "${peaksType}" == "Narrow" ]] && peaksType="--narrow_peak"

########################################################################################################################
logFile=${pathScripts}/bsub_ChIP-seq_peaks_calling_${sp}_${sample}
echo "#!/bin/bash" > "${logFile}"

if [ ${cluster} = "cluster" ]; then
  {
  echo "#SBATCH --job-name=ChIP_calling_${sp}_${sample}"
	echo "#SBATCH --output=${pathScripts}/std_output_peaks_calling_${sp}_${sample}.txt"
	echo "#SBATCH --error=${pathScripts}/std_error_peaks_calling_${sp}_${sample}.txt"
	echo "#SBATCH --partition=cpu"
	echo "#SBATCH --mem=50G"
	echo "#SBATCH --cpus-per-task=${threads}"
	echo "#SBATCH --time=20:00:00"
  } >> "${logFile}"

fi

echo "source ${pathConda}" >> "${logFile}"
echo "conda activate nextflow" >> "${logFile}"

echo "nextflow run nf-core/chipseq --input ${sampleID} --outdir ${pathResults}/${sample} --fasta ${genome} \
      ${annotations} ${blacklist} --aligner bowtie2 --macs_gsize ${genomesize} ${peaksType} -profile ${container} \
      -with-conda true ${index} --max_memory '50.GB' --max_cpus ${threads} ${skip_flags} ${resume_flag}" >> "${logFile}"

########################################################################################################################

if [ ${cluster} = "cluster" ]; then
	sbatch  "${logFile}"
else
	bash  "${logFile}"
fi

########################################################################################################################

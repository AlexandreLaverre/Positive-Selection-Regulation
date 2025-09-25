#!/bin/bash

########################################################################################################################
# Default values
threads="1"
system="local"
container="conda"
peaksType="Narrow"
skip="false"
resume="false"

# HELP
function show_help() {
    echo "Usage: ./run.nextflow.sh --sp <species_name> --sample <sample_name> [options]"
    echo
    echo "Options:"
    echo "  --sp        Species (e.g., human, mouse) [required]"
    echo "  --sample    Sample name (e.g., Wilson) [required]"
    echo "  --threads   Number of threads [default: 1]"
    echo "  --system    Execution mode: local or SLURM [default: local]"
    echo "  --container Container type: conda or singularity [default: conda]"
    echo "  --peaksType Narrow or Broad [default: Narrow]"
    echo "  --resume    Resume previous run: true/false [default: false]"
    echo "  --skip      Skip some steps: true/false [default: false]"
    echo "  --help      Show this help message"
    echo
    echo "Example:"
    echo "  ./run.nextflow.sh --sp human --sample Wilson --threads 4 --peaksType Narrow"
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sp) sp="$2"; shift ;;
        --sample) sample="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --system) system="$2"; shift ;;
        --container) container="$2"; shift ;;
        --peaksType) peaksType="$2"; shift ;;
        --resume) resume="$2"; shift ;;
        --skip) skip="$2"; shift ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown parameter: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$sp" || -z "$sample" ]]; then
    echo "Error: --sp and --sample are required."
    show_help
    exit 1
fi

########################################################################################################################
# Check if Nextflow is installed in PATH
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Nextflow not found in PATH. RegEvol_workflows environment will be activated..."
    ExistingNextFlow=false
else
    echo "[INFO] Using system-installed Nextflow: $(command -v nextflow)"
    ExistingNextFlow=true
fi

########################################################################################################################
# Check if FASTQ files exist for this species/sample, else download
export path=${path:-"$(pwd)/../../"}

fastq_dir="${path}/data/ChIP-seq/${sp}/${sample}/"
dl_list="${path}/scripts/docs/ChIP-seq/${sp}/DL_${sample}.txt"

if [ ! -d "$fastq_dir" ] || ! ls "${fastq_dir}"/*.fastq.gz >/dev/null 2>&1; then
    echo "FASTQ files not found for ${sp}/${sample}, downloading..."
    mkdir -p "$fastq_dir"
    if [ -f "$dl_list" ]; then
        while read -r url; do
            echo "Downloading $url..."
            wget -P "$fastq_dir" "$url"
        done < "$dl_list"
    else
        echo "Download list not found: $dl_list"
        exit 1
    fi
else
    echo "FASTQ files already exist for ${sp}/${sample}, skipping download."
fi
########################################################################################################################
# Define paths and variables
export pathConda="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
export pathResults=${path}/results/peaks_calling/${peaksType}Peaks/${sp}/
export pathData=${path}/data
export pathScripts=${path}/scripts/1_chipseq/logs

mkdir -p ${pathResults}

# Define parameters according to species
source ${path}/scripts/config/params.sh ${sp}

########################################################################################################################
# Define input files
export sampleID=${path}/scripts/docs/ChIP-seq/${sp}/${sample}_samples_input.csv
export genome=${pathData}/genome_sequences/${sp}/${spID}${genome_suffix}
export GTF=${pathData}/genome_sequences/${sp}/${spID}${GTF_suffix}
export GFF=${pathData}/genome_sequences/${sp}/${spID}${GFF_suffix}

# Verify required files exist
for f in "$sampleID" "$genome" "$GTF" "$GFF"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

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

if [ ${system} = "SLURM" ]; then
  {
  echo "#SBATCH --job-name=ChIP_calling_${sp}_${sample}"
	echo "#SBATCH --output=${pathScripts}/std_output_peaks_calling_${sp}_${sample}.txt"
	echo "#SBATCH --error=${pathScripts}/std_error_peaks_calling_${sp}_${sample}.txt"
	echo "#SBATCH --partition=cpu"
	echo "#SBATCH --mem=30G"
	echo "#SBATCH --cpus-per-task=${threads}"
	echo "#SBATCH --time=20:00:00"
  } >> "${logFile}"

fi

echo "source ${pathConda}" >> "${logFile}"
if [ "$ExistingNextFlow" = false ]; then
    echo "conda activate RegEvol_workflows" >> "${logFile}"
    [[ "${container}" != "singularity" ]] && mamba install conda-forge::singularity
fi
echo "[[ "$(uname)" == "Darwin" ]] && export CONDA_SUBDIR=osx-64" >> "${logFile}"


echo "nextflow run nf-core/chipseq -r 2.0.0 --input ${sampleID} --outdir ${pathResults}/${sample} --fasta ${genome} \
      ${annotations} ${blacklist} --aligner bowtie2 --macs_gsize ${genomesize} ${peaksType} --conda-frontend mamba  \
      -profile ${container} -with-conda true ${index} --max_memory '30.GB' --max_cpus ${threads} ${skip_flags} ${resume_flag}" >> "${logFile}"

########################################################################################################################

if [ ${system} = "SLURM" ]; then
	sbatch  "${logFile}"
else
	bash  "${logFile}"
fi

########################################################################################################################

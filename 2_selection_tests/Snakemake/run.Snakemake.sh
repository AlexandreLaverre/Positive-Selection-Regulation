#!/bin/bash

########################################################################################################################
# Default values
threads="1"
system="local"
dryRun="false"
unlock="false"
extend_config=""
baseDir="$(pwd)/../../../"

# HELP
function show_help() {
    echo "Usage: ./run.snakemake.sh --sp <species_name> --sample <sample_name> [options]"
    echo
    echo "Options:"
    echo "  --sp              Species (e.g., human, mouse) [required]"
    echo "  --sample          Sample name (e.g., Wilson) [required]"
    echo "  --baseDir         Path to base directory  [default: three levels up from Snakefile]"
    echo "  --threads         Number of threads us[default: 1]"
    echo "  --system          Execution mode: local or SLURM [default: local]"
    echo "  --extend_config   Add any Snakemake extended configuration parameters [default: '']"
    echo "  --dryRun          Run Snakemake in dry-run mode: true/false [default: false]"
    echo "  --unlock          Run Snakemake with the --unlock argument: true/false [default: false]"
    echo
    echo "Example:"
    echo "  ./run.snakemake.sh --sp human --sample Wilson --threads 10 --dryRun --system SLURM --extend_config TF_source='bed'"
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sp) sp="$2"; shift ;;
        --sample) sample="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --system) system="$2"; shift ;;
        --extend_config) extend_config="$2"; shift ;;
        --dryRun) dryRun=true ;;
        --unlock) unlock=true ;;
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

export Prefix=${sp}_${sample}

##################################################################
# Check if Snakemake is installed in PATH
if ! command -v snakemake &> /dev/null; then
    echo "[INFO] Snakemake not found in PATH. RegEvol_workflows environment will be activated..."
    source $(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh
    conda activate RegEvol_workflows
else
    echo "[INFO] Using system-installed Snakemake: $(command -v snakemake)"
fi

##################################################################
export path=${baseDir}
export pathLog="${path}/scripts/2_selection_tests/logs"
[[ "$(uname)" == "Darwin" ]] && export CONDA_SUBDIR=osx-64

##################################################################
echo "[INFO] Running Snakemake with the following parameters:"
echo "       Species:         ${sp}"
echo "       Sample:          ${sample}"
echo "       Threads:         ${threads}"
echo "       System:          ${system}"
echo "       Extended config: ${extend_config}"
echo "       Dry Run:         ${dryRun}"
echo "       Unlock:          ${unlock}"

cmd="snakemake -j 128 --config sp=${sp} sample=${sample} nbPart=${threads} baseDir=${baseDir} system=${system}"
[[ -n "$extend_config" ]] && cmd+=" $extend_config"

cmd+=" --rerun-triggers input --rerun-incomplete --keep-going \
       --use-conda --conda-frontend mamba --conda-prefix .snakemake/conda \
       --cluster 'sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
       -c {params.threads} --mem={params.mem} -t {params.time}'"

$unlock && cmd+=" --unlock"
$dryRun && cmd+=" --dry-run"
eval "$cmd"

##################################################################

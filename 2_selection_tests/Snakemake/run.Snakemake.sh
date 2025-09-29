#!/bin/bash

########################################################################################################################
# Default values
threads="1"
system="local"
dryRun="true"
unlock="false"
baseDir="$(pwd)/../../../"

# HELP
function show_help() {
    echo "Usage: ./run.snakemake.sh --sp <species_name> --sample <sample_name> [options]"
    echo
    echo "Options:"
    echo "  --sp        Species (e.g., human, mouse) [required]"
    echo "  --sample    Sample name (e.g., Wilson) [required]"
    echo "  --baseDir   Path to base directory  [default: three levels up from Snakefile]"
    echo "  --threads   Number of threads [default: 1]"
    echo "  --system    Execution mode: local or SLURM [default: local]"
    echo "  --dryRun    Run Snakemake in dry-run mode: true/false [default: true]"
    echo "  --unlock    Run Snakemake with the --unlock argument: true/false [default: false]"
    echo
    echo "Example:"
    echo "  ./run.snakemake.sh --sp human --sample Wilson --threads 10 --dryRun true"
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sp) sp="$2"; shift ;;
        --sample) sample="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --system) system="$2"; shift ;;
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
[[ " $* " =~ "--unlock" ]] && unlock="--unlock"
[[ " $* " =~ "--dryRun" ]] && dryRun="--dry-run"

echo "[INFO] Running Snakemake with the following parameters:"
echo "       Species:        ${sp}"
echo "       Sample:         ${sample}"
echo "       Base Directory: ${baseDir}"
echo "       Threads:        ${threads}"
echo "       System:         ${system}"
echo "       Dry Run:       ${dryRun}"
echo "       Unlock:        ${unlock}"

cmd="snakemake --rerun-triggers mtime -j 64 --config sp=${sp} sample=${sample} nbPart=${threads} baseDir=${baseDir} \
          system=${system} --rerun-incomplete --use-conda --conda-frontend mamba --conda-prefix .snakemake/conda \
          --cluster 'sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
          -c {params.threads} --mem={params.mem} -t {params.time}'"

$unlock && cmd+=" --unlock"
$dryRun && cmd+=" --dry-run"
eval "$cmd"

##################################################################

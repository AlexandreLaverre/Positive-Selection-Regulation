#!/bin/bash

export species=$1				  # i.e: human dog ...
export sample=$2				  # i.e: Wilson Schmidt12 ...
export nbThreads=$3				  # i.e: int (number of part for parallelization)
export dryRun=${4:-""}		      # i.e: -n or nothing (run snakemake in dry-run mode)"

export Prefix=${species}_${sample}

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
export path=${path:-"$(pwd)/../../"}
export pathLog="${path}/scripts/2_selection_tests/logs"
[[ "$(uname)" == "Darwin" ]] && export CONDA_SUBDIR=osx-64

##################################################################
if [ "${dryRun}" = "-n" ]; then
	echo "Run Snakemake in dry-run mode"
fi

snakemake ${dryRun} --rerun-triggers mtime -j 64 --config sp=${species} sample=${sample} nbPart=${nbThreads} --rerun-incomplete \
          --use-conda --conda-frontend mamba --conda-prefix .snakemake/conda \
          --cluster "sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
          -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################

#!/bin/bash

export species=$1				  # i.e: human dog ...
export sample=$2				  # i.e: Wilson Schmidt12 ...
export nbThreads=$3				  # i.e: int (number of part for parallelization)
export dryRun=${4:-""}		      # i.e: -n or nothing (run snakemake in dry-run mode)"

export Prefix=${species}_${sample}

##################################################################
# Ensure Conda environment exists
if ! conda env list | grep -q '^RegEvol_snakemake\s'; then
    echo "Creating Snakemake Conda environment..."
    conda env create -f envs/RegEvol_snakemake.yaml
fi

export path=${path:-"$(pwd)/../../"}
export pathLog="${path}/scripts/2_selection_tests/logs"
[[ "$(uname)" == "Darwin" ]] && export CONDA_SUBDIR=osx-64

export pathConda="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
source ${pathConda}
conda activate snakemake

##################################################################
if [ "${dryRun}" = "-n" ]; then
	echo "Run Snakemake in dry-run mode"
fi

snakemake ${dryRun} --rerun-triggers mtime -j 64 --config sp=${species} sample=${sample} nbPart=${nbThreads} --rerun-incomplete \
          --use-conda --conda-frontend mamba --conda-prefix .snakemake/conda \
          --cluster "sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
          -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################

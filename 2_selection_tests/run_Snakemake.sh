#!/bin/bash

export species=$1				  # i.e: human dog ...
export sample=$2				  # i.e: Wilson Schmidt12 ...
export cluster=$3				  # i.e: local or cluster
export nbThreads=$4				# i.e: int (number of part for parallelization)
export dryRun=${5:-""}		# i.e: -n or nothing (run snakemake in dry-run mode)"

export Prefix=${species}_${sample}

##################################################################

export path=${path:-"$(pwd)/../../"}
export pathLog="${path}/scripts/detect_positive_selection/logs"

export pathConda="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
source ${pathConda}
conda activate TestPos

##################################################################
if [ "${dryRun}" = "-n" ]; then
	echo "Run Snakemake in dry-run mode"
fi

snakemake ${dryRun} --rerun-triggers mtime -j 64 --config sp=${species} sample=${sample}  \
          cluster=${cluster} nbPart=${nbThreads} --rerun-incomplete \
          --cluster "sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
          -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################

#!/bin/bash

export species=$1				# i.e: human dog ...
export sample=$2				# i.e: CEBPA HNF4A ...
export cluster=$3				# i.e: local or cluster
export nbThreads=$4				# i.e: int (number of part for parallelization)
export nbRand=$5				# i.e: int (randomization)
export dryRun=${6:-""}				# i.e: -n or nothing (run snakemake in dry-run mode"
export alignType=${7:-"MAF"}			# i.e: MAF or pairwise
export method=${8:-"parsimony"}			# i.e: only working with alignType=pairwise

export Prefix=${species}_${sample}

if [ ${dryRun} = "-n" ]; then
	echo "Run Snakemake in dry-run mode"
fi

##################################################################

if [ ${cluster} = "cluster" ]; then
	export pathLog=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/scripts/detect_positive_selection/logs
	source /work/FAC/FBM/DEE/mrobinso/evolseq/Tools/mambaforge/etc/profile.d/conda.sh
else
	export pathLog=/Users/alaverre/Documents/Detecting_positive_selection/scripts/detect_positive_selection/logs
	source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
fi

conda activate TestPos

##################################################################

snakemake ${dryRun} --rerun-triggers input -j 50 --config sp=${species} sample=${sample} AlignType=${alignType} \
          AncMethod=${method} cluster=${cluster} nbPart=${nbThreads} nbRand=${nbRand} --rerun-incomplete \
          --cluster "sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
          -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################

#!/bin/bash

export species=$1				# i.e: human dog ...
export sample=$2				# i.e: CEBPA HNF4A ...
export cluster=$3				# i.e: local or cluster
export nbPart=$4				# i.e: int (number of part for parallelization)
export nbRand=$5				# i.e: int (randomization)
export alignType=${6:-"MAF"}		# i.e: MAF or pairwise
export method=${7:-"parismony"}	# i.e: only working with alignType=pairwise

export Prefix=${species}_${sample}

##################################################################

if [ ${cluster} = "cluster" ]; then
	export pathLog=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/scripts/detect_positive_selection
	source /users/alaverre/Tools/mambaforge/etc/profile.d/conda.sh
else
	export pathLog=/Users/alaverre/Documents/Detecting_positive_selection/scripts/detect_positive_selection
	source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
fi

conda activate snakemake

##################################################################

snakemake -j  --config Sp=${species} Sample=${sample} AlignType=${alignType} AncMethod=${method} cluster=${cluster} nbPart=${nbPart} nbRand=${nbRand} --rerun-incomplete --cluster "sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} -c {params.threads} --mem={params.mem} -t {params.time}"

##################################################################

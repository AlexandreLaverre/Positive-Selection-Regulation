#!/bin/bash

export species=$1				# i.e: human dog ...
export sample=$2				# i.e: CEBPA HNF4A ...
export cluster=$3				# i.e: local or cluster
export nbThreads=$4				# i.e: int (number of part for parallelization)
export dryRun=${5:-""}				# i.e: -n or nothing (run snakemake in dry-run mode)"

export Prefix=${species}_${sample}

if [ "${dryRun}" = "-n" ]; then
	echo "Run Snakemake in dry-run mode"
fi

##################################################################

if [ "${cluster}" = "cluster" ]; then
	export pathLog=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/scripts/detect_positive_selection/logs
	source /work/FAC/FBM/DEE/mrobinso/evolseq/Tools/mambaforge/etc/profile.d/conda.sh
else
	export pathLog=/Users/alaverre/Documents/Detecting_positive_selection/scripts/detect_positive_selection/logs
	source /Users/alaverre/miniconda3/etc/profile.d/conda.sh
fi

conda activate TestPos

##################################################################
# mtime or input or params
snakemake ${dryRun} --rerun-triggers mtime -j 64 --config sp=${species} sample=${sample}  \
          cluster=${cluster} nbPart=${nbThreads} --rerun-incomplete \
          --cluster "sbatch -p cpu -N 1 -o ${pathLog}/slurm.out_${Prefix} -e ${pathLog}/slurm.err_${Prefix} \
          -c {params.threads} --mem={params.mem} -t {params.time}" --touch

#--cleanup-metadata ../results/peaks_calling/NarrowPeaks/human/Wilson/CEBPA.consensus_summits.bed ../results/peaks_calling/NarrowPeaks/human/Wilson/CEBPA.peaks_UCSC_names.bed

##################################################################

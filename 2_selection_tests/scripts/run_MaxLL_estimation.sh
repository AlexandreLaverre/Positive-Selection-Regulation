#!/bin/bash

export species=$1				# i.e: human dog ...
export sample=$2				# i.e: Wilson/CEBPA Schmidt/HNF4A ...
export nbThread=$3				# i.e: int (number of part for parallelization)
export nbBin=$4					# i.e: int (number of bins for all deltas distribution)
export cluster=$5				# i.e: local or cluster
export method=$6				# i.e: local or cluster
export simul=${7:-False}				# i.e: int (number of bins for all deltas distribution)

export Prefix=${species}_${sample/\//_}

export path=${path:-"$(pwd)/../../"}
export pathConda="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"

##################################################################

if [ ${cluster} = "cluster" ]; then
	export pathScripts="${path}/scripts/2_selection_tests/scripts/"

	# Calculate time needed: considering 1500 peaks per hour and per thread
	export time=$(( ($(wc -l < ${path}/results/positive_selection/NarrowPeaks/${species}/${sample}/deltas/ancestral_to_observed_deltaSVM.txt) / (1500*${nbThread}))+1 ))

	# BATCH informations
	echo "#!/bin/sh" > ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --job-name=MaxLL_${Prefix}" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --output=${pathScripts}/logs/std_output_MaxLL_${Prefix}.txt" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --error=${pathScripts}/logs/std_error_MaxLL_${Prefix}.txt" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --partition=cpu" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --mem=16G" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --cpus-per-task=${nbThread}" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "#SBATCH --time=${time}:00:00" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}

	echo "source ${pathConda}" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	echo "conda activate TestPos" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}

	echo "python ${pathScripts}/MaxLL_estimation.py ${species} ${sample} -T ${nbThread} --NbBin ${nbBin} -S ${simul} --binType ${method} --cluster" >> ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}
	sbatch ${pathScripts}/logs/bsub_script_MaxLL_${Prefix}

else
	export pathScripts="${path}/scripts/2_selection_tests/scripts/"
	python ${pathScripts}/MaxLL_estimation.py ${species} ${sample} -T ${nbThread} --NbBin ${nbBin} -S ${simul}
fi

##################################################################

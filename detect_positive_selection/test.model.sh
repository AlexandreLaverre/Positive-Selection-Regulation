#!/bin/bash

export sp=$1
export sample=$2
export firstRunNb=$3
export lastRunNb=$4
export threads=$5

####################################################################################

export path=/Users/alaverre/Documents/Detecting_positive_selection/
export pathResults=${path}/results/tools_tests/${sp}/${sample}/Test-models/
export pathScripts=${path}/scripts/

export pathData=${path}/Tools/JialinTool/data/TFBS_ancestor_focal_sequences/${sp}_sequences/
export BED=${path}/Tools/JialinTool/data/${sp}/${sp}_ChIP-Seq/hsap_${sample}.bed

#########################################################################

for run in $(seq ${firstRunNb} ${lastRunNb})
do
	(
	mkdir -p ${pathResults}/run_${run}

	export log=${pathResults}/run_${run}/log.txt
	
	echo "Generate random seq: run ${run}"
	Rscript ${pathScripts}/GenerateNegativeSeq.R ${sp} ${sample} ${BED} ${pathResults}/run_${run} &> ${log}

	echo "Training the model: run ${run}"
	gkmtrain -l 10 -T ${threads} ${pathResults}/run_${run}/posSet.fa ${pathResults}/run_${run}/negSet.fa ${pathResults}/run_${run}/${sample} &> ${log}
	
	echo "Model prediction: run ${run}"
	gkmpredict -T ${threads} ${path}/results/tools_tests/${sp}/kmer.fa ${pathResults}/run_${run}/${sample}.model.txt ${pathResults}/run_${run}/prediction.txt &> ${log}

	echo "Calculate deltaSVM : run ${run}"
	perl ${pathScripts}/testPosSelec.pl ${pathData}/hsap_${sample}_filtered_ancestor.fa ${pathData}/hsap_${sample}_filtered_focal.fa ${pathResults}/run_${run}/prediction.txt 1000 ${pathResults}/run_${run}/deltaSVM.txt &> ${log}
	
	) &
	
done 


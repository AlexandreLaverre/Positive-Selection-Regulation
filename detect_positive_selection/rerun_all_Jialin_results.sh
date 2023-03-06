#!/bin/bash

sp=$1

path="/Users/alaverre/Documents/Detecting_positive_selection/"
pathData="${path}/Tools/JialinTool/data/${sp}/"
pathResults="${path}/results/positive_selection/rerun_Jialin_corrected/${sp}/"

mkdir -p ${pathResults}

if [ ${sp} = "human" ]; then
	samples=("CEBPA" "HNF4A")
	conditions=()
	
elif [ ${sp} = "mouse" ]; then
	samples=("CEBPA" "HNF4A" "FOXA1")
	conditions=("specific_loss_" "specific_gain_" "conserved_")
	
elif [ ${sp} = "fly" ]; then
	samples=("CTCF")
	conditions=("specific_loss_" "specific_gain_" "conserved_")
fi

####################################################################################

for samp in "${samples[@]}"
do
	echo ${samp}
	for cond in "${conditions[@]}"
	do
		if [ -e ${pathResults}/${samp}_${cond}deltaSVM.txt ]; then
			echo "deltaSVM for ${cond} is already done!"
		else
			echo ${cond}
			./testPosSelec.pl ${pathData}/${sp}_sequences/${samp}_${cond}filtered_ancestor.fa ${pathData}/${sp}_sequences/${samp}_${cond}filtered_focal.fa ${pathData}/${sp}_SVM_model/${samp}/kmer_10_library_weigths.txt 10000 ${pathResults}/${samp}_${cond}deltaSVM.txt
		fi
	done
done


####################################################################################

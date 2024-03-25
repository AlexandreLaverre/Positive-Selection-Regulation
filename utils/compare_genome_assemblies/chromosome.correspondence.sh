#!/bin/bash

export species=$1
export genome1=$2
export genome2=$3
export suffix=$4
export cluster=$5

####################################################################################
if [ ${cluster} = "cluster" ]; then
	export path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel
else
	export path=/Users/alaverre/Documents/Detecting_positive_selection/
fi

export pathScripts=${path}/scripts/compare_genome_assemblies

####################################################################################

perl ${pathScripts}/chromosome.correspondence.pl --pathAssembly1=${genome1} \
      --pathAssembly2=${genome2} --pathOutput=${pathGenomes}/chromosome_correspondence_${suffix}.txt

####################################################################################

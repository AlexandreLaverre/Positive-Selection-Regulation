#!/bin/bash

species=$1
genome1=$(basename -- "$2")
genome2=$(basename -- "$3")
suffix=$4

####################################################################################

export path=${path:-"$(pwd)/../../"}
export pathGenomes=${path}/data/genome_sequences/${species}
export pathScripts=${path}/scripts/utils/compare_genome_assemblies

####################################################################################

perl ${pathScripts}/chromosome.correspondence.pl --pathAssembly1=${pathGenomes}/${genome1} \
      --pathAssembly2=${pathGenomes}/${genome2} --pathOutput=${pathGenomes}/chromosome_correspondence_${suffix}.txt

####################################################################################

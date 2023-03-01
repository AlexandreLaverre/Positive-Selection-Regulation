#!/bin/bash

export species=$1
export genome1=$2
export genome2=$3

####################################################################################

export path=/Users/alaverre/Documents/Detecting_positive_selection/
export pathGenomes=${path}/data/genome_sequences/${species}
export pathScripts=${path}/scripts/compare_genome_assemblies

####################################################################################

perl ${pathScripts}/chromosome.correspondence.pl --pathAssembly1=${pathGenomes}/${genome1} --pathAssembly2=${pathGenomes}/${genome2} --pathOutput=${pathGenomes}/chromosome_correspondence.txt

####################################################################################

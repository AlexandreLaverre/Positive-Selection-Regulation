#!/bin/bash

export species=$1
export type=$2
export nbmax=$3
export cluster=$4

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/mydatalocal/IPLOSS
fi

####################################################################################

export pathGenomes=${path}/data/genome_sequences/${species}
export pathGenome=${path}/data/genome_sequences/${species}/${type}
export pathWGA=${path}/data/whole_genome_alignments/363birds/orthoblocks_withIP_bychr_ref_Cairina_moschata.maf
export pathScripts=${path}/scripts/compare_genome_assemblies

####################################################################################

perl ${pathScripts}/check.genome.assembly.maf.pl --pathGenomeSequence=${pathGenome}/genome_sequence.fa.gz --pathMAF=${pathWGA} --species=${species} --minAlnSize=20 --maxNbAln=${nbmax} --pathOutput=${pathGenomes}/chromosome_correspondence_363birds_${type}.txt

####################################################################################

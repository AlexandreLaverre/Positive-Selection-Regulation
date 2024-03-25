#!/bin/bash

export sp=$1				        # i.e: dog human mouse ...
export sample=$2			      # i.e: Wilson Schmidt Rensch ...
export cluster=$3			      # i.e: local or cluster

########################################################################################################################
if [ "${cluster}" = "cluster" ]; then
  path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/peaks_calling/NarrowPeaks
else
  path=/Users/alaverre/Documents/Detecting_positive_selection/results/peaks_calling/NarrowPeaks
fi

pathPeaks=${path}/${sp}/${sample}/bowtie2/mergedLibrary/macs2/narrowPeak
pathOutput=${path}/${sp}/${sample}

########################################################################################################################
for TF in `ls ${pathPeaks}/consensus`
do
  echo "$TF"
  # Merge all summits and overlap with consensus peaks
  cat ${pathPeaks}/${TF}_*summits.bed > ${pathPeaks}/${TF}_merge_summits.bed
  overlap.py ${pathPeaks}/consensus/${TF}/${TF}.consensus_peaks.bed ${pathPeaks}/${TF}_merge_summits.bed ${pathPeaks}/${TF}_overlap_consensus_max_summits.txt --keep_max --reference_ID

  # Format consensus summits BED file
  cut -f 5 ${pathPeaks}/${TF}_overlap_consensus_max_summits.txt > ${pathPeaks}/${TF}_consensus_summits_ID.txt
  cut -f 4 ${pathPeaks}/${TF}_overlap_consensus_max_summits.txt > ${pathPeaks}/${TF}_consensus_ID.txt
  sed -i "s/:/\t/g" ${pathPeaks}/${TF}_consensus_summits_ID.txt
  paste ${pathPeaks}/${TF}_consensus_summits_ID.txt ${pathPeaks}/${TF}_consensus_ID.txt | tail -n +2 > ${pathOutput}/${TF}.consensus_summits.bed

  # clean tmp files
  rm ${pathPeaks}/${TF}_merge_summits.bed ${pathPeaks}/${TF}_overlap_consensus_max_summits.txt ${pathPeaks}/${TF}_consensus_summits_ID.txt ${pathPeaks}/${TF}_consensus_ID.txt

done

########################################################################################################################

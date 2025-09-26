#!/bin/bash

export sp=$1				        # i.e: dog human mouse ...
export sample=$2			      # i.e: Wilson Schmidt Rensch ...

########################################################################################################################
export path=${path:-"$(pwd)/../../../"}
export pathPeaks=${path}/results/peaks_calling/NarrowPeaks/${sp}/${sample}/bowtie2/mergedLibrary/macs2/narrowPeak
export pathOutput=${path}/results/peaks_calling/NarrowPeaks/${sp}/${sample}

########################################################################################################################
if [ -d "${pathPeaks}/consensus" ]; then
  echo "Consensus directory exists"
  for file in "${pathPeaks}/consensus/"*; do
    TF=$(basename "$file")
    echo "$TF"


    # Merge all summits and overlap with consensus peaks
    cat "${pathPeaks}/${TF}_*summits.bed" > "${pathPeaks}/${TF}_merge_summits.bed"
    python ${path}/scripts/2_selection_tests/scripts/utils/overlap.py "${pathPeaks}/consensus/${TF}/${TF}.consensus_peaks.bed" \
     "${pathPeaks}/${TF}_merge_summits.bed" "${pathPeaks}/${TF}_overlap_consensus_max_summits.txt" --keep_max --reference_ID

    # Format consensus summits BED file
    cut -f 5 "${pathPeaks}/${TF}_overlap_consensus_max_summits.txt" > "${pathPeaks}/${TF}_consensus_summits_ID.txt"
    cut -f 4 "${pathPeaks}/${TF}_overlap_consensus_max_summits.txt" > "${pathPeaks}/${TF}_consensus_ID.txt"
    sed -i "s/:/\t/g" "${pathPeaks}/${TF}_consensus_summits_ID.txt"
    paste "${pathPeaks}/${TF}_consensus_summits_ID.txt" "${pathPeaks}/${TF}_consensus_ID.txt" | tail -n +2 > "${pathOutput}/${TF}.consensus_summits.bed"

    # clean tmp files
    rm "${pathPeaks}/${TF}_merge_summits.bed" "${pathPeaks}/${TF}_overlap_consensus_max_summits.txt" \
    "${pathPeaks}/${TF}_consensus_summits_ID.txt" "${pathPeaks}/${TF}_consensus_ID.txt"

  done
else
  for file in "${pathPeaks}"/*_sample_summits.bed; do
    TF=$(basename "$file" _sample_summits.bed)
    cp "$file" "${pathOutput}/${TF}_consensus_summits.bed"
  done
fi

########################################################################################################################

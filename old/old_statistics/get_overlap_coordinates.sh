#!/bin/bash

# Define the main BED file and the transcription factor BED files with their names
main_bed="all_human_peaks.bed"
tf_names=("Wilson/CEBPA" "Wilson/HNF4A" "Wilson/HNF6" "Wilson/FOXA1" "Schmidt12/CTCF")

bed_files=$(printf "../%s.peaks_UCSC_names.bed " "${tf_names[@]}")
cat $bed_files | sort -k1,1 -k2,2n | bedtools merge > all_human_peaks.bed

# Create a copy of the main BED file for accumulating results
cp "$main_bed" overlap.bed

# Process each TFBS file
for tf in "${tf_names[@]}"; do
    echo $tf
    bed_file="../${tf}.peaks_UCSC_names.bed"

    # Generate binary column: 1 if overlap, 0 otherwise
    bedtools intersect -a all_human_peaks.bed -b "$bed_file" -c | cut -f4 > temp_overlap_column.bed

    # Append this overlap column to the final results
    paste overlap.bed temp_overlap_column.bed > updated_results.bed
    mv updated_results.bed overlap.bed
done

# Add the header at the top of the file
header="Chrom\tStart\tEnd\t$(echo "${tf_names[@]}")"
echo -e "$header" | cat - overlap.bed > all_human_peaks_overlap.bed

# Cleanup temporary files
rm overlap.bed temp_overlap_column.bed

# Count overlap
cut -f4- all_human_peaks_overlap.bed | sort | uniq -c > overlap_counts.txt

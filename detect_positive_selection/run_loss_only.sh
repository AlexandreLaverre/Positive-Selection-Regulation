#!/bin/bash
sp=$1		# mouse
sample=$2	# Wilson
TF=$3		# CEBPA_loss
BED=$4		# Mmusculus_loss_peaks.txt

path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/
pathBED=${path}/results/peaks_overlap/peaks_by_type/${BED}
pathResults=${path}/results/positive_selection/${sp}/${sample}/${TF}/
pathScripts=${path}/scripts/detect_positive_selection

./extract_sequences_from_MAF.sh ${sp} ${sample} ${TF} ${pathBED} cluster


cd ${pathResults}/Alignments/
mkdir -p ${pathResults}/sequences/

# Get all non empty ancestral sequences in one file
find ancestral_sequences -name "*nogap.fa" -size +0 | xargs basename -s _nogap.fa > list_ancestral.txt
find ancestral_sequences -name "*nogap.fa" -size +0 -exec cat {} + > ../sequences/filtered_ancestral_sequences.fa

# Get all corresponding focal sequences in one file
find focal_sequences -name "*nogap.fa" -size +0 -exec cat {} + > ../sequences/all_focal_sequences.fa
seqtk subseq ../sequences/all_focal_sequences.fa list_ancestral.txt > ../sequences/filtered_focal_sequences.fa
        
# Make sequences in uppercase to remove potential soft repeat mask 
awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' ../sequences/filtered_focal_sequences.fa > ../sequences/filtered_focal_sequences_upper.fa


# Create Archive
#tar -czf ${pathResults}/alignments.archive.tar.gz -C ${pathResults} Alignments/
#rm -r ${pathResults}/Alignments/

# Run Positive Selection test
python ${pathScripts}/testPosSelec.py ${sp} ${sample} ${TF} 10000 cluster --NbThread 5

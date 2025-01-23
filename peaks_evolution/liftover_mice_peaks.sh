#!/bin/bash
source /work/FAC/FBM/DEE/mrobinso/evolseq/Tools/mambaforge/etc/profile.d/conda.sh
conda activate cactus

path="/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/"
pathPeaks="${path}/results/peaks_calling/NarrowPeaks/"
pathChain="${path}/data/chain_files/"
pathResults="${path}/results/peaks_overlap/NarrowPeaks/mouse/"

# LiftOver
refs=("musculus" "spretus" "caroli")
targets=("musculus" "spretus" "caroli")
TFs=("CEBPA" "FOXA1" "HNF4A")

for ref in "${refs[@]}"; do
    common=$([[ $ref == "musculus" ]] && echo "mouse" || echo "$ref")
    sample=$([[ $ref == "musculus" ]] && echo "Wilson" || echo "Stefflova")
    
    for target in "${targets[@]}"; do
        [[ $ref == "$target" ]] && continue
        
        for TF in "${TFs[@]}"; do
            # Check if the output file already exists
            [[ -e "${pathResults}/${TF}/lifted_peaks/M${ref}_${TF}_to_M${target}_0.9.bed" ]] && continue

            echo "Lifting ${ref} to ${target} for ${TF}..."
            mkdir -p ${pathResults}/${TF}/lifted_peaks/ 

            liftOver \
                "${pathPeaks}/${common}/${sample}/${TF}.peaks_UCSC_names.bcded" \
                "${pathChain}/Mus_${ref}_to_Mus_${target}.chain.gz" \
                "${pathResults}/${TF}/lifted_peaks/M${ref}_${TF}_to_M${target}_0.9.bed" \
                "${pathResults}/${TF}/lifted_peaks/M${ref}_${TF}_to_M${target}_0.9.unMapped" \
                -minMatch=0.9

        done
    done
done

# Overlap to mouse
targets=("musculus" "spretus")
for target in "${targets[@]}"; do
    for TF in "${TFs[@]}"; do 
        overlap.py "${pathPeaks}/${ref}/${sample}/${TF}.peaks_UCSC_names.bed" \
        "${pathResults}/${TF}/lifted_peaks/M${target}_${TF}_to_M${ref}_0.9.bed" \
        "${pathResults}/${TF}/M${target}_lifted_overlap_M${ref}.txt" --count_overlap -v --reference_ID
    done
done


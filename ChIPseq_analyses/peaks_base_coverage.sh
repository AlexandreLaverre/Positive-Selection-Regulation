#!/bin/bash

export sp=$1				        # i.e: dog human mouse ...
export sample=$2			      # i.e: Wilson Schmidt Rensch ...
export TF=$3			          # i.e: CEBPA CTCF HNF4A ...
export threads=$4

if [[ $sp == "human" ]]; then GenomeSize=2913022398; fi
if [[ $sp == "drosophila" ]]; then GenomeSize=142573017; fi
if [[ $sp == "mouse" ]]; then GenomeSize=2652783500; fi

path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/peaks_calling/NarrowPeaks
pathPeaks=${path}/${sp}/${sample}/bowtie2/mergedLibrary/
pathResults=${pathPeaks}/deepTools/coverage/
mkdir -p "${pathResults}/logs" "${pathResults}/${TF}"

for pathBam in "$pathPeaks"/*"$TF"*bam; do
  bam=$(basename "$pathBam")
  indiv=$(echo "$bam" | cut -d "." -f1 | cut -d "_" -f2)
  input="input_DNA_${indiv}.mLb.clN.sorted.bam"
  logFile="${pathResults}/logs/bsub_BAM_coverage_${sp}_${TF}_${indiv}.sh"

  { echo "#!/bin/bash"
    echo "#SBATCH --job-name=BAM_coverage_${sp}_${TF}_${indiv}"
    echo "#SBATCH --output=${pathResults}/logs/std_output_BAM_coverage_${sp}_${TF}_${indiv}.txt"
    echo "#SBATCH --error=${pathResults}/logs/std_error_BAM_coverage_${sp}_${TF}_${indiv}.txt"
    echo "#SBATCH --partition=cpu"
    echo "#SBATCH --mem=10G"
    echo "#SBATCH --cpus-per-task=${threads}"
    echo "#SBATCH --time=20:00:00"
    echo "source /work/FAC/FBM/DEE/mrobinso/evolseq/Tools/mambaforge/etc/profile.d/conda.sh"
    echo "conda activate deeptools"

    # Obtain normalized reads count
    echo "bamCompare -b1 ${pathBam} \
    -b2 ${pathPeaks}/${input} \
    -o ${pathResults}/${TF}/${indiv}_bgNorm_RPKM.bw \
    --binSize 1 \
    --scaleFactorsMethod None \
    --normalizeUsing RPKM
    --effectiveGenomeSize ${GenomeSize} \
    -p ${threads} 2> ${pathPeaks}/deepTools/coverage/logs/${TF}_${indiv}_bamCompare.log"
    #    --normalizeUsing RPKM --centerReads --sampleLength 1000 --scaleFactorsMethod readCount/None \

    # Matrix for each sample
    echo "computeMatrix reference-point --referencePoint center \
    -a 200 -b 200 \
    -bs 1 \
    -R ${pathPeaks}/macs2/narrowPeak/consensus/${TF}/${TF}.consensus_peaks.bed \
    -S ${pathResults}/${TF}/${indiv}_bgNorm_RPKM.bw \
    -o ${pathResults}/${TF}/${indiv}_matrix_RPKM_200bp.gz \
    -p ${threads} \
    --outFileSortedRegions ${pathResults}/${TF}/${indiv}_peaks_RPKM_200bp.bed"

  } > "${logFile}"

  sbatch "${logFile}"
done






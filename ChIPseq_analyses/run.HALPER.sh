#!/bin/bash

export sp=$1				        # i.e: dog human mouse ...
export TF=$2			          # i.e: CEBPA CTCF HNF4A ...
export cluster=$3 		      # i.e: local or cluster

########################################################################################################################
# Define path and variables
path=/beegfs/data/alaverre/PosSel/
source ${path}/scripts/config/params.sh "${sp}" "${cluster}"
pathData=${path}/results/peaks_calling/${sp}/
pathResults=${path}/results/homologous_peaks/${sp}/
pathHAL=/beegfs/banque/hal/241-mammalian-2020v2.hal

all_species=("Homo_sapiens" "Macaca_mulatta" "Mus_musculus" "Mus_spretus" "Mus_caroli" "Felis_catus" \
             "Canis_lupus_familiaris" "Oryctolagus_cuniculus" "Rattus_norvegicus" "Bos_taurus")

# Remove reference species from all to define target species separated by comma
other_species=($(IFS=" "; echo "${all_species[@]/$sp}"))

for target in "${other_species[@]}"; do
  OutputFile=${pathResults}/HALPER_${TF}_${sp}2${target}.bed
  if [ -e "${OutputFile}" ]; then
    echo "${TF}: ${sp} to ${target} HALPER already done!"
  else
  ######################################################################################################################
  # Write instructions for cluster
  mkdir -p "${pathResults}"
  mkdir -p "${path}/scripts/ChIPseq_analyses/logs/"
  logFile=${path}/scripts/ChIPseq_analyses/logs/bsub_HALPER_"${sp}"_"${TF}"_to_"${target}"

  echo "#!/bin/bash" > "${logFile}"
  if [ "${cluster}" = "cluster" ]; then
    {
    echo "#SBATCH --job-name=HALPER_${TF}_${sp}_to_${target}"
    echo "#SBATCH --partition=cpu"
    echo "#SBATCH --mem=5G"
    echo "#SBATCH --cpus-per-task=2"
    echo "#SBATCH --time=8:00:00"
    } >> "${logFile}"
  fi

  ######################################################################################################################
  # Write instructions for HALPER
  {
    # Lift summits
    echo "halLiftover --bedType 4 ${pathHAL} ${sp} ${pathData}/${TF}.consensus_summits.bed ${target} \
          ${pathResults}/${TF}.consensus_summits_to_${target}.bed &"

    # Lift peaks
    echo "halLiftover --bedType 4 ${pathHAL} ${sp} ${pathData}/${TF}.consensus_peaks.bed ${target} \
          ${pathResults}/${TF}.consensus_peaks_to_${target}.bed"
    echo "wait"

    # Merge lifted regions with HALPER
    echo "python3 -m orthologFind -max_frac 1.2 -min_frac 0.8 -protect_dist 5 -qFile ${pathData}/${TF}.consensus_peaks.bed \
          -tFile ${pathResults}/${TF}.consensus_peaks_to_${target}.bed -sFile ${pathResults}/${TF}.consensus_summits_to_${target}.bed \
          -oFile ${pathResults}/HALPER_${TF}_${sp}2${target}.bed -mult_keepone"

  } >> "${logFile}"

  ######################################################################################################################
  # Run script
  if [ "${cluster}" = "cluster" ]; then
    echo "sbatch ${logFile}"
  else
    bash "${logFile}"
  fi
  fi
done

########################################################################################################################

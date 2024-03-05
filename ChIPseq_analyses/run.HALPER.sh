#!/bin/bash

export sp=$1				        # i.e: dog human mouse ...
export sample=$2			      # i.e: Wilson Schmidt Rensch ...
export TF=$3			          # i.e: CEBPA CTCF HNF4A ...
export cluster=$4			      # i.e: local or cluster

########################################################################################################################
# Define path and variables
path=/beegfs/data/alaverre/PosSel/
source ${path}/scripts/params.sh "${sp}" "${cluster}"
pathData=${path}/results/peaks_calling/${sp}/${sample}/${TF}/
pathResults=${path}/results/homologous_peaks/${sp}/${sample}/${TF}/
pathHAL=/beegfs/banque/hal/241-mammalian-2020v2.hal

all_species=("Homo_sapiens" "Macaca_mulatta" "Mus_musculus" "Mus_spretus" "Mus_caroli" "Felis_catus" \
             "Canis_lupus_familiaris" "Oryctolagus_cuniculus" "Rattus_norvegicus" "Bos_taurus")

# Remove reference species from all to define target species separated by comma
other_species=$(IFS=,; echo "${all_species[@]/$sp_name}")

########################################################################################################################
# Write instruction for HALPER
logFile=${path}/scripts/homologous_peaks/logs//bsub_HALPER_"${sp}"_"${sample}"_"${TF}"

echo "#!/bin/bash" > "${logFile}"
if [ "${cluster}" = "cluster" ]; then
	{
  echo "#SBATCH --job-name=HALPER_${sp}_${sample}_${TF}"
	echo "#SBATCH --partition=cpu"
	echo "#SBATCH --mem=10G"
	echo "#SBATCH --cpus-per-task=1"
	echo "#SBATCH --time=12:00:00"
	echo "#SBATCH --array=1-9"
	} >> "${logFile}"
fi

echo "source ${pathConda}" >> "${logFile}"
echo "conda activate hal" >> "${logFile}"

echo "halper_map_peak_orthologs.sh -b ${pathData}/.narrowPeak -o ${pathResults} \
      -s ${sp_name} -t ${other_species} -c ${pathHAL}" >> "${logFile}"

########################################################################################################################
# Run HALPER
if [ "${cluster}" = "cluster" ]; then
	sbatch "${logFile}"
else
	bash "${logFile}"
fi

########################################################################################################################

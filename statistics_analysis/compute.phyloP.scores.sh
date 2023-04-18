#!/bin/bash

################ Input parameters #####################
export species=$1       			  # i.e : human or mouse
export sample=$2				  # i.e : Wilson/CEBPA
export masked_exons=${3:-"FALSE"}   # i.e : TRUE or FALSE

################ Export paths #####################

export path=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/

export pathphyloP=${path}/data/phyloP/${species}/
export pathBED=${path}/results/peaks_calling/${species}/${sample}.peaks.bed_UCSC_names
export pathResults=${path}/results/phyloP/${species}/${sample}
export pathScripts=${path}/scripts/statistics_analysis

mkdir -p ${pathResults}

################ Define aligments input files #####################

if [ "${species}" = "mouse" ]; then
  export chromosomes=(chr{1..19} chrX chrY)
  export way=60way.glire
elif [ "${species}" = "human" ]; then
  export chromosomes=(chr{1..22} chrX chrY)
  export way=17way
fi

export suffixPhylo=phyloP${way}.wigFix.gz
export coordConvention=0_open_end

#useless for the moment
if [ "${masked_exons}" = "TRUE" ]; then
    export suffixExons="_MaskedExons_Ensembl94"
    export pathExons=${path}/data/exons/${species}_exons_Ensembl94.txt
else
    export suffixExons=""
    export pathExons="NA"
fi

################ Running Compute phyloP scores #####################

for chr in "${chromosomes[@]}"
do
    export pathPhyloP_scores=${pathphyloP}/${chr}.${suffixPhylo}
    echo ${chr}
    if [ -e "${pathPhyloP_scores}" ]; then
      if [ -e "${pathResults}/${chr}_${way}${suffixExons}.txt" ]; then
	    echo "already done"
      else
    	    echo "#!/bin/bash" > ${pathScripts}/log/${sample}_${way}_${species}_${chr}${suffixExons}
	    
	    echo "perl ${pathScripts}/compute.phyloP.scores.pl --pathCoords=${pathBED} --coordConvention=${coordConvention} --pathMaskExonBlocks=${pathExons} --pathScores=${pathPhyloP_scores} --chr=${chr} --pathOutput=${pathResults}/${chr}_${way}${suffixExons}.txt" >> ${pathScripts}/log/${sample}_${way}_${species}_${chr}${suffixExons}

      #bash ${pathScripts}/log/${sample}_${way}_${species}_${chr}${suffixExons}
      sbatch -p cpu --time=1:00:00 --mem=30GB -c 1 -o ${pathScripts}/log/${sample}_${way}_${species}_${chr}${suffixExons}_std_output -e ${pathScripts}/log/${sample}_${way}_${species}_${chr}${suffixExons}_std_error ${pathScripts}/log/${sample}_${way}_${species}_${chr}${suffixExons} 
       fi
    fi
done

#####################################################################

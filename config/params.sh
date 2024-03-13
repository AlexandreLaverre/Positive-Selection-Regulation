#!/bin/bash

sp=$1
cluster=$2

if [ ${cluster} = "local" ]; then
	export pathData=/Users/alaverre/Documents/Detecting_positive_selection/data
	export pathConda=/Users/alaverre/miniconda3/etc/profile.d/conda.sh
else
	export pathData=/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/data
	export pathConda=/work/FAC/FBM/DEE/mrobinso/evolseq/Tools/mambaforge/etc/profile.d/conda.sh
fi

########################################################
EnsemblRelease=102

if [ ${sp} = "dog" ]; then
	spID="Canis_lupus_familiaris.CanFam3.1"
	sp_name="Canis_lupus_familiaris"
	close_species="Canis_lupus_familiaris,Lycaon_pictus,Vulpes_vulpes" #Canis_lupus_lupus
	anc_name="Anc04" # Anc09
	chroms=({1..38} "MT" "X")
	export genomesize=2237684358 # from 50: https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/conf/igenomes.config
fi

if [ ${sp} = "human" ]; then
	spID="Homo_sapiens.GRCh38"
	sp_name="Homo_sapiens"
	close_species="Homo_sapiens,Pan_troglodytes,Gorilla_gorilla"
	anc_name="fullTreeAnc105"
	chroms=({1..22} "X" "Y")
	export genomesize=2701262066 
fi

if [ ${sp} = "mouse" ]; then
	spID="Mus_musculus.GRCm38"
	sp_name="Mus_musculus"
	close_species="Mus_musculus,Mus_spretus,Mus_caroli"
	anc_name="fullTreeAnc35"
	chroms=({1..19} "X" "Y")
	genomesize=2307679482
fi

if [ ${sp} = "macaca" ]; then
	spID="sup2kb_GCA_000772875.3_Mmul_8.0.1"
	sp_name="Macaca_mulatta"
	close_species="Macaca_mulatta,Macaca_fascicularis,Macaca_nemestrina"
	anc_name="fullTreeAnc91"
	chroms=(CM002977.3 CM002980.3 CM002984.2 CM002983.2 CM002981.2 CM002982.3 CM002991.3 CM002985.3 CM002987.3 CM002992.3 CM002989.3 CM002979.2 CM002978.2 CM002988.2 CM002986.2 CM002994.2 CM002990.2 CM002995.2 CM002996.3 CM002993.2 CM002997.3 CM003438.1) # from NCBI (1 to 20, X and Y)
	export genomesize=2498932238
fi

if [ ${sp} = "chicken" ]; then
	spID="Gallus_gallus.GRCg6a"
	sp_name="Gallus_gallus"
	close_species="Gallus_gallus,Coturnix_japonica,Meleagris_gallopavo"
	anc_name="birdAnc333"
	chroms=({1..33} "Z" "W" "MT")
	export genomesize=974987959
fi

if [ ${sp} = "rat" ]; then
	spID="Rattus_norvegicus.Rnor_6.0"
	sp_name="Rattus_norvegicus"
	close_species="Rattus_norvegicus,Mus_musculus,Acomys_cahirinus"
	anc_name="fullTreeAnc38"
	chroms=({1..20} "X" "Y")
	export genomesize=2375372135
fi

if [ ${sp} = "cat" ]; then
	spID="sup2kb_GCA_000181335.3_Felis_catus_8.0"
	sp_name="Felis_catus"
	close_species="Felis_catus,Felis_nigripes,Puma_concolor"
	anc_name="fullTreeAnc202"
	chroms=(CM0013{78..96}.2) # from NCBI (A1 to F2 and X)
	export genomesize=2300000000
fi

if [ ${sp} = "cow" ]; then
	spID="GCA_000003205.6_Btau_5.0.1"
	sp_name="Bos_taurus"
	close_species="Bos_taurus,Bos_indicus,Bos_mutus"
	anc_name="fullTreeAnc170"
	chroms=(CM000177.6 CM000178.7 CM000179.6 CM000180.6 CM000181.7 CM000182.7 CM000183.7 CM000184.6 CM000185.7 CM000186.6 CM000187.6 CM000188.7 CM000189.7 CM000190.6 CM000191.7 CM000192.5 CM000193.6 CM000194.7 CM000195.6 CM000196.6 CM000197.6 CM000198.6 CM000199.8 CM000200.7 CM000201.6 CM000202.7 CM000203.6 CM000204.6 CM000205.6 CM000206.6 CM001061.2) # from NCBI (1 to 29, X and Y)
	export genomesize=2370644326
fi

if [ ${sp} = "rabbit" ]; then
	spID="Oryctolagus_cuniculus.OryCun2.0"
	sp_name="Oryctolagus_cuniculus"
	close_species="Oryctolagus_cuniculus,Lepus_americanus,Ochotona_princeps"
	anc_name="fullTreeAnc15"
	chroms=({1..21} "X")
	export genomesize=2300000000 
fi

if [ ${sp} = "spretus" ]; then
	spID="Mus_spretus.SPRET_EiJ_v1"
	sp_name="Mus_spretus"
	close_species="Mus_spretus,Mus_musculus,Mus_caroli"
	anc_name="fullTreeAnc35"
	chroms=(CM0040{94..99}.1 CM004{100..113}.1) #Ensembl:{1..19} "X"
	export genomesize=2307679482 # using the same as Mus_musculus
fi

if [ ${sp} = "caroli" ]; then
	spID="Mus_caroli.CAROLI_EIJ_v1.1"
	sp_name="Mus_caroli"
	close_species="Mus_caroli,Mus_musculus,Mus_spretus"
	anc_name="fullTreeAnc36"
	chroms=(LT608{229..248}.1)  #Ensembl: {1..19} "X"
	export genomesize=2307679482 # using the same as Mus_musculus
fi

########################################################
# Files prefix and suffix according to source

if [ ${sp} = "macaca" ] || [ ${sp} = "cat" ] || [ ${sp} = "cattle" ]; then
	# Anotations from NCBI
	export genome_suffix="_genomic.fna.gz"
	export GTF_suffix="_genomic.gtf.gz"
	export GFF_suffix="_genomic.gff.gz"
else
	# Anotations from Ensembl
	export genome_suffix=".dna_sm.toplevel.fa.gz"
	export GTF_suffix=".${EnsemblRelease}.gtf.gz"
	export GFF_suffix=".${EnsemblRelease}.gff3.gz"
fi

# Blacklisted regions files for ChIP-seq calling
if [ ${sp} = "human" ] || [ ${sp} = "mouse" ]; then
	export blacklist="--blacklist  ${pathData}/genome_sequences/${sp}/blacklist.txt"
fi

########################################################

# RegEvol

RegEvol is a computational framework for detecting directional selection in regulatory sequences through phenotypic predictions and phenotype-to-fitness models. This ongoing project aims to identify signatures of positive selection acting on transcription factor binding sites across mammals and drosophila, using both statistical and evolutionary approaches.

# Project Overview
The pipeline includes the following key steps:

### Regulatory Data Processing
- **TF ChIP-seq Peak Calling**: Extraction and processing of transcription factor binding sites from publicly available ChIP-seq datasets in mammals and drosophila.
- **Homologous Sequence Retrieval**: Identification of orthologous regulatory regions using whole-genome alignments from the Zoonomia Consortium.

### Predictive Modeling
- **gkm-SVM Training**: Training of gapped k-mer support vector machine (gkm-SVM) models on species-specific regulatory sequences to predict binding affinity changes.

### Positive Selection Detection
- **Permutation Test**: Application of the method from [Liu & Robinson-Rechavi, 2020](https://www.science.org/doi/full/10.1126/sciadv.abc9863) to detect extreme change in binding affinity.
- **RegSel Test**: A novel maximum likelihood-based test (RegEvol) to infer evolutionary regimes of regulatory regions from predicted binding affinity shifts.

### Comparative and Functional Analyses
- **Simulations**: Sequence evolution simulations under controlled parameters to validate the performance and robustness of selection detection methods.
- **Conservation Metrics**: Comparison of selection signals with evolutionary conservation scores (phastCons and phyloP).
- **Functional Enrichment**: Gene Ontology enrichment analysis of genes associated with positively selected regulatory regions.
- **Divergence vs. Polymorphism**: Integration of interspecies divergence and intraspecies polymorphism data to interpret evolutionary dynamics.

## Requirements

Before running the pipelines, ensure the following tools are installed on your system:

- **[Conda / Mamba](https://docs.conda.io/en/latest/)** – for managing environments and dependencies.
- **[Singularity](https://sylabs.io/docs/)** – required to use Docker/Singularity containers in Snakemake rules.

This workflow also depends on:

- **[Snakemake](https://snakemake.readthedocs.io/)** ≥ X.Y : to run the RegEvol pipeline.
- **[Nextflow](https://www.nextflow.io/)** : to run the ChIP-seq calling pipeline.

### Option 1: Use your existing installation  
If you already have **Snakemake** and **Nextflow** installed on your system, you can skip creating the bundled environment and run the workflow directly.

### Option 2: Create the provided environment  
If you don’t, you can create the environment with:
```bash
mamba env create -f config/workflows.yaml
```

import os


def get_TFs(config, pathPeaks):
    """Retrieve TFs from the config file or BED file"""
    if config.get("TF_source") == "config":
        return config["TFs"][config["sample"]]
    else:
        bed_file = os.path.join(pathPeaks, "FlyTFPeaksPrimaryTargets.tsv")
        with open(bed_file) as f:
            next(f)  # Skip header
            return sorted({line.split("\t")[3].strip() for line in f})  # Extract unique TFs


def get_localrules(config):
    """Define local rules based on cluster mode"""
    if config["cluster"] == "cluster":
        return ["all", "GetPeaks", "SubSetPeaks", "BED_split", "ConcatSeq", "ConsensusSummits",
                "ModelPrediction", "ChromosomeCorrespondence", "ConvertCoordinates",
                "DownloadVCF", "MergeAllChromosome"]
    else:
        return ["all", "GetPeaks", "SubSetPeaks", "GenerateNegativeSeq", "ModelTraining",
                "ModelValidation", "ModelPrediction", "BED_split",
                "InferAncestralPairwise", "GetSequencesMultiple", "ConcatSeq",
                "PermutationTest", "ArchiveAlignments", "MergeAllChromosome"]


def get_ruleorder(config):
    """Define rule order based on AlignType and TF_source"""
    ruleorder = []
    if config["AlignType"] == "pairwise":
        ruleorder.append("InferAncestralPairwise > GetSequencesMultiple")
    else:
        ruleorder.append("GetSequencesMultiple > InferAncestralPairwise")

    if config.get("TF_source") == "config":
        ruleorder.append("GetPeaks > SubSetPeaks")
    else:
        ruleorder.append("SubSetPeaks > GetPeaks")

    return ruleorder

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


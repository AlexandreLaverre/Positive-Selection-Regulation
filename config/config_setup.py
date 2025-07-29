import os


def get_TFs(config, pathPeaks):
    """Determine the list of TFs"""

    # If TFs is passed directly via command line (e.g. as a string or list)
    if isinstance(config.get("TFs"), (str, list)):
        if isinstance(config["TFs"], str):
            return config["TFs"].split(",")
        return config["TFs"]

    # If TF_source is set to "bed", read from the large BED file
    if config.get("TF_source") == "bed":
        bed_file = os.path.join(pathPeaks, "FlyTFPeaksPrimaryTargets.tsv")
        with open(bed_file) as f:
            next(f)  # Skip header
            return sorted({line.split("\t")[4].strip() for line in f})  # Extract unique TFs

    # Fall back to config file
    else:
        return config["TFs"][config["sample"]]


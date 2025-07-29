import os


def get_TFs(config, pathPeaks):
    """Determine the list of TFs"""
    # If TF_source is set to "bed", read from the large BED file
    if config.get("TF_source") == "bed":
        bed_file = os.path.join(pathPeaks, "FlyTFPeaksPrimaryTargets.tsv")
        with open(bed_file) as f:
            next(f)  # Skip header
            return sorted({line.split("\t")[4].strip() for line in f})  # Extract unique TFs

    else:
        # If TFs is passed directly via command line
        if isinstance(config.get("TFs"), (str, list)):
            override_value = config["TFs"].split(",") if isinstance(config["TFs"], str) else config["TFs"]
            config["TFs"] = {config["sample"]: override_value}

        return config["TFs"][config["sample"]]


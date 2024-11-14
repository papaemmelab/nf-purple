#!/usr/bin/env python3

import argparse
import shutil
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    description=(
        "Bin cobalt probes with similar LogR values "
        "together to decrease oversegmentation."
    )
)
parser.add_argument(
    "--in_pcf",
    type=str,
    required=True,
    help="Path to the input cobalt ratio .pcf file.",
)
parser.add_argument(
    "--bin_probes", 
    type=int, 
    required=True,
    help="Max probe bin size."
)
parser.add_argument(
    "--bin_log_r",
    type=float,
    required=True,
    help="Max probe logR difference to bin."
)
args = parser.parse_args()

cobalt_ratio_pcf = pd.read_csv(args.in_pcf, sep="\t")
cobalt_ratio_pcf_probes = pd.DataFrame(columns=cobalt_ratio_pcf.columns)

# First bin by probes
chrom_arm = None
last_idx = None
for idx, seg in cobalt_ratio_pcf.iterrows():
    if chrom_arm != "_".join(seg[["chrom", "arm"]].astype(str)):
        chrom_arm = "_".join(seg[["chrom", "arm"]].astype(str))
        cobalt_ratio_pcf_probes = pd.concat(
            [cobalt_ratio_pcf_probes, seg.to_frame().T], ignore_index=True
        )
        last_idx = cobalt_ratio_pcf_probes.index[-1]
        continue
    if (
        cobalt_ratio_pcf_probes.loc[last_idx, "n.probes"] <= args.bin_probes
        or seg["n.probes"] <= args.bin_probes
    ):
        means = [
            cobalt_ratio_pcf_probes.loc[last_idx, "mean"]
        ] * cobalt_ratio_pcf_probes.loc[last_idx, "n.probes"]
        means.extend([seg["mean"]] * seg["n.probes"])
        cobalt_ratio_pcf_probes.loc[last_idx, "mean"] = np.mean(means)
        cobalt_ratio_pcf_probes.loc[last_idx, "n.probes"] += seg["n.probes"]
        cobalt_ratio_pcf_probes.loc[last_idx, "end.pos"] = seg["end.pos"]
    else:
        cobalt_ratio_pcf_probes = pd.concat(
            [cobalt_ratio_pcf_probes, seg.to_frame().T], ignore_index=True
        )
        last_idx = cobalt_ratio_pcf_probes.index[-1]

# Then bin by logR mean
cobalt_ratio_pcf_probes = cobalt_ratio_pcf_probes.reset_index().drop(columns="index")
cobalt_ratio_pcf_probes_logR = pd.DataFrame(columns=cobalt_ratio_pcf_probes.columns)
chrom_arm = None
for idx, seg in cobalt_ratio_pcf_probes.iterrows():
    if chrom_arm != "_".join(seg[["chrom", "arm"]].astype(str)):
        chrom_arm = "_".join(seg[["chrom", "arm"]].astype(str))
        cobalt_ratio_pcf_probes_logR = pd.concat(
            [cobalt_ratio_pcf_probes_logR, seg.to_frame().T], ignore_index=True
        )
        last_idx = cobalt_ratio_pcf_probes_logR.index[-1]
        continue
    if (
        abs(cobalt_ratio_pcf_probes.loc[last_idx, "mean"] - seg["mean"])
        <= args.bin_log_r
    ):
        means = [
            cobalt_ratio_pcf_probes_logR.loc[last_idx, "mean"]
        ] * cobalt_ratio_pcf_probes_logR.loc[last_idx, "n.probes"]
        means.extend([seg["mean"]] * seg["n.probes"])
        cobalt_ratio_pcf_probes_logR.loc[last_idx, "mean"] = np.mean(means)
        cobalt_ratio_pcf_probes_logR.loc[last_idx, "n.probes"] += seg["n.probes"]
        cobalt_ratio_pcf_probes_logR.loc[last_idx, "end.pos"] = seg["end.pos"]
    else:
        cobalt_ratio_pcf_probes_logR = pd.concat(
            [cobalt_ratio_pcf_probes_logR, seg.to_frame().T], ignore_index=True
        )
        last_idx = cobalt_ratio_pcf_probes_logR.index[-1]

# store input with another name to replace original
shutil.move(args.in_pcf, args.in_pcf.replace(".pcf", ".original.pcf"))
cobalt_ratio_pcf_probes_logR.to_csv(args.in_pcf, sep="\t", index=False)

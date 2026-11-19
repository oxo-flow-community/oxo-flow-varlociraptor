"""Write the somalier groups file (one line per sample group: its comma-joined
sample names).

Port of the upstream workflow/scripts/somalier_groups.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.11.0), adapted to argv:
  somalier_groups.py --samples '{"columns": [...], "data": [...]}' --output GROUPS.txt --log LOG

Upstream writes `snakemake.params.samples.groupby("group").apply(
lambda df: ",".join(df["sample_name"]), include_groups=False).to_csv(
index=False, header=False, sep="\t")`; the samples sheet is passed as the
frame-encoded JSON (same encoding as the datavzrd variables files).
"""

import argparse
import json
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, help="JSON-encoded samples DataFrame")
    parser.add_argument("--output", required=True, help="Output groups file")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    if args.log:
        sys.stderr = open(args.log, "w")

    samples_data = json.loads(args.samples)
    samples = pd.DataFrame(samples_data["data"], columns=samples_data["columns"])

    samples.groupby("group").apply(
        lambda df: ",".join(df["sample_name"]), include_groups=False
    ).to_csv(args.output, index=False, header=False, sep="\t")


if __name__ == "__main__":
    main()

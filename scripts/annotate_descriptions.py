"""Merge COSMIC signature descriptions into the joined signature table.

Port of the upstream workflow/scripts/annotate_descriptions.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
  annotate_descriptions.py SIGNATURES.tsv DESCRIPTIONS.tsv --output OUT.tsv --log LOG

The upstream snakemake script API is replaced by positional inputs + flags;
the left merge on Signature and the "Signature: Description" relabeling are
unchanged. Diagnostics go to the log file (sys.stderr redirect, as upstream).
"""

import argparse
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sig", help="Joined signature frequency table")
    parser.add_argument("desc", help="COSMIC signature description table")
    parser.add_argument("--output", required=True, help="Output annotated table")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    sys.stderr = open(args.log, "w")

    signatures_df = pd.read_csv(args.sig, sep="\t")
    description_df = pd.read_csv(args.desc, sep="\t")
    signatures_df = pd.merge(signatures_df, description_df, how="left", on="Signature")
    signatures_df["Signature"] = signatures_df.apply(
        lambda row: f"{row['Signature']}: {row['Description']}", axis=1
    )

    signatures_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()

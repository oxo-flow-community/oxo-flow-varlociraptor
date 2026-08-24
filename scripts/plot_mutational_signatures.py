"""Plot mutational signature exposures and mutation counts (altair).

Port of the upstream workflow/scripts/plot_mutational_signatures.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
  plot_mutational_signatures.py SIGNATURES.tsv COUNTS.tsv --output OUT.html --log LOG

The upstream snakemake script API is replaced by positional inputs + flags.
The area chart of signature frequencies (reversed Minimum VAF axis) layered
with the mutation-count line is unchanged. Diagnostics go to the log file
(sys.stderr redirect, as upstream). NOTE: the upstream script uses sys
without importing it (NameError when the stderr redirect executes); the port
adds the missing import.
"""

import argparse
import sys

import altair as alt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("signatures", help="Annotated signature frequency table")
    parser.add_argument("counts", help="Mutation counts table")
    parser.add_argument("--output", required=True, help="Output HTML chart")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    sys.stderr = open(args.log, "w")

    signatures_df = pd.read_csv(args.signatures, sep="\t")

    signatures = (
        alt.Chart(signatures_df)
        .mark_area(interpolate="monotone")
        .encode(
            x=alt.X("Minimum VAF:Q", scale=alt.Scale(reverse=True)),
            y=f"Frequency:Q",
            color="Signature:N",
            tooltip="Description:N",
        )
    )

    mut_counts_df = pd.read_csv(args.counts, sep="\t")

    counts = (
        alt.Chart(mut_counts_df)
        .mark_line(interpolate="basis", color="black")
        .encode(
            x=alt.X("Minimum VAF:Q", scale=alt.Scale(reverse=True)),
            y="Mutation Count:Q",
        )
    )

    final_chart = alt.layer(signatures, counts).resolve_scale(y="independent")

    final_chart.save(args.output)


if __name__ == "__main__":
    main()

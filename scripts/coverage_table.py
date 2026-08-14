"""Compute the per-gene average coverage table.

Port of the upstream workflow/scripts/coverage_table.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
  coverage_table.py COVERAGE.bed... --output OUT.csv --min-cov 5 --log LOG

Each input BED line is chromosome, start, end, gene, coverage (produced by
bedtools_merge). Genes below the minimum average coverage are skipped; the
output is a wide TSV with one row per (chromosome, gene).
"""

import argparse
import sys

import numpy as np
import pandas as pd


def add_missing_columns(df, samples):
    missing_columns = set(samples).difference(df.columns)
    df[list(missing_columns)] = np.nan
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("beds", nargs="+", help="Filtered region BED files")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--min-cov", type=float, default=0, help="Minimum average coverage")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    if args.log:
        sys.stderr = open(args.log, "w")

    group_regions = dict()
    samples = []
    for bed in args.beds:
        sample = bed.split("/")[-1].split(".")[0]
        samples.append(sample)
        with open(bed, "r") as covered_regions:
            for line in covered_regions:
                line = line.strip().split("\t")
                chromosome = line[0]
                gene = line[3]
                coverage = line[4]
                if float(coverage) < args.min_cov:
                    continue
                if (chromosome, gene) not in group_regions:
                    group_regions[(chromosome, gene)] = dict()
                group_regions[(chromosome, gene)][sample] = coverage

    if bool(group_regions):
        df = pd.DataFrame.from_dict(group_regions).T
        df.index.names = ("chromosome", "gene")
        df.reset_index(inplace=True)
        df = add_missing_columns(df, samples)
    else:
        df = pd.DataFrame(columns=["chromosome", "gene"] + samples)

    with open(args.output, "w") as csv_file:
        df.to_csv(csv_file, index=False, sep="\t")


if __name__ == "__main__":
    main()

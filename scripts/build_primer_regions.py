"""Build the primer regions table from a primer BED/BEDPE file.

Port of the upstream workflow/scripts/build_primer_regions.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
  build_primer_regions.py PRIMERS.bed --output PRIMER_REGIONS.tsv --log LOG

The upstream snakemake script API is replaced by positional input + flags;
the chunked parsing (chunksize = 10**6) and the output column layout
(chrom, left_start, left_end, right_start, right_end) are unchanged. BEDPE
inputs are converted to per-chromosome pairs (inter-chromosomal primer pairs
are reported to the log).
"""

import argparse


def parse_bed(input_path, log_file, out):
    print("chrom\tleft_start\tleft_end\tright_start\tright_end", file=out)
    for data_primers in pd.read_csv(
        input_path,
        sep="\t",
        header=None,
        chunksize=chunksize,
        usecols=[0, 1, 2, 5],
    ):
        for row in data_primers.iterrows():
            row_id = row[0]
            row = row[1]
            if row[5] == "+":
                print(
                    "{chrom}\t{start}\t{end}\t-1\t-1".format(
                        chrom=row[0], start=row[1] + 1, end=row[2]
                    ),
                    file=out,
                )
            elif row[5] == "-":
                print(
                    "{chrom}\t-1\t-1\t{start}\t{end}".format(
                        chrom=row[0], start=row[1] + 1, end=row[2]
                    ),
                    file=out,
                )
            else:
                print("Invalid strand in row {}".format(row_id), file=log_file)


def parse_bedpe(input_path, log_file, out):
    for data_primers in pd.read_csv(
        input_path,
        sep="\t",
        header=None,
        chunksize=chunksize,
        usecols=[0, 1, 2, 3, 4, 5],
    ):
        valid_primers = data_primers[0] == data_primers[3]
        valid_data = data_primers[valid_primers].copy()
        valid_data.iloc[:, [1, 4]] += 1
        valid_data.drop(columns=[3], inplace=True)
        valid_data.dropna(how="all", inplace=True)
        valid_data.to_csv(
            out,
            sep="\t",
            index=False,
            header=["chrom", "left_start", "left_end", "right_start", "right_end"],
        )
        print(
            data_primers[~valid_primers].to_csv(sep="\t", index=False, header=False),
            file=log_file,
        )


chunksize = 10**6


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primers_bed", help="Primer BED or BEDPE file")
    parser.add_argument("--output", required=True, help="Output primer regions TSV")
    parser.add_argument("--log", required=True, help="Log file (invalid strands / BEDPE pairs)")
    args = parser.parse_args()

    with open(args.output, "w") as out:
        with open(args.log, "w") as log_file:
            if args.primers_bed.endswith("bedpe"):
                parse_bedpe(args.primers_bed, log_file, out)
            else:
                parse_bed(args.primers_bed, log_file, out)


if __name__ == "__main__":
    import pandas as pd

    main()

"""Transform Ensembl GTF gene annotations into a BED file.

Port of the upstream workflow/scripts/transform_gene_regions.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
  transform_gene_regions.py ANNOTATION.gtf --output GENE_ANNOTATION.bed --log LOG

Only gene features on chromosomes 1-22, X and Y are kept; coordinates are
1-based as in the upstream script.
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotation", help="Ensembl GTF annotation file")
    parser.add_argument("--output", required=True, help="Output BED file")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    if args.log:
        sys.stderr = open(args.log, "w")

    chromosomes = list(map(str, range(1, 23))) + ["X", "Y"]

    with open(args.annotation, "r") as annotations, open(args.output, "w") as out_file:
        for line in annotations.readlines():
            if line.startswith("#"):
                continue
            line = line.split("\t")
            chromosome = line[0]
            feature = line[2]
            if feature != "gene" or chromosome not in chromosomes:
                continue
            start = str(int(line[3]) - 1)
            end = line[4]
            desc = dict([x.split(" ") for x in line[8].split("; ")])
            gene = desc.get("gene_name", desc["gene_id"])[1:-1]
            print("\t".join([chromosome, start, end, gene]), file=out_file)


if __name__ == "__main__":
    main()

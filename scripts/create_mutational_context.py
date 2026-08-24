"""Create the trinucleotide context and mutation-count tables.

Port of the upstream workflow/scripts/create_mutational_context.py (MIT,
snakemake-workflows dna-seq-varlociraptor v6.10.0), adapted to argv:
  create_mutational_context.py CALLS.bcf REF.fasta --sample-alias ALIAS
    --group GROUP --vafs 5,10,...,100 --output-context CONTEXT.tsv
    --output-counts COUNTS.tsv --log LOG

The upstream snakemake script API is replaced by positional inputs + flags
(the vaf thresholds are comma-joined, matching the rule's --vafs argument).
Diagnostics are written to the log file exactly as upstream (sys.stderr is
redirected to it). Upstream behavior kept verbatim: records with multiple
alternate alleles abort with exit code 1; non-SNV records are skipped with a
log message; a reference-base mismatch prints a log message and exits with
code 0 (an upstream quirk, preserved for fidelity).
"""

import argparse
import sys

import numpy as np  # noqa: F401  (upstream imports it; kept for parity)
import pandas as pd
from Bio import SeqIO
import pysam


def get_ref_triplet(ref_seq, variant_pos):
    left_idx = variant_pos - 1 if variant_pos > 0 else 0
    right_idx = variant_pos + 2
    return ref_seq[left_idx:right_idx]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bcf", help="Variant callset (BCF/VCF)")
    parser.add_argument("ref", help="Reference FASTA")
    parser.add_argument("--sample-alias", required=True, help="Sample column to read AF from")
    parser.add_argument("--group", required=True, help="Group name stamped into the context table")
    parser.add_argument("--vafs", required=True, help="Comma-separated VAF thresholds in percent")
    parser.add_argument("--output-context", required=True, help="Output context TSV")
    parser.add_argument("--output-counts", required=True, help="Output mutation counts TSV")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    sys.stderr = open(args.log, "w")

    reference = SeqIO.parse(args.ref, "fasta")
    sample_alias = args.sample_alias
    # TODO rather use a combination of sample alias and group
    sample_name = args.group
    bcf = pysam.VariantFile(args.bcf)
    current_chrom_id = None
    current_chrom_seq = None
    single_base_substitutions = []
    mutation_counts = []

    for bcf_record in bcf:
        variant_chrom = bcf_record.chrom
        variant_pos = bcf_record.pos - 1
        ref_base = bcf_record.ref
        alt_bases = bcf_record.alts
        # alternatively check ANN field for VARIANT_CLASS == SNV?
        if len(alt_bases) > 1:
            print("Record has mutliple alterations", file=sys.stderr)
            print(f"{variant_chrom}\t{variant_pos}", file=sys.stderr)
            exit(1)
        if len(ref_base) != 1 or len(alt_bases[0]) != 1:
            print(
                f"Record skipped - No SNV: {variant_chrom}\t{variant_pos}",
                file=sys.stderr,
            )
            continue
        while variant_chrom != current_chrom_id:
            current_chrom = next(reference)
            current_chrom_id = current_chrom.id
            current_chrom_seq = current_chrom.seq
        ref_triplet = get_ref_triplet(current_chrom_seq, variant_pos)
        if ref_triplet[1] != ref_base:
            print(
                f"Error: Missmatching reference base: {variant_chrom}\t{variant_pos}\t{ref_base}\t{ref_triplet}",
                file=sys.stderr,
            )
            exit()
        allele_frequency = float(bcf_record.samples[sample_alias]["AF"][0])
        single_base_substitutions.append((ref_triplet, alt_bases[0], allele_frequency))

    df = pd.DataFrame(single_base_substitutions, columns=["Triplet", "Alt", "AF"])
    df["Sample"] = sample_name

    df.to_csv(args.output_context, sep="\t", index=False)

    # Count mutations
    for min_vaf in [int(v) for v in args.vafs.split(",")]:
        min_vaf = min_vaf / 100
        temp_df = df[df["AF"] >= min_vaf]
        mutation_counts.append((min_vaf, len(temp_df.index)))

    mutation_count_df = pd.DataFrame(
        mutation_counts, columns=["Minimum VAF", "Mutation Count"]
    )
    mutation_count_df.to_csv(args.output_counts, sep="\t", index=False)


if __name__ == "__main__":
    main()

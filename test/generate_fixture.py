#!/usr/bin/env python3
"""Regenerate the synthetic varlociraptor fixture reads from the real GRCh38
reference (first 100kb of chr21).

Why: the original hand-written reads (random sequences) mapped through vg
giraffe with MAPQ 0-11; GATK BaseRecalibrator's default
MappingQualityNotZeroReadFilter then dropped them and the BQSR table came
out empty (live: ApplyBQSR "The covariates table is missing ReadGroup
SRR702070 in RecalTable0"). Reads sampled from the actual reference map
with realistic mapping quality, which is what the BQSR rules need.

10k pairs x 75bp, ~1% error, 500bp insert, seed 42 (deterministic).
Usage: python3 gen-varl-fixture.py [path-to-GRCh38-fasta]
"""
import gzip
import random
import sys

N = 10000
LEN = 75
INSERT = 500
SEED = 42
FASTA = sys.argv[1] if len(sys.argv) > 1 else "resources/genome.dna.homo_sapiens.GRCh38.111.fasta"
OUT1 = "test/fixtures/raw/SRR702070_1.fastq.gz"
OUT2 = "test/fixtures/raw/SRR702070_2.fastq.gz"

seq = {}
name, buf = None, []
with open(FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                seq[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seq[name] = "".join(buf)

if "21" not in seq:
    raise SystemExit(f"contig 21 not in {FASTA}: {sorted(seq)[:5]}...")
# chr21[15_000_000:15_100_000]: sampled windows were tested live —
# 5.0-5.1Mb maps to a paralog at MAPQ 3 (repeat-rich p-arm), while
# 15.0-15.1Mb maps uniquely (192/200 probe reads at MAPQ 60).
window = seq["21"][15_000_000:15_100_000]
n = len(window)

rng = random.Random(SEED)
COMP = str.maketrans("ACGT", "TGCA")


def rc(x):
    return x.translate(COMP)[::-1]


def mutate(x):
    return "".join(b if rng.random() > 0.01 else rng.choice("ACGT") for b in x)


with gzip.open(OUT1, "wt") as f1, gzip.open(OUT2, "wt") as f2:
    for i in range(N):
        pos = rng.randrange(0, n - INSERT - LEN)
        left = mutate(window[pos:pos + LEN])
        right = mutate(rc(window[pos + INSERT:pos + INSERT + LEN]))
        q = f"SRR702070.{i + 1}"
        f1.write(f"@{q} 1 length={LEN}\n{left}\n+\n{'I' * LEN}\n")
        f2.write(f"@{q} 2 length={LEN}\n{right}\n+\n{'I' * LEN}\n")

print(f"wrote {N} pairs x {LEN}bp to {OUT1} / {OUT2}")

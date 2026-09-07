#!/usr/bin/env bash
# Acceptance test for the oxo-flow-varlociraptor port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
# oxo-flow v0.11.0 prints the plan to stderr; capture both streams
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards} ({log} stays literal)"
"$OXO" debug main.oxoflow 2>&1 | grep -E '\{(config\.|sample\}|group\}|input\[|output\[|threads\}|memory\})' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "==> UMI branch: dry-run with umi_read set"
# Exclusive-gate sanity: setting umi_read must activate the annotate_umis flow
# (bam_index_sorted -> umi_tools group -> mark_duplicates with --BARCODE_TAG BX)
# instead of the plain dedup; the BQSR chain keeps the DNA path.
sed -e 's/^umi_read = ""$/umi_read = "some-umi"/' \
    main.oxoflow > .umi-test-tmp.oxoflow
grep -q '^umi_read = "some-umi"$' .umi-test-tmp.oxoflow
trap 'rm -f .umi-test-tmp.oxoflow' EXIT
"$OXO" dry-run .umi-test-tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-umi-$$.txt 2>&1
# Plan lines instantiate wildcards, so match rule-name prefixes and require
# [run state for the positive checks. The two-space anchor before "[run"
# keeps `mark_duplicates` from matching its `mark_duplicates_umi` twin.
grep -qE "^  [0-9]+\. mapping::bam_index_sorted[^ ]*  \[run" /tmp/oxo-dryrun-umi-$$.txt \
    || { echo "UMI branch: bam_index_sorted not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. mapping::annotate_umis[^ ]*  \[run" /tmp/oxo-dryrun-umi-$$.txt \
    || { echo "UMI branch: annotate_umis not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. mapping::mark_duplicates_umi[^ ]*  \[run" /tmp/oxo-dryrun-umi-$$.txt \
    || { echo "UMI branch: mark_duplicates_umi not scheduled"; exit 1; }
if grep -qE "^  [0-9]+\. mapping::mark_duplicates  \[run" /tmp/oxo-dryrun-umi-$$.txt; then
    echo "UMI branch: plain mark_duplicates unexpectedly scheduled"; exit 1
fi
grep -qE "^  [0-9]+\. mapping::recalibrate_base_qualities[^ ]*  \[run" /tmp/oxo-dryrun-umi-$$.txt \
    || { echo "UMI branch: DNA recalibrate_base_qualities not scheduled"; exit 1; }
if grep -qE "^  [0-9]+\. mapping::recalibrate_base_qualities_rna[^ ]*  \[run" /tmp/oxo-dryrun-umi-$$.txt; then
    echo "UMI branch: rna recalibrate_base_qualities unexpectedly scheduled"; exit 1
fi
rm -f .umi-test-tmp.oxoflow
trap - EXIT
echo "  annotate_umis -> mark_duplicates_umi on; plain mark_duplicates off; DNA BQSR kept"

echo "==> RNA branch: dry-run with datatype=rna"
# get_recalibrate_quality_input RNA branch: SplitNCigarReads is inserted and
# the BQSR chain re-sources from results/split/{sample}.bam; the DNA twins
# and the consensus twins must stay off.
sed -e 's/^datatype = "dna"$/datatype = "rna"/' \
    main.oxoflow > .rna-test-tmp.oxoflow
grep -q '^datatype = "rna"$' .rna-test-tmp.oxoflow
trap 'rm -f .rna-test-tmp.oxoflow' EXIT
"$OXO" dry-run .rna-test-tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-rna-$$.txt 2>&1
grep -qE "^  [0-9]+\. mapping::splitncigarreads[^ ]*  \[run" /tmp/oxo-dryrun-rna-$$.txt \
    || { echo "RNA branch: splitncigarreads not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. mapping::bam_index_split[^ ]*  \[run" /tmp/oxo-dryrun-rna-$$.txt \
    || { echo "RNA branch: bam_index_split not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. mapping::recalibrate_base_qualities_rna[^ ]*  \[run" /tmp/oxo-dryrun-rna-$$.txt \
    || { echo "RNA branch: recalibrate_base_qualities_rna not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. mapping::apply_bqsr_rna[^ ]*  \[run" /tmp/oxo-dryrun-rna-$$.txt \
    || { echo "RNA branch: apply_bqsr_rna not scheduled"; exit 1; }
if grep -qE "^  [0-9]+\. mapping::recalibrate_base_qualities  \[run" /tmp/oxo-dryrun-rna-$$.txt; then
    echo "RNA branch: DNA recalibrate_base_qualities unexpectedly scheduled"; exit 1
fi
if grep -qE "^  [0-9]+\. mapping::apply_bqsr  \[run" /tmp/oxo-dryrun-rna-$$.txt; then
    echo "RNA branch: DNA apply_bqsr unexpectedly scheduled"; exit 1
fi
if grep -qE "^  [0-9]+\. consensus::(recalibrate_base_qualities|apply_bqsr)_consensus[^ ]*  \[run" /tmp/oxo-dryrun-rna-$$.txt; then
    echo "RNA branch: consensus BQSR twins unexpectedly scheduled"; exit 1
fi
rm -f .rna-test-tmp.oxoflow
trap - EXIT
echo "  splitncigarreads -> *_rna BQSR on; DNA and consensus twins off"

echo "==> testcase branch: dry-run with testcase_activate + testcase_locus"
# Debug module: both ported callers must schedule their
# gather_observations/testcase pair; nothing consumes results/testcases/
# (opt-in like upstream).
sed -e 's/^testcase_activate = false$/testcase_activate = true/' \
    -e 's/^testcase_locus = ""$/testcase_locus = "chr1:1000-2000"/' \
    main.oxoflow > .testcase-test-tmp.oxoflow
grep -q '^testcase_activate = true$' .testcase-test-tmp.oxoflow
grep -q '^testcase_locus = "chr1:1000-2000"$' .testcase-test-tmp.oxoflow
trap 'rm -f .testcase-test-tmp.oxoflow' EXIT
"$OXO" dry-run .testcase-test-tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-testcase-$$.txt 2>&1
grep -qE "^  [0-9]+\. testcase::gather_observations_[^ ]*  \[run" /tmp/oxo-dryrun-testcase-$$.txt \
    || { echo "testcase branch: gather_observations rules not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. testcase::testcase_[^ ]*  \[run" /tmp/oxo-dryrun-testcase-$$.txt \
    || { echo "testcase branch: testcase rules not scheduled"; exit 1; }
[ "$(grep -cE "^  [0-9]+\. testcase::(gather_observations|testcase)_[^ ]*  \[run" /tmp/oxo-dryrun-testcase-$$.txt)" -eq 4 ] \
    || { echo "testcase branch: expected 4 testcase rules scheduled"; exit 1; }
rm -f .testcase-test-tmp.oxoflow
trap - EXIT
echo "  gather_observations_* + testcase_* scheduled for both callers"

echo "==> fusions branch: dry-run with calling=fusions + fusion_activate"
# Fusions-continuation gate: the calling module re-wires to consume arriba
# candidates (convert_fusions_to_vcf) and the arriba preprocess/call twins
# join the plan; the fusion:: STAR/arriba family activates.
sed -e 's/^calling = "variants"$/calling = "fusions"/' \
    -e 's/^fusion_activate = false$/fusion_activate = true/' \
    main.oxoflow > .fusions-test-tmp.oxoflow
grep -q '^calling = "fusions"$' .fusions-test-tmp.oxoflow
grep -q '^fusion_activate = true$' .fusions-test-tmp.oxoflow
trap 'rm -f .fusions-test-tmp.oxoflow' EXIT
"$OXO" dry-run .fusions-test-tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-fusions-$$.txt 2>&1
grep -qE "^  [0-9]+\. fusion::arriba[^ ]*  \[run" /tmp/oxo-dryrun-fusions-$$.txt \
    || { echo "fusions branch: fusion::arriba not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. fusion::convert_fusions[^ ]*  \[run" /tmp/oxo-dryrun-fusions-$$.txt \
    || { echo "fusions branch: fusion::convert_fusions not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. calling::varlociraptor_preprocess_arriba[^ ]*  \[run" /tmp/oxo-dryrun-fusions-$$.txt \
    || { echo "fusions branch: varlociraptor_preprocess_arriba not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. calling::varlociraptor_call_arriba[^ ]*  \[run" /tmp/oxo-dryrun-fusions-$$.txt \
    || { echo "fusions branch: varlociraptor_call_arriba not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. filtering::merge_calls_fusions[^ ]*  \[run" /tmp/oxo-dryrun-fusions-$$.txt \
    || { echo "fusions branch: merge_calls_fusions not scheduled"; exit 1; }
if grep -qE "^  [0-9]+\. fusion::star_index[^ ]*  \[skip" /tmp/oxo-dryrun-fusions-$$.txt; then
    echo "fusions branch: fusion::star_index unexpectedly skipped"; exit 1
fi
rm -f .fusions-test-tmp.oxoflow
trap - EXIT
echo "  fusion:: family + arriba preprocess/call twins on; STAR index activated"

echo "PASS"

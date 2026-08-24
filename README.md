# oxo-flow-varlociraptor — Small and structural variant calling with Varlociraptor

[![CI](https://github.com/oxo-flow-community/oxo-flow-varlociraptor/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-varlociraptor/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`snakemake-workflows/dna-seq-varlociraptor`](https://github.com/snakemake-workflows/dna-seq-varlociraptor) @ `v6.10.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

Align paired-end short reads against the 1000 Genomes human pangenome with vg
giraffe, call candidate small and structural variants with freebayes and
delly, and get a scenario-driven somatic call set from Varlociraptor: alignment
properties are estimated, variants are called under a tumor scenario, FDR is
controlled per variant type (SNV/INS/DEL/MNV/BND/INV/DUP/REP), and the calls
are annotated with VEP and dbSNFP, filtered, and rendered as an interactive
datavzrd report with an oncoprint-style label-sorting table. You get
pangenome-aligned, deduplicated and recalibrated BAMs, per-gene coverage
tables, FDR-controlled and annotated variant calls (BCF/VCF), the
postprocessed variant table, and the variant and gene-coverage HTML reports.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda: `conda install -c bioconda oxo-flow-cli` (note: the
conda package may lag behind releases; other platform binaries are available
on the [releases page](https://github.com/Traitome/oxo-flow/releases)).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-varlociraptor.git
cd oxo-flow-varlociraptor
```

### 3. Requirements

- **Input reads (user-provided)** — paired-end FASTQ at
  `reads_dir/<sample>_1.fastq.gz` and `reads_dir/<sample>_2.fastq.gz`
  (paired-end only). The sample cohort is declared in `[[sample_groups]]` with
  group metadata (alias, platform, purity, datatype, calling mode) expressed
  through the sample group name — one sample group = one tumor sample with a
  fixed scenario. The repository ships tiny synthetic FASTQ fixtures
  (`test/fixtures/raw/`), so the default config validates and dry-runs
  cleanly without any downloads.
- **Reference data (none to prepare, or provide your own)** — by default all
  reference resources are downloaded by the workflow itself into `resources/`,
  exactly like upstream: the GRCh38 primary assembly FASTA (Ensembl release
  111) plus `.fai`/`.dict` indices, the Ensembl release 111 GTF annotation,
  the VEP cache and plugins (release 111), REVEL scores, the Ensembl
  known-variants VCFs, and the HPRC v1.1 human pangenome graph (~5 GB total;
  requires unimpeded network access to ensembl.org, zenodo.org and the AWS
  human-pangenomics bucket). To use pre-downloaded databases instead, place
  them at the paths the `ref::` rules declare (see `modules/ref.oxoflow`
  `output =` lines) and run with `skip_ref_downloads=true`.
- **Compute** — up to 96 CPUs / 32 GB per rule (freebayes candidate calling:
  96 threads; vg giraffe mapping: 64 threads; samtools sort: 16 threads / 32G;
  Varlociraptor calling: 8G).
- **Tools** — conda environments with pinned versions, one env per tool pin
  set under `envs/` (declared per rule via `[rules.environment]` in the
  module files). No containers are used; conda/mamba is required at runtime.
- **Disk** — the reference downloads under `resources/` are large (pangenome
  graph, VEP cache, known-variants VCFs), and `results/` grows with BAMs,
  BCFs, tables and reports.

## Usage

```bash
# Point OXO at your oxo-flow binary (>= 0.12.0)
export OXO=oxo-flow

# 1. Validate and lint the workflow
"$OXO" validate main.oxoflow
"$OXO" lint main.oxoflow

# 2. Preview the execution plan (the shipped fixtures are the default config)
"$OXO" dry-run main.oxoflow --samples first:1

# 3. Run with your own data: set reads_dir and the sample group list, then
"$OXO" run main.oxoflow -j 1

# Acceptance test (validate + lint + dry-run + debug):
bash test/run.sh
```

The port covers the **default-parameter execution path** for a single tumor
sample group (`SRR702070_group`, sample `SRR702070`, alias `tumor`, purity
1.0, calling mode `variants`): pangenome mapping, candidate calling,
Varlociraptor calling + FDR control over 8 variant types, VEP/dbSNFP
annotation, filtering, and the datavzrd report. The upstream branch modules
(trimming, primers, MAF export, population-db filtering, plugins, CHM1
benchmarking, bwa alignment, mutational burden/signatures, fusion calling)
are ported as gated modules, off by default — see "Gated branch modules".
`oxo-flow validate` counts every rule unconditionally: 135 rules / 259
dependencies (88 rules / 163 on the default path alone; the executor skips
gated rules at run time when their key is `false`).

### Configuration

| config key | upstream | default | notes |
|---|---|---|---|
| `reads_dir` | `config/units.tsv` | `test/fixtures/raw` | input FASTQs; upstream uses absolute /projects/... paths, the port reads from the repo fixtures |
| `skip_ref_downloads` | — | `false` | `true` when reference databases are pre-placed at the `ref::` resource paths |
| `trimming_activate` | `trimming` section | `false` | SRA download (fasterq-dump) + fastp trimming rules |
| `primers_activate` / `primers_fa1` / `primers_fa2` | `primers/trimming` | `false` / `""` / `""` | fgbio primer assignment/trimming flow; `primers_fa2` empty = single-end primer fasta |
| `maf_activate` | `maf/activate` | `false` | vcf2maf.pl MAF export (variants + fusions) |
| `population_db_activate` / `population_db_path` / `population_db_alias` / `population_db_fdr` / `population_db_events` | `population/db/*` | `false` / `resources/population_db.variants.bcf` / `tumor` / `0.05` / `somatic_tumor_high,somatic_tumor_medium` | population-db FDR control + db update (events comma-joined) |
| `bwa_align_activate` | linear-reference branch of `mapping.smk` | `false` | `map_reads_bwa` + `bwa_index` (the default path aligns with vg giraffe) |
| `mutational_burden_activate` / `mutational_burden_events` | `mutational_burden/*` | `false` / `somatic_tumor_low,somatic_tumor_medium,somatic_tumor_high` | covered-coding-sites count + burden curve/hist (events comma-joined) |
| `mutational_signatures_activate` | `mutational_signatures` section | `false` | COSMIC v3.4 signature fitting + plots |
| `benchmarking_activate` | `benchmarking.smk` | `false` | CHM1 benchmark rules that need no third-party downloads |
| `plugins_activate` / `cadd_build` / `cadd_version` / `cadd_variant_type` | `plugins.smk` / `wildcards` | `false` / `GRCh38` / `v1.7` / `snv` | CADD score download for VEP |
| `fusion_activate` | `fusion_calling.smk` (star_arriba meta wrapper) | `false` | STAR + Arriba fusion candidate calling |

The remaining upstream parameters of the ported path are pinned to the
upstream default values (see `config/scenario.yaml` and the module files):
event id `some_id` (variants + fusions types, events present /
somatic_tumor_high / somatic_tumor_medium, FDR threshold 0.05, mode
`local-smart`), VEP cache/plugins directories `resources/vep/{cache,plugins}`,
REVEL score file `resources/revel_scores.tsv.gz`, dbSNFP
`resources/variation.vcf.gz`, and the datavzrd report templates.

### Gated branch modules

Nine upstream branch modules are ported as separate files under `modules/`,
each gated by `when = "config.<key>_activate"` (default `false`, so the
default path is unchanged). Each module header documents its rule map,
frozen wildcards, deviations, and excluded upstream rules:

| module | upstream rules | notes |
|---|---|---|
| `modules/trimming.oxoflow` | SRA download + fastp | no `envs/sra-tools.yaml`/`envs/fastp.yaml` in the repo — base env |
| `modules/primers.oxoflow` | fgbio primer flow | primer panel frozen to `uniform`; needs `bwa_index` |
| `modules/maf.oxoflow` | vcf2maf.pl export | fusions pair has no upstream producer in the port (no fusions calling) |
| `modules/population.oxoflow` | population-db filter/update + `gather_annotated_calls` | db read as-is, no input edge (upstream `before_update` flag) |
| `modules/plugins.oxoflow` | CADD download | REVEL rules already ported in `ref.oxoflow` |
| `modules/benchmarking.oxoflow` | CHM1 benchmark (3 rules) | `chm_eval_sample`/`chm_eval_kit`/`chm_eval` excluded (13 GB + third-party downloads) |
| `modules/mapping_bwa.oxoflow` | bwa index + align | no `envs/bwa.yaml` in the repo — base env |
| `modules/burden_signatures.oxoflow` | mutational burden + signatures | 20 VAF thresholds (5..100 step 5); `gather_annotated_calls` feeds the burden input |
| `modules/fusion.oxoflow` | star_arriba candidate calling | ends at the group candidate BCF; no `envs/star.yaml` in the repo — base env |

## Source

Upstream: **[snakemake-workflows/dna-seq-varlociraptor](https://github.com/snakemake-workflows/dna-seq-varlociraptor)** @ `v6.10.0` (commit `b65c3350b31d68f7fd36a497b9d11b37f2d03df3`), MIT license. Created 2026-08-15; this workflow may lag behind upstream releases. See [NOTICE.md](NOTICE.md) for attribution.

## Fidelity

The port aims for byte-identical commands on the default path. Known,
deliberate deviations:

| upstream | port | reason |
|---|---|---|
| `scatter.calling(16)` (rules run 16x, once per scatter item) | single chunk, `scatteritem=0` | oxo-flow has no scatter construct; with one small sample the 16 chunks are identical work |
| rule outputs that are directories (VEP cache/plugins, oncoprint `label_sortings/`/`variant-oncoprints/` dirs) | directory + `.completed` marker file output | oxo-flow targets files, not directories |
| scenario rendered at run time from `config/scenario.yaml` (yte template) | pre-rendered `resources/scenarios/SRR702070_group.yaml` for the default sample group; the template is kept verbatim at `config/scenario.yaml` | one scenario (purity 1.0) in the default path |
| `download_vep_plugins.py` with a hard-coded Ensembl variation FTP list and fallback | the `--release`/`--output`/`--log` argv variant of the same wrapper port | one release (111), one output dir; the FTP fallback list was dropped as dead code in the default path |
| wrapper-utils based rules (calls, tables, report) | plain `python scripts/*.py` argv ports of the same wrappers | wrapper-utils is a Snakemake runtime; the ported scripts keep the wrapper logic verbatim |
| `filter_odds` | not ported | not reachable in the default path (`filter: present` only); the population/burden branches consume `gather_annotated_calls` instead (ported in `population.oxoflow`) |
| template oncoprint views (`gene_oncoprint` / `variant_oncoprints` datasets) | empty (upstream defaults with a single group) | `prepare_oncoprint` itself runs and feeds the label-sorting table, exactly like upstream |
| vembrane filter/table expressions evaluated from Python at run time | precomputed literal expression/header (34 columns) | same semantics, evaluated once |
| upstream `config/units.tsv` absolute `/projects/...` read paths | `config.reads_dir` + sample group fixture paths | portability |
| Snakemake `temp()` outputs | `temporary = true` | engine equivalent |
| per-rule conda environments | one env per tool pin set (`envs/`) | same packages, same pins, consolidated |
| snakemake `before_update`/`update` flags (population db) | no input edge; the db path is read/written as-is | oxo-flow has no such flags; a declared input would create a DAG cycle (`validate` rejects it) |
| snakemake `temp()` outputs of the gated branch modules | plain outputs (`temporary = true` where the default path used it) | see the module headers; `join_mutational_signatures` writes with `>` instead of the upstream `>>` because the engine does not pre-delete outputs |
| snakemake script API (`snakemake.input/output/params`) in the 6 branch scripts | argv ports (`--output`/`--log` flags, comma-joined lists) | same logic verbatim, cf. the default-path script ports |
| `chm_eval_sample`/`chm_namesort`/`chm_to_fastq`/`chm_eval_kit`/`chm_eval` (benchmarking) | excluded | 13 GB CHM1 alignment + third-party CHM-eval kit downloads; the ported rules keep the orphan inputs that would only be produced by them |

## Test

```bash
bash test/run.sh
```

Runs `oxo-flow validate`, `oxo-flow lint`, a `dry-run` smoke check and a debug
scan for unexpanded wildcards; CI runs the same script on every push.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md); the upstream MIT license is included verbatim at
[LICENSE.upstream](LICENSE.upstream).

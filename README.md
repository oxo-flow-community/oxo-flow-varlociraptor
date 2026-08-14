# oxo-flow-varlociraptor

Small and structural variant calling with Varlociraptor.
This repository is a port of
[snakemake-workflows/dna-seq-varlociraptor](https://github.com/snakemake-workflows/dna-seq-varlociraptor)
v6.10.0 (MIT) to [oxo-flow](https://github.com/Traitome/oxo-flow) (>= 0.11.0).

## Quick start

```bash
# Point OXO at your oxo-flow binary (>= 0.11.0)
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

## About

The upstream workflow aligns short reads against the 1000 Genomes human
pangenome with vg giraffe, calls candidate variants with freebayes and delly,
estimates alignment properties, calls variants with Varlociraptor under a
scenario, performs FDR control per variant type, annotates with VEP and
dbSNFP, filters, and renders a datavzrd report with an oncoprint-style
label-sorting table.

The port covers the **default-parameter execution path** for a single tumor
sample group (`SRR702070_group`, sample `SRR702070`, alias `tumor`, purity
1.0, calling mode `variants`): pangenome mapping, candidate calling,
Varlociraptor calling + FDR control over 8 variant types, VEP/dbSNFP
annotation, filtering, and the datavzrd report. 88 rules in total.

### Inputs

- Reads: `reads_dir/<sample>_1.fastq.gz` and `reads_dir/<sample>_2.fastq.gz`
  (paired-end only). The sample cohort is declared in `[[sample_groups]]` with
  group metadata (alias, platform, purity, datatype, calling mode) expressed
  through the sample group name — the port follows the upstream convention
  that one sample group = one tumor sample with a fixed scenario.
- All reference resources (1000 Genomes pangenome GBZ, GRCh38 FASTA, VEP
  cache/plugins, REVEL scores, dbSNFP/dbSNP VCFs, annotation GTF) are
  downloaded by the workflow itself into `resources/`, exactly like upstream.

The repository ships tiny synthetic FASTQ fixtures (`test/fixtures/raw/`), so
the default config validates and dry-runs cleanly without any downloads.

### Configuration

| config key | upstream | default | notes |
|---|---|---|---|
| `reads_dir` | `config/units.tsv` | `test/fixtures/raw` | input FASTQs; upstream uses absolute /projects/... paths, the port reads from the repo fixtures |

The remaining upstream parameters of the ported path are pinned to the
upstream default values (see `config/scenario.yaml` and the module files):
event id `some_id` (variants + fusions types, events present /
somatic_tumor_high / somatic_tumor_medium, FDR threshold 0.05, mode
`local-smart`), VEP cache/plugins directories `resources/vep/{cache,plugins}`,
REVEL score file `resources/revel_scores.tsv.gz`, dbSNFP
`resources/variation.vcf.gz`, and the datavzrd report templates.

## Outputs

All paths are relative to the run directory and follow the upstream layout
(`results/...`), so rendered commands are byte-identical to the upstream
snakemake-wrappers ports.

- `results/calls/...` — Varlociraptor raw calls, VEP-annotated calls
  (`results/calls/vep_annotated/`), dbSNFP-annotated calls
  (`results/calls/db_annotated/`), annotation-filtered calls
  (`results/calls/filtered/`), per-variant-type FDR-controlled calls
  (`results/calls/fdr-controlled/{group}/some_id/*.variants.bcf`) and the
  merged final call set
  (`results/final-calls/{group}/{group}.some_id.variants.fdr-controlled.bcf`).
- `results/tables/...` — the postprocessed variant table
  (`{group}.some_id.variants.postprocessed.fdr-controlled.tsv`) and the
  label-sorting outputs for the oncoprint table.
- `results/datavzrd-report/all.some_id.variants.fdr-controlled/index.html` —
  the interactive variant report (datavzrd).
- `results/datavzrd-report/{group}.coverage/index.html` and
  `results/coverage/{group}.csv` — the gene-coverage report and table.
- `results/mapped/`, `results/dedup/`, `results/recal/`, `results/qc/`,
  `results/regions/`, `results/observations/`,
  `results/candidate-calls/`, `results/scenarios/` — intermediate artifacts
  mirroring upstream.

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
| `gather_annotated_calls` / `filter_odds` | not ported | not reachable in the default path (benchmarking off + `filter: present` only) |
| template oncoprint views (`gene_oncoprint` / `variant_oncoprints` datasets) | empty (upstream defaults with a single group) | `prepare_oncoprint` itself runs and feeds the label-sorting table, exactly like upstream |
| vembrane filter/table expressions evaluated from Python at run time | precomputed literal expression/header (34 columns) | same semantics, evaluated once |
| upstream `config/units.tsv` absolute `/projects/...` read paths | `config.reads_dir` + sample group fixture paths | portability |
| Snakemake `temp()` outputs | `temporary = true` | engine equivalent |
| per-rule conda environments | one env per tool pin set (`envs/`) | same packages, same pins, consolidated |

## Not ported (metadata `excluded`)

Everything outside the default-parameter main path of the upstream workflow:
bwa mapping, read trimming (fastp), fusion calling (arriba), MAF conversion,
mutational burden and mutational signature analyses, the population database
(germline AF annotation), dgidb druggability, CADD annotation, primer design,
benchmarking, consensus reads, target regions, and the template oncoprint
views. None of these are reachable with the upstream default configuration.

## License

Apache-2.0 (this port), see `LICENSE` and `NOTICE.md`.
The upstream workflow is MIT; its license text is kept verbatim at
`LICENSE.upstream` (Apache-2.0 §4(d): attribution notices from the Source
form must be retained).

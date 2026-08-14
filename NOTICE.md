oxo-flow-varlociraptor
Copyright (c) 2026 oxo-flow-community

This pipeline is a port of snakemake-workflows/dna-seq-varlociraptor
(https://github.com/snakemake-workflows/dna-seq-varlociraptor), version
6.10.0 (commit b65c3350b31d68f7fd36a497b9d11b37f2d03df3), authored by the
dna-seq-varlociraptor workflow developers (Felix Mölder, Johannes Köster,
and contributors).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---------------------------------------------------------------------
Upstream license

This port is derived from snakemake-workflows/dna-seq-varlociraptor under
the MIT license. The upstream LICENSE is included verbatim in this
repository at LICENSE.upstream (fetched from the upstream repository at
the ported commit b65c3350b31d68f7fd36a497b9d11b37f2d03df3, tag v6.10.0).
(Apache-2.0 §4(d): attribution notices from the Source form must be
retained.)
---------------------------------------------------------------------

Files copied or adapted from the upstream repository (MIT, kept as close to
the originals as the engine allows):

- config/scenario.yaml — verbatim copy of the upstream yte scenario template
- config/super_interesting_genes.tsv — verbatim copy of the upstream default
- resources/datavzrd/variant-calls-template.datavzrd.yaml,
  gene-coverage-template.datavzrd.yaml, linkouts.js — verbatim copies of
  workflow/resources/datavzrd/* at the ported commit
- resources/scenarios/SRR702070_group.yaml — scenario rendered from
  config/scenario.yaml for the default sample group (purity 1.0)
- test/fixtures/raw/*.fastq.gz — synthetic minimal FASTQs (not from
  upstream; generated for this port so the default config dry-runs cleanly)

Scripts in scripts/ fall into two groups:

1. Adapted from the upstream repository's own workflow/scripts/ (MIT),
   converted from Snakemake wildcard/param references to plain argv
   interfaces with the same logic and outputs:
   - scripts/oncoprint.py (upstream workflow/scripts/oncoprint.py)
   - scripts/process-call-tables.py (upstream
     workflow/scripts/process-call-tables.py; the population-db and
     join_short_obs branches were dropped as unreachable in the default
     path)
   - scripts/coverage_table.py (upstream workflow/scripts/coverage_table.py)
   - scripts/transform_gene_regions.py (upstream
     workflow/scripts/transform_gene_regions.py)

2. Reimplementations of snakemake-wrappers
   (https://github.com/snakemake/snakemake-wrappers, MIT) invoked by the
   upstream rules, as plain argv scripts, keeping the commands and tool
   versions identical:
   - scripts/download_vep_plugins.py — port of the bio/vep/plugins wrapper
     (v8.0.0); the inline get_vep_cache shell mirrors the bio/vep/cache
     wrapper (v8.0.0)
   - scripts/datavzrd.py — port of the utils/datavzrd wrapper (v9.14.0)

Conda environment files in envs/ pin the same packages (and versions) as the
upstream workflow's envs and the wrapper environments of the invoked
snakemake-wrappers versions.

# Reproduce GSE84161 normalization and mapping

Run from the repository root. This workflow prepares expression measurements and descriptive assay QC. It does not calculate program scores or treatment/donor contrasts. The effective [execution plan](../config/gse84161-normalization-plan.json) fixes the input hashes, sample order, software and RMA settings. Earlier execution plans are preserved when implementation or runtime repairs require a new freeze.

## Inputs and runtime

The large inputs remain in ignored local storage. The [source manifest](../config/gse84161-gene-map-sources.json) records the complete GPL570 SOFT snapshot, original Stankey program ZIP, Ensembl 116 crosswalk and independent BioMart count response. The [raw capture](../config/gse84161-raw-capture.json) and [array inventory](../data/derived/gse84161-raw-arrays.csv) identify the 30 CEL files and their compressed/decompressed hashes. Recover the exact pinned bytes before replay; a changed public response must be reviewed separately instead of accepting its new hash automatically.

The Ensembl export uses the June 2026 archive, `hsapiens_gene_ensembl`, no filters, and the three attributes `ensembl_gene_id`, `entrezgene_id`, `hgnc_symbol`. Both exact XML queries are retained in the source manifest. Require the terminal `[success]` marker, 91,748 relationship rows, and 86,411 distinct Ensembl IDs matching the separately returned count. Do not use a live or target-only lookup to replace this full relation graph. [Official Ensembl BioMart documentation](https://jun2026.archive.ensembl.org/info/data/biomart/index.html?redirect=no).

Follow the separate [runtime reproduction instructions](gse84161-runtime-reproduction.md). The runtime lives in `work/gse84161-normalization/runtime/`; the older coloc runtime is preserved. The base package lock and any documented same-version build override together describe the actual executed environment. A package version alone is insufficient to reproduce a local build override.

Extract exactly the 30 expected `.CEL.gz` members into `work/gse84161-normalization/cel/`. Preserve their compressed names and bytes. The runner checks every extracted member against its saved sample accession, byte count and SHA-256 before opening intensities. It also verifies the complete archive and rejects missing, duplicate or unexpected members.

## Replay

With the pinned inputs and runtime present:

```bash
.venv/bin/python scripts/map_gse84161.py
.venv/bin/python scripts/verify_gse84161_mapping.py
work/gse84161-normalization/runtime/bin/Rscript scripts/normalize_gse84161.R --plan config/gse84161-normalization-plan.json --check-only
work/gse84161-normalization/runtime/bin/Rscript scripts/normalize_gse84161.R --plan config/gse84161-normalization-plan.json --execute
work/gse84161-normalization/runtime/bin/Rscript scripts/verify_gse84161_normalization.R
.venv/bin/python -m unittest discover -s tests -v
```

Check-only validates provenance and ordering without reading CEL intensities. Execute requires the frozen plan, reads all 30 arrays together, and computes full-CDF `affy::rma` with background correction, quantile normalization and `bgversion=2`. There is no additional log transform or normalization after per-gene median aggregation. [Official affy reference](https://www.bioconductor.org/packages/release/bioc/manuals/affy/man/affy.pdf).

Review execution success before proceeding to the independent verifier. Preserve previous local outputs if a separate replay record is required: a normal replay writes the same output paths. Mapping tables use deterministic gzip serialization. Time-stamped audit records and PDF metadata can change across otherwise numerically identical replays; compare their substantive contents and numeric matrices rather than asserting universal byte identity.

## Outputs and interpretation

Full matrices are RDS numeric-double objects in ignored `work/gse84161-normalization/output/` storage: 54,675 probe sets by 30 arrays, and 20,486 Entrez genes by 30 arrays. The probe matrix preserves the full platform; the Entrez matrix retains every eligible array gene even if its Ensembl crosswalk is absent or ambiguous. The latter uses the median of every eligible probe for that gene, with explicit probe-ID joins and fixed sample order. RDS preserves double precision without rounded text export.

Versioned outputs include the complete probe map, gene crosswalk, per-member mapping/reasons, all nine program coverage rows, contributing probe counts, per-array raw/RMA/RLE summaries, diagnostic plots and independent verification records. Missing annotations never become zero measurements. Program coverage uses the original source-member denominators. The nonoverlap sensitivity uses its parent's coverage gate and a separate minimum on remaining genes.

The full vector diagnostic PDFs contain millions of plotted points and stay ignored at their original recorded paths. Their audit hashes are preserved. Compact PNG renders of those exact PDFs are versioned for review; this storage choice does not modify expression values, the QC summaries or point selection.

Raw PM summaries describe CDF-defined perfect-match entries. RLE subtracts each probe's across-array median for diagnostics only; it does not replace the expression matrix. No signal-based exclusions are automatic. Passing numeric and identity checks does not establish donor correspondence, remove batch effects or demonstrate a biological treatment response.

The fixed method was recorded before new program coverage and array processing. Program coverage was known before the RMA execution freeze. Later source-manifest wording clarified missing-ID reason labels and the nonoverlap denominator without changing source bytes or mapped membership; the effective plan records that chronology. The [earlier study investigation](../reports/gse84161-study-investigation.md) and both validation protocols remain historical snapshots, not descriptions of this later completed preparation.

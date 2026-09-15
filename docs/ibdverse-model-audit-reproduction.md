# Reproduce the IBDverse original-model audit

This stage audits the public inputs for the two previously selected blood-monocyte contexts and PRKD2. It does not fit a new association, infer new biological causal components, request model predictions or calculate a shared-cause posterior. The [report](../reports/ibdverse-model-audit.md) gives the findings and exact remaining inputs.

## Frozen scope and source boundary

- [Audit plan](../config/ibdverse-model-audit-plan.json): fixed gene, two contexts, original GRCh38 gene window, finite source families, qualification rules and 150-MB initial acquisition ceiling.
- [Baseline freeze](../config/ibdverse-model-audit-freeze.json): the plan hash, preceding commit and hashes of all 1,256 pre-existing tracked files.
- [Acquisition ledger](../config/ibdverse-model-audit-sources.json): source URLs, retrieval dates, bounded request parameters, pinned repository revisions, response identities, errors and successful-content hashes. BioStudies and PMC requests and skill-script hashes are recorded separately from ordinary file requests. Skill input JSON is retained under `config/ibdverse-model-audit-requests/`.
- [Four-member ledger](../config/ibdverse-model-audit-member-queries.json): exact compressed ranges and complete source CRC, size, SHA-256 and row-count results. The complete preceding 6,047-member ZIP directory is reused by hash; the 393-GB parent is not downloaded.
- [H5AD access plan](../config/ibdverse-model-audit-h5ad-plan.json), [schema-triggered mapping amendment](../config/ibdverse-model-audit-h5ad-mapping-amendment.json) and [range ledger](../config/ibdverse-model-audit-h5ad-queries.json): only observation metadata is decoded. The mapping amendment was frozen after the schema was inspected and before context counts were read.
- [Coloc range ledger](../config/ibdverse-model-audit-coloc-queries.json): complete three-member directory and bounded nested-gzip headers only. These are not complete result-stream scans or complete-member CRC checks.

Retained successful sources and auxiliary ranges total **41,343,821 bytes**. This counts retained bytes, not a precise network billing total: failure-response bodies and skill receipt overhead are not all counted. The HDF5 reader used 19 blocks of 262,144 bytes. A fetched block can contain neighboring bytes; no expression array was decoded or analyzed.

All full raw files, observation identifiers, sample/individual maps, rendered paper pages and HTTP caches remain in ignored `data/raw/ibdverse-model-audit/`. Versioned derived tables contain aggregate counts, source-association rows and provenance. The code from the study and its linked repositories was inspected as data, not executed.

## Offline reproduction

Run from the repository root using the existing `.venv`. The analysis uses Python, NumPy, pandas, h5py and requests already available in the project environment; no external model service is needed. Exact runtime versions are in the [stage verification](../reports/ibdverse-model-audit-stage-verification.json).

The local ignored cache is required for exact offline reproduction. A fresh Git clone alone does not contain the full raw source data or cached remote ranges. Restore these from a preserved cache and verify every recorded hash. For reacquisition, use the exact recorded URLs/ranges and source identity checks; a changed remote file must fail the frozen-byte checks rather than silently replace this snapshot. Mutable public metadata, repository-head responses and the expired anonymous browser-download link may require the preserved snapshot for exact historical replay. This guide does not claim a one-command fresh-download reconstruction of those historical services.

```bash
.venv/bin/python scripts/replay_ibdverse_model_audit.py
```

The replay performs these steps:

```bash
.venv/bin/python scripts/fetch_ibdverse_model_members.py --offline
.venv/bin/python scripts/read_ibdverse_model_metadata.py --offline
.venv/bin/python scripts/analyze_ibdverse_model_audit.py
.venv/bin/python scripts/verify_ibdverse_model_audit.py
```

It regenerates four exact PRKD2 source excerpts, nine audit CSVs and three JSON analysis/metadata reports. Those **16 generated files** must match their saved bytes. The additional coloc directory CSV is reconstructed independently from the cached central directory and must also match its frozen hash. All acquisition ledgers remain unchanged.

CSV output retains the standard writer's CRLF record endings. The original gene-summary TSV excerpts retain the trailing separator for the empty `qvals` field. A directory-scoped `.gitattributes` declares those two data-format conventions for Git's whitespace check; source fields and recorded hashes are not altered to remove meaningful separators.

The primary ZIP implementation checks the local filename, complete raw-DEFLATE boundary, uncompressed size and CRC. The independent verifier uses the standard ZIP reader to scan all four complete members and independently extract PRKD2 rows. It also validates all 6,047 central-directory entries. The HDF5 count verifier uses categorical-code tuple counts followed by the public mapping, independent of the primary expanded-column pandas joins and groupings. It confirms 2,196,874 mapped cells and the fixed two-context aggregates without exporting individual identifiers.

The verifier checks code-file Git blob hashes against four complete pinned repository trees, every retained source/receipt hash, source numeric strings and all differences between the earlier nominal exports and new lead summaries. It requires preservation of every preceding scientific file; the only permitted prior-file edits are insertions in README and the evidence log. RNA eligibility and conditional PC arithmetic must remain explicitly unverified as executed-model metadata.

## Repository checks

```bash
.venv/bin/python scripts/build_catalog.py
.venv/bin/python -m unittest discover -s tests -v
git diff --check
```

The catalog rebuild must leave its previous outputs unchanged. New failure-case tests cover truncated/appended/corrupt compressed members, malformed rows, exact gene identifiers, preserved blank/source values, missing categorical labels, duplicate/incomplete sample joins, bounded offline access and the distinction between inferred counts and a qualified new causal comparison.

The [independent verification](../reports/ibdverse-model-audit-verification.json), [offline replay](../reports/ibdverse-model-audit-replay-verification.json) and [stage verification](../reports/ibdverse-model-audit-stage-verification.json) record completed checks. Exact model N and selected PCs, executed-run provenance, signed study covariance or complete RNA components, and PSC covariance/sample-scale provenance remain missing. Reproduction of this audit does not resolve those scientific limitations.

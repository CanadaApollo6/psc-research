# PSC BACH2: complete filtered-cohort input qualification

Date: 16 September 2026 UTC. **Eight released archives; no expression effects.** This completes the authorized payload pass. The four [metadata](psc-bach2-input-qualification.md) and four [first-archive](psc-bach2-count-qualification.md) deliverables remain unchanged.

## Decision

**The complete released RNA-feature source-count-format interface is qualified for all 55,460 fixed retained cells across eight source donors.** All eight full ordered RNA ID/name/type axes are identical: **36,601 `Gene Expression` features**. All retained sample-plus-barcode keys match. No gene-axis intersection, union, renaming, collapse or imputation is needed or performed. No expression matrix assembly, normalization, score, gene/state selection or model fit was performed.

Every archive also contains the same **38 `Antibody Capture` features**, kept separate from RNA. **No deposited antibody feature explicitly names BACH2.** These are readable author-named marker labels, not opaque tag-only identifiers. The release therefore supplies no explicit direct BACH2-protein feature; RNA BACH2 is not an ADT measurement.

This qualifies a **released-source quantity and cell interface**, not an untouched-original-UMI or count-likelihood certificate. Every MEX file explicitly declares `cellranger-6.1.1`, while the source Methods/IDF describe Cellranger 3.0.2 and GRCh38 release 93. Repetition across eight archives does not resolve the export/reference history or supply missing historical RNA/SCT margins. Those limits remain, without triggering another acquisition search.

## Bounded acquisition and integrity

Only the remaining seven already inventoried filtered archives were downloaded: **298,421,791 new body bytes**, exactly the approved allowance. **Zero new metadata body bytes** were requested or read, out of the separate 1,000,000-byte allowance. Existing sample05 was reused without a body redownload. All eight compressed archives total **325,315,169 bytes**.

The seven downloads ran through the ordinary [BioStudies E-MTAB-14013 file route](https://www.ebi.ac.uk/biostudies/files/E-MTAB-14013/sample01_filtered_feature_bc_matrix.tar.gz), with one BioStudies → EBI FTP HTTPS redirect per file. Acquisition ran from 03:01:20 to 03:10:38 UTC on 16 September 2026. Fixed filenames, exact Content-Length/body caps, at most three workers, exclusive source-file creation and a shared stop-on-uncertainty event were used. No retry or unlisted archive request occurred. Receipts and the initial implementation snapshot are retained. The maintained script additionally restricts URL paths/ports and checks individual receipts against the acquisition summary; the actual recorded source routes satisfy these tighter checks.

| Source sample | Archive bytes | SHA-256 |
|---|---:|---|
| `sample01` | 54,449,959 | `cbede10441fecfc1013077088829ddb9916140d9f8489c371c98efb3492837a6` |
| `sample02` | 44,540,169 | `b47a0d5322ff6d560e385d5b205b4081b9c88dc56bcf3a7e702863d9a41c64e8` |
| `sample03` | 43,425,341 | `f0a957fadccf7998a8c366e537483d23bceba1b095b5cbf90e3920402dcf2091` |
| `sample04` | 36,719,700 | `0a7bc817823e121b929a1716abb1f32b132c909b19b0b159a17f2484393f7d52` |
| `sample05` | 26,893,378 | `5c57cc1c8eb48ce0208062b92fcd056b39e0f90e8eb6cb712ccd0d76f95dda29` |
| `sample06` | 36,674,386 | `e251f2df1148f9ae06c831360295782b681bece649395d2b2bbfe2aaf4b8e826` |
| `sample07` | 29,629,283 | `631ef4435687e966ce3fe743691bd35d601928460c096c1e42ff0701cdf9c9a9` |
| `sample08` | 52,982,953 | `1808db1c3333642526eedd0b880c98cd4df4adcdd5d62866f17dae6d56666e03` |

For each sample the archive is named `<source sample>_filtered_feature_bc_matrix.tar.gz`. Exact original URLs, resolved URLs and retrieval timestamps are in the [aggregate JSON](../data/derived/psc-bach2-cohort-qualification.json). Phase plan, annotation addendum, seven individual receipts and acquisition summary are under ignored `data/raw/psc-bach2-qualification/cohort/`. The sample05 body remains in its closed first-phase raw directory.

All **32 gzip streams**—eight outer archives and 24 nested members—pass complete CRC/ISIZE/EOF validation. Each TAR contains one correctly sample-prefixed directory and exactly three regular gzip members: `features.tsv.gz`, `barcodes.tsv.gz` and `matrix.mtx.gz`. Member checks cover checksums, sizes, paths, aliases, types, bounds and zero end blocks. There are no links, traversal paths, sparse/special files or unparsed nonzero appendages. An independent manual 512-byte TAR parser also checks every header, payload padding and terminal block. Original member paths are never used as output destinations. Every compressed and decoded member is hashed, with bounded output and no source overwrite.

## Full feature and source-cell interfaces

Every source matrix has **36,639 features**: 36,601 RNA and 38 antibody features. The complete ordered three-column axes match across all eight archives, not just their dimensions or feature-ID sets. IDs are unique. Ten duplicated source name values remain unchanged; no name-based merge occurs. The RNA IDs have unversioned human ENSG syntax, but no independent gene/reference reconciliation was attempted.

The canonical ordered RNA-axis SHA-256 is `aef32a0841aec54ce17a661c064501697b2b7885bb53dff8a930eef74230ae2d` (UTF-8 compact JSON of original ID/name/type triples). Original decoded feature-file hashes are separately recorded. Full axes, including implicit all-zero feature rows, remain intact. This is the complete **deposited** axis, not a reconstructed pre-filter reference universe.

| Source sample | Filtered source cells | Fixed retained cells, all matched | Source-only extras |
|---|---:|---:|---:|
| `sample01` | 9,876 | 9,327 | 549 |
| `sample02` | 7,596 | 7,034 | 562 |
| `sample03` | 8,240 | 7,771 | 469 |
| `sample04` | 6,864 | 6,477 | 387 |
| `sample05` | 4,703 | 4,278 | 425 |
| `sample06` | 6,519 | 6,090 | 429 |
| `sample07` | 5,299 | 4,933 | 366 |
| `sample08` | 10,218 | 9,550 | 668 |

Total: **59,315 source cells; 55,460 fixed retained cells matched; zero retained cells missing; 3,855 source-only extras.** Those extras are not given invented exclusion reasons, states or clinical identities. Their count does not reconcile the separate paper 56,209-cell versus annotation 55,460-cell discrepancy.

The source SDRF is read with repeated file columns preserved. The published filtered-archive basename gives the metadata sample key and its accession-local source donor; source genotype characteristic/factor, sex and retained metadata condition agree. There are eight unique source donors and the already qualified four/four source genotype labels. Detailed sample/donor/genotype and cell joins remain ignored, with no candidate clinical join.

Only the exact matching `sample_name + '_'` prefix is removed from retained metadata barcodes. Sample remains part of every cell key. **2,006 bare source barcode strings recur across archives**; 1,769 recur in the fixed retained metadata. A global bare-barcode join would be wrong. No state label or expression value changes membership.

## Exact sparse quantities; RNA and ADT stay separate

All eight headers are `matrix coordinate integer general`, with the exact source comment:

```text
%metadata_json: {"software_version": "cellranger-6.1.1", "format_version": 2}
```

All **82,812,719** observed sparse coordinates match their declarations and axes. Values are finite, positive exact integers. There are no duplicate coordinates, out-of-range indices, fractional values or explicit stored zeros. Numeric parsing uses integers or exact decimal/rational values before any float conversion. Nothing is rounded or repaired into a count.

| Original feature type | Positive source coordinates | All source cells: exact source units | Fixed retained cells: exact source units |
|---|---:|---:|---:|
| `Gene Expression` | 81,380,756 | 232,942,493 | 209,220,999 |
| `Antibody Capture` | 1,431,963 | 110,669,357 | 99,655,981 |

These are technical aggregates, not genotype/state gene contrasts. RNA and antibody totals are not combined. Individual type-specific margins remain ignored. Source units are not certified untouched RNA UMIs, and antibody tags are **not absolute protein molecules**.

## Antibody target annotation: no explicit BACH2 protein feature

All 38 rows below appear in **every archive**, with the same original ID, name and type. For this release each feature ID equals its feature name; the first column preserves both exact strings. Order follows the original antibody rows. The panel contains 21 CD-style names and 17 other readable source marker names. No opaque target label is present in this enumeration. This source-label assessment does not validate antibody clone identity or specificity, reinterpret TCR notation, infer targets from RNA, or test translation effects.

| Original feature ID **and** identical original feature name | Original feature type |
|---|---|
| `ADT_CD45RA` | `Antibody Capture` |
| `ADT_CD45RO` | `Antibody Capture` |
| `ADT_CD62L` | `Antibody Capture` |
| `ADT_CCR7` | `Antibody Capture` |
| `ADT_CD31` | `Antibody Capture` |
| `ADT_CD27` | `Antibody Capture` |
| `ADT_CD28` | `Antibody Capture` |
| `ADT_CD95` | `Antibody Capture` |
| `ADT_CD122` | `Antibody Capture` |
| `ADT_CXCR3` | `Antibody Capture` |
| `ADT_CD11a` | `Antibody Capture` |
| `ADT_CD49d` | `Antibody Capture` |
| `ADT_CD127` | `Antibody Capture` |
| `ADT_CD25` | `Antibody Capture` |
| `ADT_HLA-DR` | `Antibody Capture` |
| `ADT_GARP` | `Antibody Capture` |
| `ADT_CD39` | `Antibody Capture` |
| `ADT_CD73` | `Antibody Capture` |
| `ADT_CTLA-4` | `Antibody Capture` |
| `ADT_GITR` | `Antibody Capture` |
| `ADT_ICOS` | `Antibody Capture` |
| `ADT_CCR4` | `Antibody Capture` |
| `ADT_CCR6` | `Antibody Capture` |
| `ADT_CCR5` | `Antibody Capture` |
| `ADT_CXCR4` | `Antibody Capture` |
| `ADT_CXCR5` | `Antibody Capture` |
| `ADT_CD44` | `Antibody Capture` |
| `ADT_CD69` | `Antibody Capture` |
| `ADT_CD161` | `Antibody Capture` |
| `ADT_TCR_Va72` | `Antibody Capture` |
| `ADT_CD279` | `Antibody Capture` |
| `ADT_CD223` | `Antibody Capture` |
| `ADT_KLRG1` | `Antibody Capture` |
| `ADT_CD57` | `Antibody Capture` |
| `ADT_CD366` | `Antibody Capture` |
| `ADT_TIGIT` | `Antibody Capture` |
| `ADT_TCRgd` | `Antibody Capture` |
| `ADT_TCR_Va24-Ja18` | `Antibody Capture` |

**None is BACH2.** Absence of a named BACH2 antibody feature is not a claim that BACH2 protein is absent biologically. No protein value or RNA–protein correlation was analyzed.

## Source-quantity limits and the next scientific boundary

The explicit Cell Ranger exporter metadata is direct source evidence beyond a filename or integer test. Together with the declared RNA feature type, complete source matrix and exact retained-cell interface, it supports the released RNA-feature source-count-format target above.

The Methods/IDF instead name **Cellranger 3.0.2** and **GRCh38 human reference genome (release 93)**. No executed reference bundle, export command or separate layer object is supplied in these three-file payloads. The historical gene/QC/Scrublet/SCTransform workflow remains a source description, not an instruction to force its filters onto this release. The fixed metadata has no RNA/SCT count or feature margins. Repeating the same MEX schema neither corrects those source versions nor creates missing export history. No untouched-original-UMI or count-likelihood gate is promoted as passed.

**No further payload search is started or needed for this completed qualification target.** A later donor-level descriptive question could be defined on the full 36,601-feature released RNA universe and fixed 55,460 retained cells. Its gene/state definitions, transformations, donor-level estimand and uncertainty checks still need an explicit freeze and authorization before effects. No such method was selected here. Direct BACH2-protein/translation measurement is not supplied by this named ADT panel. The earlier biological, clinical-key and variant-identification limits are unchanged.

## Validation and reproducibility

- **32 cohort-specific native synthetic tests pass.** They cover seven-file/byte restrictions, sample05 exclusion from download scope, no retries or overwrite, stop behavior, safe redirects, exact sample/donor/barcode joins, repeated source columns, axis/order/name differences, no implicit feature reconciliation, RNA/ADT separation, exact values beyond float precision, fractional values and source-label-only BACH2 annotation.
- The closed safe parser is imported without edits. Its 30 tests and the 30 preceding metadata tests also pass in native discovery.
- **18 independent synthetic tests pass.** The independent implementation imports neither main qualification parser. Two complete independent source scans and all 28 ordered pairwise RNA-axis comparisons pass; all ADT axes are also identical. Exact main/independent comparisons pass: **696/696 aggregates; 1,172,448 full ID/name/type/position fields; 177,945 source barcode/position/retention fields; 118,630 RNA/ADT column margins; 38 antibody triples and 304 sample presences**. Independent positive-feature-per-cell margins have no main counterpart and are not included in that comparison count. The independent implementation and comparison records are under `work/psc-bach2-cohort-qualification/independent/`.
- The final main script reproduces the aggregate JSON **byte for byte offline**. The two complete independent passes also reproduce 33 full axis/join/margin/difference files byte for byte. Combined native BACH2 discovery runs **92 tests, all passing**.
- Full native repository discovery: **837 tests executed; seven errors**, all the already known missing historical ETS2/GSE84161 caches or historical Rscript runtime. No BACH2 failure occurred. The exact log is `work/psc-bach2-cohort-qualification/full-tests.log`; unrelated branches were not repaired. Catalog import logic is unchanged.
- All eight closed metadata/first-archive deliverable byte and SHA-256 pins are unchanged.

The independent first launch had a relative receipt-path bookkeeping error and stopped before archive parsing. It was fixed with a regression test. Its subsequent annotation wording was expanded from CD-only recognition to all readable source marker names. A comparison-key repair records that the main TAR API strips a directory trailing slash while the manual header reader preserves it. Safe canonical names, types, offsets and byte hashes all agree; both displays remain recorded. None of these repairs changed source values, RNA axes, retention or analytical scope.

```text
.venv/bin/python scripts/qualify_psc_bach2_cohort.py qualify
.venv/bin/python -m unittest discover -s tests -p test_psc_bach2_cohort.py -v
```

Qualification replays offline from pinned sources. Decoding caches are rehashed and checked against original TAR payloads. Matrices, full axes, individual technical margins, source joins and independent details remain ignored under `work/psc-bach2-cohort-qualification/`. Only aggregate qualification and feature-target annotation are public outputs. No new network search, raw-feature/read/liver-count body, normalization, score, model, gene contrast, genotype inference, clinical merge, contact, private access, installation or Git operation occurred.

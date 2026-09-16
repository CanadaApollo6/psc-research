# GSE239283 PSC organoid input qualification

## Status and scope

Source acquisition, numeric checks and donor–library joins are complete. No gene-wise test, target inspection, program score, treatment contrast or model fit was performed. An analysis freeze remains the root agent’s responsibility.

The design contains **8 donor keys** (4 PSC and 4 non-PSC), **8 complete donor-treatment pairs**, and **16 source libraries**. Sixteen libraries are not sixteen donor pairs. Each donor contributes CNT and IL17 at passage 3. GEO and supplemental Methods specify 100 ng/mL recombinant human IL-17A for 24 hours versus vehicle; vehicle composition is not stated.

## Critical assay distinction

The merged integer matrix has 21,562 features × 46,343 cells and 46,282,529 entries. Every cell’s total and positive-feature count exactly matches **nCount_SCT and nFeature_SCT**. Only 73 cell totals match nCount_RNA. Its total is 89,760,200; the metadata RNA total is 1,128,472,969. The merged matrix must **not** be used as unnormalized RNA input to a count model. Integer values do not establish raw RNA provenance.

Methods state SCTransform (v2) and Harmony (v0.1.1). The margin match supports an SCT-assay-consistent quantity, not recovery of the exact historical assay/slot export command. The complete merged input and its separately labeled library sums are preserved.

## Original processed RNA exports

Under the expanded cumulative 1-GB processed-input authorization, the complete original-archive-members.json ledger accounts for 48 approved GEO TAR members: one matrix, feature list and barcode list per source library. The archive contains processed 10x count exports, not FASTQ reads. All members are regular basename-only files. Extraction rejects extra/duplicate members, path traversal, absolute paths, symlinks, hardlinks and unexpected sizes.

All 16 matrices have the same **36,601 unique Ensembl gene IDs**. Their 36,591 distinct gene labels are retained separately; duplicate symbols are not collapsed. The original exports contain 48,164 cells and 195,537,244 stored entries. Every released cell maps exactly once within its stated source library after a separately recorded technical merge-suffix removal. The analysis-preparation subset contains all **46,343 released cells**. The other 1,821 source cells remain outside that published subset, without a newly invented exclusion rule.

The all-feature retained-cell count total is **1,128,543,232**, versus **1,128,472,969** in historical nCount_RNA metadata: a difference of **70,263 counts** and **68,460 positive entries**. Totals and positive-feature counts both match exactly in 20,540 cells and differ in 25,803. Original-minus-metadata ranges are 0 to 132 counts and 0 to 64 features per cell. No source value was changed. No guessed historical gene filter was applied to force agreement. This is an unresolved historical feature-universe/export difference, not evidence that the integer SCT export was raw RNA.

Before effects, the root agent accepted the **complete original RNA feature universe for the same published retained cells**. This is a new released-raw-universe estimand, not exact reconstruction of the historical filtered RNA matrix. The small metadata discrepancy does not prove its cause. All 36,601 unique source IDs are retained; 10 duplicated symbol labels cover 20 separate feature rows and are not collapsed. The JSON eligibility decision binds the exact count, feature, library and cell-join hashes and independent checks. A separate final root analysis freeze is still required; this qualification does not authorize effects.

## Global and library QC

The merged matrix and all original matrices pass gzip CRC/EOF, dimensions, finite nonnegative integer entries, coordinate bounds and duplicate checks. Neither contains explicit stored zeros. Every feature and released cell is retained. The original-RNA subset has 28,960 positive-total features and 7,641 zero-total features; the merged SCT matrix has 1,208 zero-total features. These are availability counts, not gene filters or biological nulls.

| Source library | Released cells | Original cells | Full-feature RNA total in released cells |
|---|---:|---:|---:|
| PSC_4_CNT | 5,737 | 5,822 | 88,703,431 |
| PSC_4_IL17 | 1,886 | 2,020 | 67,949,998 |
| PSC_6_CNT | 4,076 | 4,120 | 75,718,718 |
| PSC_6_IL17 | 2,115 | 2,227 | 62,858,224 |
| PSC_8_CNT | 3,409 | 3,425 | 80,228,717 |
| PSC_8_IL17 | 7,051 | 7,078 | 86,158,409 |
| PSC_9_CNT | 4,859 | 4,909 | 95,818,240 |
| PSC_9_IL17 | 1,527 | 1,585 | 56,103,957 |
| nonPSC_2_CNT | 1,055 | 1,364 | 60,249,667 |
| nonPSC_2_IL17 | 1,965 | 2,217 | 73,469,208 |
| nonPSC_3_CNT | 464 | 521 | 40,434,995 |
| nonPSC_3_IL17 | 1,524 | 1,723 | 70,337,228 |
| nonPSC_4_CNT | 1,745 | 1,843 | 59,528,685 |
| nonPSC_4_IL17 | 2,323 | 2,420 | 60,928,708 |
| nonPSC_5_CNT | 1,625 | 1,862 | 63,847,999 |
| nonPSC_5_IL17 | 4,982 | 5,028 | 86,207,048 |

| Condition | Treatment | Donors | Cells | Full-feature RNA total |
|---|---|---:|---:|---:|
| PSC | CNT | 4 | 18,081 | 340,469,106 |
| PSC | IL17 | 4 | 12,579 | 273,070,588 |
| nonPSC | CNT | 4 | 4,889 | 224,061,346 |
| nonPSC | IL17 | 4 | 10,794 | 290,942,192 |

Full released-matrix library totals are computed before any future feature filter. Cell barcodes, source joins and feature-by-library expression remain in ignored `work/psc-organoid-inputs/`.

## Identity, annotation and release limits

Non-PSC donors are clinically indicated ERC procedure controls, not healthy volunteers. The larger clinical Table 1 covers 10 PSC and 7 non-PSC patients. The inspected GEO, article Methods/Table 1 and supplemental Methods do not independently crosswalk its rows to the eight sequencing labels. Age, sex, IBD, cirrhosis, exact source site and brush/bile identity therefore remain unassigned. Numeric donor suffixes were not treated as table-row identifiers. One bile-derived organoid is described without its sequencing key. The text/table brush-versus-bile tension is preserved, not resolved by assumption.

GEO declares hg38 and Cell Ranger 5.0.0. Methods declare Seurat v4, SCTransform v2 and Harmony v0.1.1. The exact genome reference bundle, annotation release, chemistry subversion, feature prefilters and executed export commands are not established. Source labels remain unchanged; the original Ensembl IDs are not replaced by current aliases. The merged gene list contains two identical symbol-style columns rather than stable Ensembl IDs.

The supplemental cell-QC rule is **nCounts < 200 or mitochondrial reads > 40%**, not nFeatures < 200. It is not assumed to explain all absent cells or the historical RNA feature difference. No new cell, gene or donor exclusions were applied. Code is offered on author request in the inspected sources; no request or contact was made.

## Reproduction and files

```text
.venv/bin/python scripts/qualify_psc_organoid_inputs.py --acquire
.venv/bin/python scripts/qualify_psc_organoid_inputs.py --acquire-original
.venv/bin/python scripts/qualify_psc_organoid_inputs.py --qualify --qualify-original --finalize
.venv/bin/python -m unittest tests.test_psc_organoid_inputs -v
```

- Aggregate machine-readable result: `data/derived/psc-organoid-input-qualification.json`.
- Source URLs, original retrieval times, byte counts and SHA-256: `data/raw/psc-organoid-il17/source-manifest.json` and the public JSON’s `sources` field.
- Safe member inventory/hashes: `data/raw/psc-organoid-il17/original-archive-members.json`.
- Raw-source matrix: `work/psc-organoid-inputs/raw-rna-pseudobulk.tsv.gz`, with `source_feature_id` then 16 exact source-library columns.
- Axis and library/QC tables: `raw-rna-features.tsv` and `raw-rna-libraries.tsv` in the same ignored directory.
- Cell crosswalk and independent verification stay ignored. The public report contains only study/library aggregates.

Acquired/copied source bodies total 879,495,941 bytes under the 1,000,000,000-byte limit. No raw reads, controlled-access data, researcher contact, model prediction or clinical recommendation is involved.

The sample size remains four donors per group. More cells do not create more independent donors. This stage establishes input provenance and numeric integrity, not an IL-17 disease mechanism or treatment benefit.

## Validation

All 46 organoid synthetic tests pass. The independent SciPy sparse implementation reproduces all 585,616 feature-by-library values and both per-cell RNA margins, with 0 mismatches across 2,806,150 comparisons. A second full pass uses different chunk boundaries; 9 count/axis/join tables reproduce byte for byte.

The required full repository suite ran 402 tests and reported 7 errors from unavailable pre-existing ETS2/GSE84161 source caches or its separate Rscript runtime. These files were not downloaded, installed or altered. The organoid tests themselves pass; a claim that the entire repository suite passed would be incorrect.

# ETS2 macrophage context-dataset search and readiness audit

This is the bounded public-data audit for generic human macrophage stimulation/context calibration of the fixed ETS2 program transfer. It follows the effective prospective amendment in [`config/ets2-validation-amended-rules.json`](../config/ets2-validation-amended-rules.json), SHA-256 `cbbaa6e626dad9f494367cf1a8c45c0db270f4c330f9bcd61a72280fe1b8ef6b`; the original rules remain preserved in [`config/ets2-validation-rules.json`](../config/ets2-validation-rules.json), SHA-256 `a46c1378a186127bcc784cc5fea4dff48d963e4f693692bb7a4a43637429a6e7`.

The search stopped after a five-entry shortlist. The versioned lookup boundary was accession-controlled GEO metadata and file listings for GSE255234, GSE193336, GSE46903, GSE47189 and GSE61298, plus Zenodo release 11563707 (MacroMap); for GSE193336 and GSE255234 it also included the advertised SOFT and processed-matrix URLs, and for GSE255234 the linked primary article and Europe PMC XML. The exact URLs and stopping boundary are retained in [`work/ets2-context-search/source-notes.md`](../work/ets2-context-search/source-notes.md). Initial exploratory keyword queries and discarded-hit records were not fully archived, so this is a curated shortlist with no exhaustive negative result. No FASTQ/SRA reads, ETS2 values, frozen-program values, differential-expression results, target effects or program scores were inspected.

The exploratory hits not promoted into the five-entry inventory were GSE85333 (array/drug-arm route), GSE79077 (array route with three visible donor labels), GSE235897 (replicate labels without explicit donor identity) and GSE35449 (array M0/M1/M2 route). These names and the observed exclusion reasons are retained for transparency, but their initial keyword query strings, pagination/caps and full discarded-hit records were not archived. They must not be interpreted as exhaustive exclusions.

## Readiness decision

GSE255234 is the strongest current fallback for a **Tier-C generic context calibration**: it has an explicit four-donor design with control, LPS and two heme arms, and its small public workbook passes the file-level nonnegative integer-count check. It remains pending for the unchanged primary procedure because the methods describe quantification against hg38 RefSeq genes while the supplied matrix contains Ensembl identifiers. The full count-universe annotation and denominator provenance must be resolved in the dataset-specific freeze. The second heme arm is labeled only “Heme (lower concertration)” in GEO; its exact dose must remain source-label-only until frozen.

GSE193336 is a useful paired-donor reserve, with four explicit unstimulated/LPS pairs, but the supplied file is not an integer matrix as stored. It contains 66,295 numerically fractional cells (65,146 decimal values and 1,149 scientific-notation values), including values such as `3.366` and `5.2904e-09`. It requires a separately frozen noninteger-expression procedure and cannot enter the unchanged full-universe `log2(CPM+1)` workflow.

MacroMap is a potentially valuable large iPSC-derived context resource, but the raw expression matrix and exact donor/condition pairing were not downloaded or verified in this bounded pass. GSE61298 is a normalized-array reserve and needs a separate assay procedure plus donor-label confirmation. GSE46903/GSE47189 is a definite source-reuse audit because GSE46903 is the expression SubSeries of GSE47189, which the Stankey paper uses for TPP model selection and ETS2 co-expression; it cannot be used as independent validation.

The shortlist is calibration-only. Generic inflammatory or metabolic responses do not establish ETS2 dependence, disease activity, treatment efficacy or causal direction. Tier C cannot replace the open direct ETS2 validation question in Tier A.

## Candidate inventory

| Candidate | Biological design | Assay/file status | Independence and maximum claim | Current disposition |
|---|---|---|---|---|
| [GSE255234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255234) | Primary human monocyte-derived macrophages; four donors; 16 samples; each donor has control, LPS, H10 and H5 | Bulk RNA-seq workbook; 64,253 unique Ensembl IDs; 1,028,048 nonnegative integer cells; source methods cite RefSeq, so full annotation/count-universe provenance is unresolved | No reuse identified in the inspected Stankey narrative; individual cross-study donor linkage is unavailable; maximum claim is donor-paired generic LPS/heme context calibration | Recommended Tier-C fallback, pending dataset-specific freeze |
| [GSE193336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193336) | Primary human macrophages from CD14-positive monocytes; four donors; eight samples; unstimulated/LPS pairs | Bulk RNA-seq CSV; 58,825 unique Ensembl IDs; finite/nonnegative but 66,295 numerically fractional cells | No reuse identified in the inspected Stankey narrative; individual cross-study donor linkage is unavailable; maximum claim is separate noninteger-expression LPS calibration | Reserve only; raw-integer gate failed |
| [MacroMap Zenodo 11563707](https://zenodo.org/records/11563707) | iPSC-derived macrophages from 209 healthy unrelated HipSci lines; 4,698 libraries; 24 precursor/control/stimulation-time conditions | Public raw-expression route advertised at 174.1 MB (MD5 recorded), but matrix/header type and complete donor map not checked | No reuse found in the inspected Stankey narrative; iPSC transfer limits primary-blood generalization; maximum claim is iPSC stimulation-context calibration | Metadata shortlist only |
| [GSE61298](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61298) | Primary human monocyte-derived macrophages; 21 array samples; visible rep1–3 labels for GM-CSF/M-CSF lineages | PrimeView Affymetrix array; processed/normalized data and CEL files; no public raw integer count matrix | Separate study; individual cross-study donor linkage is unavailable; maximum claim is within-study normalized-array context calibration | Separate-array reserve; not ready for unchanged count procedure |
| [GSE46903 / GSE47189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46903) | Xue macrophage activation series; 384 GSE46903 samples; 299 macrophage transcriptomes across 28 conditions | HumanHT-12 array text, log2/quantile-normalized and non-normalized; no raw integer counts | GSE46903 is an expression SubSeries of GSE47189; Stankey uses GSE47189's 314 expression files/67 activation conditions for TPP selection and ETS2 co-expression | Tier R reuse/reproduction audit; excluded from independent selection |

## GSE255234: four-donor matrix and sample map

The GEO/SOFT series is titled “Distinct metabolic responses to heme in inflammatory human and mouse macrophages - Role of nitric oxide” and reports primary hMDM bulk RNA-seq. The linked article is [Pradhan et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11130737/) ([DOI](https://doi.org/10.1016/j.redox.2024.103191)). GEO exposes the processed workbook [`GSE255234_raw_counts_for_rna_seq_submission.xlsx`](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255234/suppl/GSE255234_raw_counts_for_rna_seq_submission.xlsx), 3,938,694 bytes at retrieval.

The source narrative describes four individual human donors. The bulk transcriptome comparison in Figure 1 uses heme 10 µM or LPS 1 µg/ml alone for 6 hours. GEO labels the lower heme arm only as “Heme (lower concertration)”; the inventory calls that arm H5 as a stable internal label, but no exact unit is assigned. All three treatment-versus-control arms must be enumerated and frozen before any downstream calculation. The culture methods report human monocytes differentiated with 25 ng/ml M-CSF in RPMI with 5% AB serum.

| donor | control | LPS | H10 heme | H5 / lower-heme label |
|---|---|---|---|---|
| D1 | key 1; GSM8067098; `D1CON [1_S1]` | key 2; GSM8067099; `D1LPS [2_S2]` | key 3; GSM8067100; `D1H10 [3_S3]` | key 5; GSM8067101; `D1H5 [5_S5]` |
| D2 | key 6; GSM8067102; `D2CON [6_S6]` | key 7; GSM8067103; `D2LPS [7_S7]` | key 8; GSM8067104; `D2H10 [8_S8]` | key 10; GSM8067105; `D2H5 [10_S10]` |
| D3 | key 11; GSM8067106; `D3CON [11_S11]` | key 12; GSM8067107; `D3LPS [12_S12]` | key 13; GSM8067108; `D3H10 [13_S13]` | key 15; GSM8067109; `D3H5 [15_S15]` |
| D4 | key 16; GSM8067110; `D4CON [16_S16]` | key 17; GSM8067111; `D4LPS [17_S17]` | key 18; GSM8067112; `D4H10 [18_S18]` | key 20; GSM8067113; `D4H5 [20_S20]` |

The workbook has one sheet, `raw_counts_for_rna_seq`, with dimension `A1:Q64254`. Cell A1 is blank. B1:Q1 contain exactly the internal keys `1,2,3,5,6,7,8,10,11,12,13,15,16,17,18,20`. Column A contains 64,253 unique Ensembl stable identifiers without version suffixes; no duplicate source IDs were observed. All 1,028,048 sample cells are present, finite, nonnegative and numerically integral, and every full-column library total is positive. The deposited matrix retains 31,285 all-zero feature rows. This is a file-level structural pass only and does not prove that the deposited feature universe exactly matches the article's original quantification export.

The methods and SOFT say HT-seq 0.11.1 with hg38 RefSeq reference genes, while the supplied workbook uses Ensembl identifiers. That mismatch does not prove the matrix is invalid, but it leaves the source annotation, exported full gene universe and denominator provenance unresolved. The candidate is therefore biologically suitable as Tier C and structurally promising, while assay readiness remains pending. The matrix may be used only after the prospective freeze records the identifier namespace, annotation version, all source-counted rows, any summary rows and the duplicate policy.

## GSE193336: explicit pairs but noninteger supplied matrix

The GEO series is titled “Bulk RNA-seq on unstimulated and LPS-stimulated human macrophages.” Its [GEO SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193336/soft/GSE193336_family.soft.gz) reports four anonymous Sanofi blood donors, CD14-positive monocyte differentiation with 50 ng/ml M-CSF for five days, then 50 ng/ml *E. coli* O111:B4 LPS for 24 hours. The matrix is [`GSE193336_GEO_Processed_Data_Raw_Count_Macrophage.csv.gz`](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193336/suppl/GSE193336_GEO_Processed_Data_Raw_Count_Macrophage.csv.gz), 757,738 bytes at retrieval.

| donor | unstimulated header/GSM | LPS header/GSM |
|---|---|---|
| donor 1 | `donor 1_unstim`; GSM5792562 | `donor 1_LPS`; GSM5792563 |
| donor 2 | `donor 2_unstim`; GSM5792564 | `donor 2_LPS`; GSM5792565 |
| donor 3 | `donor 3_unstim`; GSM5792566 | `donor 3_LPS`; GSM5792567 |
| donor 4 | `donor 4_unstim`; GSM5792568 | `donor 4_LPS`; GSM5792569 |

The gzip CSV has nine columns (`geneid` plus the eight headers above), 58,825 unique Ensembl stable identifiers and 470,600 sample cells. All values are finite and nonnegative, but 66,295 are numerically fractional: 65,146 ordinary decimal tokens and 1,149 scientific-notation tokens. Examples are `3.366`, `13.2141` and `5.2904e-09`. This verifies a measurement-type failure rather than merely a lexical-format concern. Fractional values alone do not prove that the matrix is normalized or identify how it was produced. The file is retained for a separately frozen noninteger-expression calibration and is excluded from the unchanged raw-integer procedure.

## MacroMap: large iPSC-derived context route

The [Zenodo 11563707 release](https://zenodo.org/records/11563707), associated with the [analysis repository](https://github.com/andersonlab/macromap_eqtl) and [publication context](https://pmc.ncbi.nlm.nih.gov/articles/PMC12391537/), advertises real expression files rather than only eQTL summaries. The release includes `MacroMap_raw_expression.txt.gz` with rounded display size 174.1 MB and MD5 `3e585dfeb097368db7eaab9889f76f58`; all release files total about 1.7 GB.

The publication metadata describe 209 healthy unrelated HipSci/iPSC donor lines, 4,698 RNA-seq libraries and 24 conditions: D0/D2 precursor states, Ctrl_6/Ctrl_24 and IFNB, IFNG, IL4, R848, PIC, P3C, sLPS, CIL, LIL10 and MBP at 6 and 24 hours. Condition-specific donor counts are reported as 177–202. The route is scientifically useful for stimulation-context transfer, but the exact matrix headers, integer status, full feature universe and complete control pairing were not checked here. It also uses iPSC-derived macrophages and two library protocols, so it cannot by itself establish transfer to primary blood macrophages. No expression file was downloaded under this bounded pass.

## GSE61298: normalized-array reserve

GSE61298 provides 21 visible samples from primary human monocyte-derived macrophages with GM-CSF and M-CSF lineage conditions. The visible metadata expose rep1–3 labels and 48-hour activation arms, including lineage controls, LPS+IFN-gamma, IL-4 and M-CSF IL-10. The public route is a PrimeView Affymetrix array with processed/normalized data and CEL files. The current metadata do not make individual donor identity explicit beyond replicate labels. It is retained as a separate normalized-array calibration option and is not assay-ready for the unchanged raw-count procedure.

## GSE46903/GSE47189: required overlap audit

GSE46903 is the expression SubSeries of GEO superseries GSE47189, which includes the Xue macrophage activation resource and other assay subseries. GEO reports 384 GSE46903 samples and 299 macrophage transcriptomes across 28 conditions. The public expression files are array text (re-analyzed about 19.8 MB and non-normalized about 103.4 MB); no raw integer count matrix is available in the inspected listing.

The Stankey Nature Methods text explicitly states that human monocyte/macrophage expression files from GSE47189 (314 expression files across 67 activation conditions) were quantile-normalized and used for TPP model selection, and that ETS2 co-expression was evaluated across those same 67 conditions. This makes the GSE46903/GSE47189 route definite source reuse, even though it is scientifically relevant macrophage context. It is retained as Tier R for reproduction/reuse calibration only and excluded from independent selection.

## Cache and provenance manifest

The two small matrices and their SOFT metadata are cached in [`data/raw/ets2-validation/context/`](../data/raw/ets2-validation/context/). The detailed cache README records direct routes and byte checks. The current machine-readable inventory is [`data/derived/ets2-validation-context-candidates.csv`](../data/derived/ets2-validation-context-candidates.csv), and the structured manifest is [`config/ets2-validation-context-sources.json`](../config/ets2-validation-context-sources.json).

| asset | retrieved (UTC) | remote bytes | local bytes | local SHA-256 |
|---|---:|---:|---:|---|
| GSE193336 matrix | 2026-09-12T16:08:56Z | 757738 | 757738 | `2fe82f9b39abd3e4826c74b9280d20fbd3a175d0110b875a35cb84210c7bb98d` |
| GSE255234 matrix | 2026-09-12T16:09:09Z | 3938694 | 3938694 | `61d23d009b9b9a70a6127e0ae024df18f83fa3ed1a8079b7cbe86854ac5aac18` |
| GSE193336 family SOFT | 2026-09-12T16:13:33Z | 2770 | 2770 | `f5ff63adebf17e2319741479c0fbb22e1a332407515371936e4ee4220f127441` |
| GSE255234 family SOFT | 2026-09-12T16:13:33Z | 3518 | 3518 | `9cdf727be9c0d82edecf38fbb36aa970c48e9dd55b0b6a4eae75005636002e8b` |
| Pradhan 2024 full-text XML | 2026-09-12T16:14:57Z | not pinned | 138676 | `110db4d098ddcb055201aa987b4b7e0c50fbbc5e2b18c63bd1d6fe98c379cd5` |
| attempted Pradhan PDF response | 2026-09-12T16:14:57Z | 1817-byte HTML preparation stub | 1817 | `e7ebea61478bb8ac6d9afc370074a39b75bc5995628b7f81aa81d39a0be6a982` |

The attempted PDF response is retained as a failed retrieval record and was not used as evidence; the XML is the usable article source cache. No publisher SHA-256 was advertised for the GEO assets. No raw sequencing reads were downloaded.

## Next gate before any calculation

If direct Tier-A/B evidence remains unavailable, GSE255234 is the candidate to freeze first as Tier C. The freeze must record the four-donor map above, the exact control-only reference stratum, all three treatment arms, the H5 source label and dose, the Ensembl/RefSeq annotation reconciliation, the full exported feature universe including measured zeros and any non-gene rows, and the duplicate policy. Related contrasts may reuse the control reference within the frozen stratum, but a treatment sample must never enter that reference. Only after this lock can the fixed-program mapping or any descriptive calibration values be read. No result from GSE255234 could establish ETS2 causality.

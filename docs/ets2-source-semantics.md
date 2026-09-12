# ETS2 source semantics audit

Snapshot: 2026-09-12. The deterministic audit in
[`scripts/audit_ets2_sources.py`](../scripts/audit_ets2_sources.py) verifies the
pinned public release and publisher supplements, then extracts source-preserving
tables. It does not calculate an atlas score, select genes by atlas performance,
or reconstruct donor-level count matrices.

## Pinned sources

The base manifest is
[`config/ets2-benchmark-sources.json`](../config/ets2-benchmark-sources.json).
Its 34,636,310-byte public ZIP has SHA-256
`f1da61f2d58b31b1ce1cc418486222a8701698adda76324a10f5f48ecc53239c` and the
publisher MD5 `9d6b7ddcf377cb194bc72732d807ae80`; both checks pass. The archive
root is `JamesLeeLab-ETS2manuscript_Stankey_CT_et_al_2024-34a9f49/`.

The supplementary manifest is
[`config/ets2-program-sources.json`](../config/ets2-program-sources.json).
The cached publisher workbook is 5,270,597 bytes with SHA-256
`5390818bf0eb81a027ac145ac65661ed1b6333299d21e257f7d7c3e4b7ad8c74`; source
data Fig. 5 is also verified (7,907,111 bytes,
`b55d3998747dbfd45384f8c11116965029b29c24d8b764f8e4b40811d8594b9a`). The
Fig. 5 workbook is retained in the manifest, but its nonperturbation donor
sheets are outside this audit.

The primary article XML is cached from
[Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11168933/fullTextXML)
(285,979 bytes,
`3abb35514eb3584180b5bf07c25445d74331cb66e32ac319a7050b32cbbb68c4`). In the
cached XML, paragraph 76 states that differential-expression results are in
Supplementary Tables 1 and 2, paragraph 77 states that GSEA lists are ranked by
t-statistic, paragraph 24 reports the CRISPR RNA-seq experiment as n = 8,
paragraph 82 describes the gRNA1 downregulated list, paragraph 87 names the
gRNA1 modular-expression program, and paragraph 123 gives EGA accessions.

## Gene-set definitions

`RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv` contains six ordered
columns of Ensembl gene IDs without version suffixes. The source counts are:

| source column | members |
|---|---:|
| ETS2g1_UP | 668 |
| ETS2g1_DN | 928 |
| CHR21_DN | 141 |
| CHR21_UP | 215 |
| ETS2g2_UP | 716 |
| ETS2g2_DN | 875 |

There are 2,184 unique IDs in the union, 1,168 IDs in more than one column,
and 1,359 repeated memberships across columns. No column has within-column
duplicates.

Every master column reproduces exactly from S1 using the strict rule
`adj.P.Val < 0.05` plus the sign of `logFC` for the matching contrast. This is
an FDR 0.05 definition; the audit reports FDR 0.10 directional counts
separately and does not substitute them. All 12,144 S1 rows have matching
`logFC` and `t` signs with no zero-sign exception. The master-to-S1 crosswalk
preserves the original Ensembl ID, symbol, row number, `logFC`, `t`, and
adjusted P value. The source-defined `ETS2g1_DN` list therefore has 928 members
before any downstream decision to exclude ETS2 itself.

## Publisher differential tables and identifiers

The S1 worksheet body has 12,149 row slots after its two header rows (worksheet
dimension `A1:T12151`), of which 12,144 are valid Ensembl rows; S2 has 10,201
body slots (dimension `A1:T10203`), of which 10,196 are valid Ensembl rows. The
audit identifies valid rows by an `ENSG` value in column A and keeps footer text separately in
[`ets2-source-notes.csv`](../data/derived/ets2-source-notes.csv). S1 contains
481 literal `None` symbols and duplicate non-`None` symbols `POLR2J3` (2) and
`Y_RNA` (3). S2 contains 199 literal `None` symbols and duplicate non-`None`
symbols `LINC03025` (2), `U2` (10), `POLR2J3` (2), `NPIPA9` (2), `NPEPPSP1`
(2), and `SIGLEC5` (2). `None` is preserved as a literal source value and is
not treated as a valid symbol or an ID mapping.

S1 contrasts are the published groups `ETS2 gRNA 1`, `ETS2 gRNA 2`, and
`ch21q22 deletion`, with the footer stating that values are relative to
unedited cells. S2 contrasts are 250 ng, 500 ng, and the combined dose-adjusted
analysis, with the footer stating that values are relative to an equivalent
amount of reverse-complement ETS2 mRNA. The full original rows and numeric
strings are in
[`ets2-source-differential.csv.gz`](../data/derived/ets2-source-differential.csv.gz);
no count matrix is inferred from these summaries.

## Rank statistics and sign provenance

Rank files are two-column, symbol-keyed vectors with finite, unique symbols and
strict nonincreasing statistics. The source code reads them as fGSEA statistics;
the article identifies the statistic as the limma t-statistic. The primary
numeric comparison uses only a one-to-one non-`None` symbol match. Symbols with
multiple source rows are unavailable for primary identity matching; nearest
numeric rows and all-candidate sign counts are retained as diagnostics only.

* The KO gRNA1 rank has 11,658 one-to-one S1 matches. All 11,658 are within
  Decimal tolerance `1e-9` of the S1 gRNA1 `t` value, with zero sign mismatches.
  Some raw xlsx values contain longer floating-point serialization than the
  rank CSV, so this is an within-tolerance match rather than an exact-string
  claim.
* Each OE rank has 9,908 one-to-one S2 matches. None of those 9,908 is within
  `1e-9` of the corresponding S2 `t` value. One-to-one sign mismatches are 36
  for 250 ng and 42 for 500 ng. Five shared symbols are ambiguous and 66 rank
  symbols have no S2 symbol row. For context only, counting every duplicate
  source candidate gives 41 and 50 sign mismatches; those counts do not choose
  an identity.
* The MEK inhibitor ranks have no pinned publisher differential table in this
  source set. The R code defines `MEK_100nM = MEK_100 - ctrl` and
  `MEK_500nM = MEK_500 - ctrl`, then reads the released rank files, but the
  matching `write.csv` lines are commented or use a different filename. The
  saved rank direction is therefore unresolved. The audit does not assign a
  sign from expected inflammatory or drug-response biology.

The full rank rows, original statistic strings, source rows and comparison
status are in
[`ets2-source-ranks.csv.gz`](../data/derived/ets2-source-ranks.csv.gz).

## Donor and design boundaries

The R code gives the following design metadata; these counts describe the
experiments and do not turn a rank vector into a donor matrix.

| experiment | code samples | model samples | donor design | code contrast |
|---|---:|---:|---|---|
| ETS2 CRISPR | 38 names | 27 selected | 10 donor labels; 8 complete gRNA1/NTC and gRNA2/NTC pairs (donors 2, 3, 5, 6, 7, 8, 9, 10) | `g1-NTC`, `g2-NTC` |
| ETS2 overexpression | 32 | 32 | 8 complete paired donors, two doses | `ets2_250-rev_250`, `ets2_500-rev_500`, and dose-adjusted `ets2-rev` |
| MEK inhibitor | 9 | 9 | 3 complete paired donors | `MEK_100-ctrl`, `MEK_500-ctrl` |

The CRISPR code labels the second guide factor `g2` while the literal sample
names are `g4`; that source mismatch is retained. The CRISPR code also carries
partial donor labels 1 and 4 in the selected model, so “10 donor labels” does
not mean 10 complete pairs. The overexpression rank files are split-dose
outputs from one eight-donor experiment, not independent cohorts.

The CRISPR, overexpression and MEK inhibitor folders do not contain the
`raw_counts.txt` files referenced by their scripts. The public biopsy folder's
count file is a different experiment and is not substituted. The article lists
controlled-access EGA routes for CRISPR (`EGAD00001011338`), overexpression
(`EGAD00001011341`) and MEK inhibitor (`EGAD00001011337`). Thus the audit
supports summary-level provenance and source reproduction, not donor-level
limma re-fitting or held-out-donor validation.

## Output tables

[`reports/ets2-source-audit.json`](../reports/ets2-source-audit.json) records
all verified hashes, source member hashes, set counts and overlaps, workbook
row/duplicate/`None` metrics, rank comparisons, design metadata, caveats and
the exact output columns. The compact extracted outputs are:

* `ets2-source-gene-sets.csv`: ordered master memberships and S1 crosswalk.
* `ets2-source-differential.csv.gz`: all valid S1/S2 rows for all three
  published contrasts per worksheet, with original values and signs.
* `ets2-source-ranks.csv.gz`: all five rank vectors with source-match status
  and code/hash provenance.
* `ets2-source-design.csv`: every literal code sample name and model inclusion
  assignment.
* `ets2-source-notes.csv`: exact worksheet footer text.

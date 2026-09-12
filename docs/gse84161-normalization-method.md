# GSE84161 normalization and gene mapping

This is the implementation specification for the normalization and mapping requested September 12, 2026. It supplements the preserved [array framework](gse84161-array-protocol.md). It permits expression preparation and assay QC, without calculating program scores, treatment contrasts or paired donor effects. The [execution plan](../config/gse84161-normalization-plan.json) must pin the complete sources, scripts and runtime before execution.

## A necessary namespace clarification

Inspection of the original source columns before new mapping showed that the four ETS2/CHR21 sets use Ensembl gene identifiers, while the four GOBP comparison sets use gene symbols. The earlier framework's blanket wording about Ensembl members was incomplete. The source sets themselves and prior analyses are unchanged.

| Program | Source namespace | Unique source members | Intended after ETS2 exclusion |
|---|---|---:|---:|
| ETS2 gRNA1 down | Ensembl gene ID | 928 | 927 |
| ETS2 gRNA1 up | Ensembl gene ID | 668 | 668 |
| ETS2 gRNA2 down | Ensembl gene ID | 875 | 875 |
| Chr21 enhancer deletion down | Ensembl gene ID | 141 | 141 |
| Inflammatory response | Gene symbol | 815 | 815 |
| Response to interferon gamma | Gene symbol | 139 | 139 |
| Response to oxidative stress | Gene symbol | 436 | 436 |
| Apoptotic signaling | Gene symbol | 588 | 588 |

Use exactly these source columns from the pinned Stankey release. Remove ETS2 by its exact Ensembl identifier or literal symbol, then deduplicate exact source identifiers. These source columns already contain unique identifiers. Only numeric version suffixes on Ensembl identifiers may be removed.

## Fixed mapping

Use the already pinned GPL570 snapshot: all 54,675 probe sets, with October 6, 2014 row annotation dates. Split Entrez fields on `///`, strip whitespace, discard empty/unresolved tokens, and count unique IDs within each probe. Retain a probe only if exactly one valid Entrez identifier remains. Preserve the original field and exclusion reason for every probe. The preceding metadata audit found 41,834 eligible probes representing 20,486 Entrez genes; discrepancies from that accounting must be resolved before aggregation.

Use a complete human Ensembl release 116 archive export containing Ensembl gene ID, HGNC symbol and Entrez Gene ID, preserving the query, response completeness evidence and complete byte hash. This is a mapping of identifiers within humans, not cross-species orthology. [Official Ensembl ID-conversion guidance](https://jun2026.archive.ensembl.org/Help/Faq?id=125).

Build the full Ensembl-to-Entrez and Entrez-to-Ensembl relations before restricting to the array or source programs. An eligible relation has exactly one Entrez ID for that Ensembl gene and exactly one Ensembl gene for that Entrez ID. Repeated identical source rows are duplicate observations of the same relation, not multiple mappings. Blank values never become identifiers.

For the four symbol-defined programs, first require exactly one Ensembl gene for the exact, case-sensitive HGNC symbol across the full human export; then apply the same reciprocal one-to-one Ensembl–Entrez rule. No aliases, previous names, symbol-case changes, GPL gene-name shortcuts or expression-based choice are allowed. Preserve missing symbols, ambiguous symbols, missing Ensembl IDs, ambiguous crosswalk relations and absence of an eligible array probe as separate outcomes. Absence from this snapshot alone does not prove that an identifier was retired.

An eligible mapped program member must end at a gene represented by at least one eligible array probe. Report all losses against the original intended denominator. Retain the 80% coverage and ten-gene minimum gates. Derive the ninth program by subtracting the mapped comparator union from the mapped primary program in Entrez space; retain the primary source coverage gate and a separate ten-gene minimum for the remaining program. A failed program remains unavailable; no threshold relaxation or replacement is planned.

## Array normalization and gene aggregation

Authenticate the original archive and every extracted CEL against the saved inventory. Use all 30 monocyte arrays in the saved sample order, with explicit GSM-to-file correspondence. The candidate replicate labels remain metadata; this step neither confirms donor identities nor performs donor adjustment.

Run official Bioconductor `affy::rma` on all arrays together, with the full HG-U133 Plus 2 CDF, background correction and quantile normalization enabled, and `bgversion=2`. The expected result is 54,675 unique probe-set rows by 30 unique GSM columns with finite log2 RMA values. Preserve the entire normalized probe matrix, including probes excluded from gene aggregation. [Official affy reference](https://www.bioconductor.org/packages/release/bioc/manuals/affy/man/affy.pdf).

For each of the 20,486 eligible Entrez genes, take the per-array median of all its eligible normalized probe-set values. Do not select a high-variance probe or weight genes by probe count. Keep genes without an eligible Ensembl crosswalk in the full Entrez matrix, with mapping status recorded separately. Do not apply CPM, a second logarithm, donor residualization or additional normalization after median aggregation. Record output precision, matrix dimensions, feature/sample ordering, contributing probe counts and complete file hashes.

Full normalized matrices and extracted CELs stay in ignored workspace storage. Small mapping, coverage, QC and provenance summaries are versioned. Runtime installation uses a separate local environment; the earlier coloc environment and all previous report/protocol snapshots remain unchanged.

## QC defined before expression inspection

Hard stops are corrupted or unexpected source files, mismatched sample/probe identities, a non-finite or negative raw PM intensity, and non-finite or incorrectly shaped normalized matrices. Zero raw PM values are counted and disclosed; they are not automatically discarded. Do not turn missing values into zero.

1. **Raw input:** count CDF-defined perfect-match (PM) probes per array and summarize their intensity distribution using quantiles 0, 1, 25, 50, 75, 99 and 100 percent. Record zero, negative and non-finite counts. PM summaries use the perfect-match cells passed to RMA, rather than all CEL grid cells or mismatch cells. Use R's quantile type 7. A display may use `log2(PM + 1)`; RMA receives the original PM quantities.
2. **Normalized distribution:** record the same quantiles across all normalized probe sets and show per-array distributions.
3. **Relative log expression (RLE):** subtract each probe set's median across all 30 arrays from its normalized log2 RMA values. Report each array's RLE median and interquartile range, and show per-array boxplots. This is a diagnostic calculation only; its residuals do not replace the expression matrix.
4. **Exploratory sample structure:** if included, PCA uses all normalized probe sets, centered per probe and without variance scaling or selection of high-variance probes. Any plot uses GSM labels and remains descriptive; it cannot establish donor correspondence, batch effects or treatment efficacy.

There is no automatic signal-based exclusion or universal QC pass threshold in this run. Valid arrays are retained; distribution shifts can have technical or biological causes. Report unusual patterns as unresolved diagnostic observations. Any proposed new exclusion needs an explicit method amendment before dependent program effects are inspected. Passing these checks establishes usable numeric inputs, not absence of technical confounding or biological validation.

## Completion criteria

Deliver a complete 54,675-by-30 normalized probe matrix, a 20,486-by-30 Entrez matrix, full probe/gene mapping accounting, all nine program coverage outcomes, per-array QC, software/source hashes and an independently checked reproduction record. Report missing or ambiguous mappings rather than manufacturing agreement. Program scoring and donor-paired treatment interpretation remain later steps under their preserved readiness requirements.

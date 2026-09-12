# MEK-inhibitor rank reconstruction

Snapshot: 2026-09-12. This is a bounded forensic reconstruction of the two
released MEK-inhibitor rank vectors. It asks which production steps are directly
supported by the cached R code, which properties can be checked from the
released files, and which steps cannot be recovered without controlled-access
counts. It does not change the benchmark's unresolved MEK direction and does
not calculate a treatment score or clinical response.

The reproducible entry point is
[`scripts/reconstruct_ets2_inhibitor.py`](../scripts/reconstruct_ets2_inhibitor.py).
The report is
[`reports/ets2-inhibitor-reconstruction.json`](../reports/ets2-inhibitor-reconstruction.json);
the step-level evidence, rank summaries, overlap checks, pathway checks and
source notes are in the five `ets2-inhibitor-reconstruction-*.csv` files under
[`data/derived`](../data/derived).

## Evidence chain

The pinned release ZIP passes its SHA-256 and publisher MD5 checks. The relevant
MEK code, ranks, nine-column gene-set file and archived 500 nM fGSEA output are
hashed in the report. The primary article XML states that differential lists
were ranked by t-statistic and provides the MEK EGA accession. Source Data Fig.
5 workbook MOESM9 is pinned separately; sheet `5j` reports nine pathway NES
values and states that they are for MEK inhibitor-treated versus vehicle-treated
cells. Its notes are preserved verbatim in
[`ets2-inhibitor-reconstruction-notes.csv`](../data/derived/ets2-inhibitor-reconstruction-notes.csv).

The code directly supports the following intended pipeline:

1. Read `raw_counts.txt` with row names and headers.
2. Use sample labels `s1`, `s2`, `s3`, each with `100`, `500` and `ctrl`.
3. Retain genes with `cpm(eset) > 0.5` in all nine samples.
4. Apply edgeR `calcNormFactors`.
5. Fit a no-intercept design with drug level and donor terms.
6. Form `MEK_100nM = MEK_100 - ctrl` and `MEK_500nM = MEK_500 - ctrl`.
7. Run `voom`, `lmFit`, `contrasts.fit` and `eBayes`.
8. Request `topTable(... number=12801, p.value=1, lfc=0,
   adjust.method="BH")` for each coefficient.
9. Read the released headerless symbol/statistic rank files into named vectors
   for fGSEA with `eps=0.0`, `minSize=15` and `maxSize=500`.

The exact source line, line number, confidence and limitation for each step are
in `ets2-inhibitor-reconstruction-steps.csv`. The MEK folder contains no
`raw_counts.txt`, so steps 1–8 are intended code semantics rather than a
re-runnable fit. The code creates topTable objects but the matching `write.csv`
lines are commented or use a different filename. It then reads the rank files.
The MEK script has no symbol-mapping call; analogous CRISPR code uses biomaRt,
but the mapping release, one-to-many handling, duplicate collapse, tie handling
and rank serialization for MEK cannot be identified.

## Observed rank vectors

Both files are finite, unique-symbol, strictly descending vectors with no
adjacent ties and 12,095 rows. The 100 nM vector has 6,080 positive and 6,015
negative values, ranging from −15.84996246 to 7.841975336. The 500 nM vector
has 6,150 positive and 5,945 negative values, ranging from −19.60553592 to
9.062092728. Full original rank strings are already retained in
[`ets2-source-ranks.csv.gz`](../data/derived/ets2-source-ranks.csv.gz); the
reconstruction rank table provides hashes, tails, quantiles and scale checks.

Symbol overlap is descriptive only. The two MEK vectors share all 12,095
symbols (Pearson statistic correlation 0.916). MEK 100 nM shares 11,041
symbols with KO gRNA1 and 9,350 with each OE vector; MEK 500 nM has the same
counts. Sign agreement across shared symbols is 0.673 and 0.688 against KO,
and approximately 0.503, 0.492 and 0.493 against the OE comparisons in the
reported pairings. These comparisons cannot validate treatment direction or
biological similarity because rank signs are unresolved and experiments differ.

The values are plausible as signed moderated statistics because the article
describes t-statistic ranking, the code consumes the files as fGSEA statistics,
and the vectors have a two-sided finite range with both signs. A t statistic is
a coefficient divided by a moderated standard error. Without the underlying
counts, fitted standard errors and logFC columns, its magnitude cannot recover
logFC, expression, variance, dose response or an IC50. The ranks cannot support
a treatment-efficacy claim.

## Signed ES compatibility check

The released MEK gene-set file contains nine original GO columns. I calculated
the signed weighted running-sum ES from the 500 nM rank vector using absolute
statistic weights with exponent 1 and equal miss steps. I did not calculate a
permutation null or NES. For every pathway:

* the computed ES matches the archived `fgsea_PD_0.5_mac.pathways.csv` ES within
  `1e-9` (maximum absolute error `1.63e-13`);
* the computed ES sign agrees with the source-data Fig. 5 sheet `5j` NES sign;
* the sheet `5j` NES matches the archived NES within `1e-9` (maximum absolute
  difference `2e-16`); and
* reversing the rank order while negating the statistics reverses all nine ES
  signs, as a deterministic orientation-sensitivity check.

This chain makes the released vector, archived ES and publisher NES numerically
compatible. Combined with the code-defined contrast and the sheet 5j
treated-versus-vehicle note, it supports a drug-minus-vehicle working
interpretation for the 500 nM pathway result. It does not recover the missing
export step or independently certify that the saved per-gene statistic sign
means inhibitor minus vehicle. The exact historical per-gene direction and the
100 nM export linkage remain unavailable from the public release.

## Non-identifiable quantities

The public release does not permit recovery of the filtered count matrix,
normalization factors, voom weights, fitted coefficients, moderated standard
errors, adjusted P values, donor-level residuals, symbol mapping version or
the 706-row difference between the 12,801-row topTable request and each
12,095-row saved rank. The article's EGA route is controlled access. The
biopsy count file in another archive folder cannot substitute for this MEK
experiment. No R package installation or controlled-data request was made.

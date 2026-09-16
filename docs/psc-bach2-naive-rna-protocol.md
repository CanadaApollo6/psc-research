# BACH2: one fixed naive-CD4 released-RNA observation check

## Question and prospective boundary

**This protocol does not authorize expression access.** The primary implementation, independent method review, native-output auditor and complete source/runtime seals must be accepted and frozen by root before either real extraction or real-output verification. A preparation decision or successful synthetic test is not that authorization.

Question: what difference and uncertainty in donor-level released BACH2 relative RNA are supported between the source SNP and noSNP groups in the paper's broad naive-CD4 C0/C1 compartment? The [prospective claim review](../reports/psc-bach2-followup-design-review.md) selected this single, modest-value observation check before new endpoint values. It is not a gene/state search or an independent replication.

Figure S2B explicitly names the C0/C1 compartment. This protocol does not reconstruct the four later activation subclusters in Figure 3 or the healthy-donor miR4464/mock transfection experiment in Figure 5C. Historical normalization, pooling and test units are insufficiently specified for exact Figure S2B method replication. Published results were available during design; this is not a blinded, untouched validation dataset.

## Fixed source, population and feature

Use only the eight qualified filtered archives from [E-MTAB-14013](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14013), from [Poch/Bahn et al.](https://doi.org/10.1016/j.xcrm.2024.101620). The [metadata qualification](../reports/psc-bach2-input-qualification.md), [first-archive qualification](../reports/psc-bach2-count-qualification.md) and [complete cohort qualification](../reports/psc-bach2-cohort-qualification.md) preserve the original values, acquisition dates, URLs, retrieval hashes and unresolved source history. Later qualification supersedes the earlier phase's not-yet-downloaded status, not its historical record.

- Keep all eight accession-local donor units: four source `SNP` and four `noSNP`. All are male PSC cases. Groups are unpaired. Cells and genes are not independent replicates.
- Source donor independence is an assumption supported by the deposited distinct-individual labels and study design, not an independent kinship or identity investigation.
- Start with the same **55,460 retained annotation rows**. Do not add any of the 3,855 source-only cells or reconstruct undocumented historical filtering.
- Select only `paper_clusters == "C0: CD4+ TN RTE"` or `"C1: CD4+ TN mature"`. Require the deposited short-label/Seurat crosswalks C0/0 and C1/1. This fixes **9,320 + 7,339 = 16,659 cells**, with both states represented in every donor.
- Join on sample plus barcode after stripping only that exact sample's underscore prefix. A bare-barcode join is invalid: 1,769 barcode strings recur across samples.
- Select the exact, unique source triple **`ENSG00000112182 / BACH2 / Gene Expression`**. It is at **source feature row 12,314, one-based**, in all eight identical ordered feature axes. This row number is not a genomic coordinate.
- Use all **36,601 released RNA-feature rows** for denominators in these same selected cells. Exclude all 38 ADT rows. Preserve implicit zeros and duplicate symbols elsewhere; do not collapse or remap the chosen feature.
- The reciprocal BACH2/Ensembl identity was checked against the cached full Ensembl 116 reference. This does not recover the executed source annotation, alter the released feature, or identify a variant coordinate or allele.

The original material was peripheral CD3-positive T cells. That does not make all CD3 cells eligible for this endpoint. E-MTAB-14103 liver metadata have no qualified genotype or participant bridge and supply no replication samples here. No other gene, state, cohort, protein proxy, genotype-state composition test or clinical endpoint is added.

## Released quantity and estimand

For each donor d, pool all fixed C0/C1 cells before transformation:

```text
B_d = sum of released BACH2 RNA-feature units in those cells
R_d = sum of released units over ALL 36,601 RNA rows in those SAME cells
Y_d = log2(1 + 1,000,000 * B_d / R_d)
Delta = mean(Y_d over the four SNP donors) - mean(Y_d over the four noSNP donors)
```

Give donors equal weight. Do not weight by cell count, depth, precision or clinical measurements. Do not equalize C0/C1 fractions or downsample cells. Pooling intentionally retains the observed state mixture and RNA-content weighting. It cannot distinguish population redistribution from within-state regulation.

This is a difference of donor-average transformed relative RNA quantities. It is not a raw BACH2 RNA fold change, absolute RNA amount, transcription rate, molecule count or protein endpoint. The pseudocount is exactly one after CPM scaling; no alternative pseudocount or normalization is permitted.

Released source-count-format units are not certified untouched original UMIs. The Cellranger 6.1.1 MEX header versus Methods/IDF Cellranger 3.0.2 and reference-release history remain unresolved. The original-UMI/count-likelihood gate remains unpassed. This protocol therefore does not use NB or another count likelihood.

## One nominal Welch inference

For each group g, calculate the donor sample variance with divisor three. With v_g = s_g^2/4:

```text
SE = sqrt(v_SNP + v_noSNP)
df = (v_SNP + v_noSNP)^2 / (v_SNP^2/3 + v_noSNP^2/3)
CI95 = Delta +/- t_quantile(0.975, actual df) * SE
P = 2 * t_survival(abs(Delta/SE), actual df)
```

This is **one gene × one fixed population × one contrast**. Report the nominal two-sided 95% model-based interval and P, with actual Satterthwaite df. There is no genome-wide claim, secondary discovery family or project-wide error-control claim. Do not substitute infinite-df normal inference.

The working model assumes independent donors and approximately Gaussian Y within source groups. Four donors per group cannot establish calibration, representative sampling, random genotype allocation or clinical exchangeability. No normality-test switch, rank-test rescue, cell bootstrap, permutation/randomization claim or clinical covariate imputation is allowed.

### Zeros, missingness and numerical evaluation

- Retain B=0 as Y=0 when R is valid and positive. No abundance/detection filter is applied.
- Require all eight donors, complete fixed membership/axes and 0 <= B <= R with positive R. A missing donor, invalid denominator or changed input invalidates the whole endpoint. Do not replace a donor, drop to n=7 or switch sources.
- With exactly one constant arm and positive total SE, use ordinary Welch inference; df is three. With both constant arms, retain a defined point but leave CI/P unavailable. A nonzero constant-group difference does not justify a zero-width population interval.
- Other nonfinite or invalid inference is unavailable, not fabricated P=0/P=1. A positive tail that underflows the native probability representation is unavailable; a finite interval may remain available. Preserve the precise status rather than conflating missing and null.
- Use the fixed exact-integer parser bounds and stable numerical method in the [implementation design](../reports/psc-bach2-naive-rna-method-design.md). No epsilon floor reclassifies tiny variation as zero. Exact rational products evaluate the same mean-log point and omissions without cross-arm cancellation signs. Decimal working precision is computational, not biological precision; native Student-t functions remain finite-precision.
- Any pre-execution numerical correction must preserve this estimand and pass independent artificial-data tests. Once frozen and real values are accessed, do not repair methods, tune thresholds or silently replace pins to obtain a preferred result.

## Fixed diagnostics and disclosure

Retain privately, for every donor: selected-cell counts; B and R; BACH2-positive selected-cell count/fraction; Y; C0/C1 annotation fractions; all relevant source checks; and all eight donor-indexed omissions.

For each omission, recompute only the equal-donor point with 3/4 or 4/3 donors. Do not compute omission P values or intervals, drop a donor from the primary analysis, or treat the overlapping omissions as independent replication.

The public candidate is restricted to the validated aggregate schema:

1. Group n=4, mean Y and sample SD.
2. Delta, SE, actual df, nominal CI/P and availability status.
3. Eight-omission count, complete minimum/maximum, positive/negative/zero counts and strict nonzero sign reversals relative to the full point. A zero full point has no strict reversal; both omission directions remain visible.
4. Aggregate fixed source/selection QC, number of zero-B donors retained and cohort-wide BACH2-positive selected cells. Detection is not an additional genotype-group endpoint.
5. Plan, selection, source-manifest, implementation, runtime, schema and authorization digests.

No donor IDs, source joins, donor expression/denominators, donor-indexed omissions or individual plots are published. Aggregate-only is a publication boundary, not a claim of formal differential privacy or zero possible inference from small groups. Every public nested field must pass strict key/type/finite-value/CI/status checks; allowed outer keys alone do not prevent private data leakage.

Early independent annotation auditing opened a closed cohort-qualification JSON that also contained previously computed donor technical totals. Those totals were not extracted, compared, used for selection or analyzed. No real matrix or technical-margin payload was opened. The reviewer also reports root-authorized opaque hashing of all eight compressed archives without decoding members. Production preparation and validation do not read those archive bodies or the quantity-bearing qualification JSON. These different access scopes remain disclosed; the final comparator uses individual integrity/acquisition receipts. The fixed BACH2 numerator and selected-population RNA denominators have not been measured during method preparation.

## Prospective execution, verification and replay

The root freeze must bind the accepted protocol, author and independent code/tests/designs, source-only selection/plan/schema, closed input receipts, exact runtime, native-output auditor, fixed agreement criteria and acceptance receipts. Native authorization does not cover extra scientific choices. Preserve the outer root seals before and after execution.

Use the native project environment. Pin Python, installed NumPy/SciPy distributions and executable/runtime fingerprints. No dependency installation, new network source or clinical join is part of execution. Private inputs and outputs remain ignored. Default and validate-only commands must not open matrices, calculate endpoint quantities or write result outputs.

Future extraction requires an explicit execute command, a separate exact hash-bound root authorization outside the author's preparation tree, and a named fresh ignored destination. Reject path aliases, traversal, symlinks and reuse. Validate source bytes/axes/selection and stream the complete eight fixed matrices under authorization, retaining no unselected-gene or ADT endpoints.

Write the complete native inventory and last-written completion record. A partial directory, failure marker or CLI success message without an authenticated completion is not a completed result. The only publication candidate initially remains inside the ignored run directory. Root must validate it before any aggregate export.

The separately prepared native-output auditor must be tested and pinned before values. Under a later explicit audit authorization, independently stream all eight complete MatrixMarket sources and reconstruct selected-cell B/R/detection quantities, Y, means/SD, the sole contrast/inference, all omissions and the aggregate projection. Authenticate the exact completed four-file inventory and all three payload byte counts/hashes before parsing numerical results. Existing metadata/synthetic-only checks are not relabeled as real-data verification.

The [frozen verification protocol](../config/psc-bach2-naive-rna-verification-protocol.json) fixes 4,000-digit independent direct-ln/exact-rational point calculations and centered moments. Nonzero points, moments, SD, SE, df and P require relative agreement within `1e-10`; exact zero, sign, constancy and availability must agree, with no absolute epsilon. For interval endpoints, compare `(bound - Delta)/SE` to the independent standardized critical value within relative `1e-9`. The reference tail uses a positive 100-digit Decimal incomplete-beta series with float log-gamma normalization. The reference critical value uses `scipy.special.betainc` inverted with `scipy.optimize.brentq`. Shared `scipy.special.stdtr/stdtrit` corroborate native availability only. The auditor does not call the author's numerical/source helpers. Shared installed libraries remain a limit to software independence; thousands of working digits do not imply biological precision. Preserve these accepted criteria, including bounded iteration/underflow statuses, without post-value relaxation.

The [public aggregate schema](../config/psc-bach2-naive-rna-public-schema.json) and verification protocol are byte-identical copies of their ignored preparation artifacts. Reproduction must restore the exact original ignored paths and other sealed source/selection/acceptance inputs; a clone alone is incomplete. Native root authorizations must be simple named JSON files under `work/psc-bach2-naive-rna-freeze/`, at most 8,192 bytes. Each authenticates the exact plan, selection, source manifest, implementation, runtime, schema, candidate and fresh output destination. A separate per-run audit seal uses that same dedicated namespace, is capped at 16,384 bytes, and additionally binds the audit code/tests/design/protocol and the completed run's actual completion SHA. Completion hashes cannot exist prospectively: only that run-derived binding is filled after completion under the already frozen schema and rules. The [outer root freeze](../config/psc-bach2-naive-rna-execution-freeze.json) requires its pre-effect Git checkpoint to be committed and pushed before either native authorization is used.

Mechanical replay keeps sources, selection, scientific plan, code and runtime unchanged. Separate primary/replay authorizations differ only in destination binding; corresponding authorization/provenance hashes are declared differences, not scientific changes. Require exact corresponding scientific payloads and explain any other mismatch. Verification and replay are computational checks, not independent biological observations.

## Interpretation and stop rule

The published same-cohort broad RNA observation motivates this check. The original genotype labels do not identify a causal nucleotide effect. No REF/ALT, coordinate, strand or risk nucleotide is inferred from the assay code or RNA reference. LD, ancestry, batch, clinical selection and retention remain possible explanations.

The published clinical groups differ in IBD labels, treatment exposure and severity. There is no verified participant-to-archive clinical row key. A slash is not evidence of no disease or no medication; birth year is not age at sampling. Do not merge or impute those covariates. Large cell counts do not resolve the eight-donor confounding limit.

A real RNA shift can qualify an aggregate RNA-similarity reading but cannot refute or prove the separate protein/translation experiments. A null or narrow interval does not establish equivalence without an independently justified margin. There is no mature-miR4464 or named BACH2 ADT endpoint. Pooled-state RNA cannot establish cell-intrinsic regulation, a PSC mechanism, clinical benefit or a treatment recommendation.

After one primary extraction/contrast, fixed diagnostics, independent verification, mechanical replay and complete aggregate reporting, close this branch regardless of direction, P, interval width, sparsity or influence. A failed fixed gate is unavailable/PARK, not a reason to select another gene/state, change normalization, delete donors, recluster, expand the cohort or search for a favorable result.

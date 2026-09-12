# ETS2 validation protocol amendment 1

Recorded September 12, 2026, after methodological review and candidate metadata discovery, before candidate-dataset program scores, differential-expression results or donor target effects were inspected. The original [protocol](ets2-validation-protocol.md) and [rules](../config/ets2-validation-rules.json) are preserved. The [amended rules](../config/ets2-validation-amended-rules.json) are the effective specification for future runs and pin this document, the original freeze and the unchanged source definitions.

This is a prospective local amendment, not a claim that the literature search was blinded or that these choices were registered before candidate metadata were known. It changes the planned intervention analysis; it does not revise the completed liver-atlas analysis or its results.

## Control-only expression matching

Replace the original reference that averaged all conditions with a reference derived only from the appropriate matched controls. Treatment-induced expression must not determine the reference bins used to assess treatment contrasts.

For each dataset, assay, cell population, control background and time stratum, average transformed expression equally across distinct control biological cultures within each donor, then equally across donors who contribute at least one eligible planned pair in that stratum. A shared control enters once, irrespective of the number of drug doses or stimulation contrasts using it. Controls in different backgrounds are not pooled: for example, a MEK inhibitor under LPS stimulation is compared with the matched LPS-plus-vehicle control; a generic LPS experiment is compared with its unstimulated control. Freeze these strata and control sample IDs with the sample map before scoring.

All related contrasts within the same stratum reuse the same reference bins and matched gene sets. Only donors with both members of a given contrast contribute to that contrast; a missing treatment arm does not create a replacement pair. Report reference donors separately from each contrast's paired donors. The 20 bins, 1,000 draws, deterministic seeds, program exclusions and score formula remain unchanged.

## Biological relevance and assay readiness

Keep separate fields for biological evidence tier, assay readiness, donor-design status, independence status, target-engagement status and feature-universe status. A primary human MEK-inhibitor microarray study remains biological tier B even when it is not ready for the raw-count procedure. Arrays, normalized RNA measurements or other quantities require a separate prospectively fixed assay procedure; they are not recoded as integer counts.

Selection within an assay-ready tier retains the original order: complete biological donor pairs, documented freedom from batch/assay confounding, newer publication, then lexical accession. A readily available tier C, D or R dataset can support explicitly labeled calibration while the primary-human direct-ETS2 question remains open. It does not become tier A because it is easier to analyze.

For a direct genetic perturbation, preserve every guide, construct, control, dose and time arm before scoring, including relevant author-reported exclusions. Record published evidence for perturbation fidelity or target engagement and how it was measured. A construct label alone is not evidence of successful target modulation. Do not select guides, donors or conditions based on their observed program response or their ETS2 RNA change. Missing independent target-engagement evidence limits a null or positive result and remains explicit. Pharmacological studies similarly record the compound, background, matched vehicle and available pathway-engagement evidence without assuming ETS2 mediation.

## Feature universe, duplicate policy and structural QC

The unchanged CPM procedure requires the source's full exported gene-count universe before expression-abundance, differential-expression or signature filtering. Preserve measured zeros and all source-counted gene features in library totals; do not restrict totals to mapped program genes. Source-defined gene annotation choices are documented separately and are not themselves expression filtering. If prior abundance filtering or the denominator's provenance is unresolved, mark the unchanged all-feature procedure pending rather than silently assuming it is satisfied.

Record the source's quantification method, identifier namespace, annotation version where available, matrix dimensions and whether the measurement is raw gene counts. Do not infer raw counts merely from a filename. Library summary rows or transcript-level quantities require a documented source-specific handling rule before readiness. No normalized values are converted back to counts.

Each sample column must map exactly to a declared biological sample, donor, condition and assay. Require unique sample identifiers, finite nonnegative integer gene counts and positive full-universe library totals. Missing counts are not zeros. Record all missing and extra matrix/sample-map entries. Unresolved duplicated feature IDs, including collisions after numeric Ensembl-version removal, stop the affected raw-count analysis until their source meaning is resolved. Do not sum, average or choose one duplicated feature based on expression. Distinct stable gene IDs with the same symbol remain separate in the count universe; symbol-only program mapping is unavailable for that ambiguous symbol. All unmapped and ambiguous source program members remain in the coverage audit.

Preserve author-reported sample QC, technical splits and exclusions with their reasons. Sum only documented technical sequencing splits of one biological sample. For separate biological cultures, normalize each culture independently, then average transformed expression equally within donor/condition. No additional sample exclusion is selected from a program score, expected effect direction or an invented library-size cutoff. Structural failures and source QC exclusions make pairs unavailable; non-fatal library-size differences are reported. Any new QC threshold requires a dated amendment before dependent effects are inspected.

## Single-cell culture accounting

If a later selected source is single-cell, form raw pseudobulk counts for each distinct biological culture/sample, donor, condition, assay and author-defined cell population first. Combine only documented technical library splits. Apply the 20-cell primary gate and 50-cell sensitivity independently to each culture in each condition, normalize eligible culture pseudobulks separately, and average their transformed expression equally within donor/condition. Each member of a donor-paired contrast must independently pass its gate; cells in one condition cannot make another condition eligible. Report culture counts as well as donor counts and retain at least three complete donors for a paired summary. Freeze the culture and cell-population mapping before any such run. This resolves the original ambiguity between culture averaging and donor-wide count summation.

## Search boundary and provenance request

The initial phrase “bounded search” did not specify database queries or page limits in advance. The final search reports must retain the exact sources, queries, retrieval date, caps or pagination, candidate exclusions and actual stopping boundary. Do not call the inventory exhaustive or retrospectively preregistered. Finish the current primary-source metadata search and the already selected small-file checks; a future expanded search is a new recorded round.

The researcher draft distinguishes recomputed raw ES from compared archived NES. It requests the pre-symbol-collapse per-gene inputs, contrast definitions, filtering/mapping/duplicate/tie/export rules, row counts and software versions that could resolve the historical rank-export gap. It is an unsent draft and author feedback is not a prerequisite for public-data preparation.

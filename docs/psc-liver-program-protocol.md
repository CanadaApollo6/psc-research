# Pediatric liver: fixed-program tissue-context protocol

## Question and prospective boundary

This specifies one finite analysis of GSE303271 before any new liver count normalization, program scores or effects. A separately hash-bound root freeze and explicit execution command are required. Source-only qualification, identifier mapping and synthetic tests are not execution permission.

**Question:** how do fixed ETS2 perturbation-derived, comparator and unsigned extracellular-matrix transcript-abundance composites differ among the deposited PSC, ASC and AIH groups?

The focal comparison is **ETS2 gRNA1-down, PSC minus AIH**. All ten Core programs and all three contrasts remain in one complete 30-hypothesis family. The other results cannot replace the focal comparison after effects.

This is a narrower tissue-context question. The epithelial IL17A-alone response-transfer endpoint remains **HOLD and outside this family**, with no substitute genes or zero score. Original program preparation preceded MacroMap effects, but this study is not independent mechanistic confirmation of MacroMap. It tests whole-biopsy observational abundance, not ETS2 activity, cell-intrinsic regulation, cytokine response, causal disease mechanisms or clinical benefit.

## Fixed source and independent-unit assumption

- Use the complete qualified **41,096 released rows and 64 source columns** from GSE303271. Keep every row in the normalization universe, including zeros and unresolved identifiers.
- Keep deposited diagnosis groups separate: **17 PSC, 17 ASC, 30 AIH**. Do not reinterpret the article's combined PSC wording as permission to merge PSC and ASC.
- The units are source patient-samples. Accept the paper's patient-independence assertion as an explicit working assumption, not an independently verified person/visit crosswalk. Distinct accessions alone do not prove independent people. Known unresolved repeats would block execution.
- The assay is whole diagnostic liver biopsy, not sorted cells. Cell composition, fibrosis, disease stage, inflammation, medications, age, sex and technical batch can affect scores. No sample-linked covariate adjustment is available or invented.
- The source describes GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4. Historical GTF and counting settings remain unresolved. Current identifier reconciliation is not recovery of that annotation.
- Use the pinned released raw gzip, accepted float64 count cache, source feature order and sample join. Never repair, round, rescale, drop or recode source values.
- GSE159676 remains expression HOLD. The NoPSC atlas remains PARK. Neither is a rescue or replication cohort for this protocol.

Source URLs, retrieval records and transformations are in the closed [input qualification](../reports/psc-liver-input-qualification.md), [identifier qualification](../reports/psc-liver-program-inputs.md) and [external-source qualification](../reports/psc-liver-external-program-qualification.md).

## Fixed programs and identifier tiers

The exact original memberships and exclusions are those of the accepted identifier preparation. Published labels are retained. All ten executed scores use **positive, equal mapped-member weights**. In particular, do not invert gRNA1-up.

| Fixed program ID | Fixed intended coverage denominator | Core mapped | Strict mapped |
|---|---:|---:|---:|
| `ets2_g1_dn` | 927 | 848 | 825 |
| `ets2_g1_up` | 668 | 612 | 536 |
| `ets2_g2_dn` | 875 | 809 | 789 |
| `chr21_dn` | 141 | 126 | 122 |
| `inflammation` | 815 | 750 | 744 |
| `interferon_gamma` | 139 | 121 | 118 |
| `oxidative_stress` | 436 | 403 | 400 |
| `apoptotic_signaling` | 588 | 544 | 542 |
| `ets2_g1_dn_without_comparators` | Parent 927 gate, then at least 10 retained | 695 | 673 |
| `reactome_ecm_organization` | 321 official source symbols | 301 | 299 |

Core `exact_symbol_ensembl_bijection` is primary. Strict `plus_reciprocal_entrez` is an identity sensitivity only. Use the globally qualified Ensembl 116 relationships and the existing complete ambiguity audit, not a measured-only Entrez bijection. Do not rescue aliases, collapse ambiguous genes, use an alternate favorable namespace or zero-fill missing members.

For each ordinary program require `5*mapped >= 4*intended` and at least ten unique mapped members. Nonoverlap removes only the original four-comparator union: inflammation, interferon gamma, oxidative stress and apoptotic signaling. It retains its parent gate and a separate ten-member minimum. ECM is not added to that subtraction.

ECM is the complete official MSigDB M610 / Reactome R-HSA-1474244 set, MSigDB 2026.1.Hs / Reactome 95. Its denominator is 321 official symbols, not the 340 preserved source Ensembl records. It is unsigned constitutive organization membership, not a directional fibrosis-response signature. Literal ETS2 is absent; no new target exclusion was made. The qualified manifest, complete CSV and source-worker handoff must all authenticate. The held epithelial IL17 slot is neither a member of this panel nor an empty endpoint.

All twenty tier/program specifications, mapping statuses and source references are pinned before effects. Their canonical digest is `17a291f47fbb0ec0e815ac2f592ba4d0f657f3f5e3cbe4689dae5ecc44c90698`. This canonical digest differs from the serialized-file SHA and must not be substituted for it.

## Normalization and estimand

For each source column, sum counts over **all 41,096 released rows**. Compute `log2(1 + 1000000 * count / column_total)` (the stable implementation is `log1p(CPM)/log(2)`). A program score is the arithmetic mean of these values across its fixed unique mapped genes. No matching, learned weights, variance standardization, expression filtering or historical-universe reconstruction is used.

Counts must have the frozen shape, be finite, nonnegative and exactly integral. Each count and exact column total must be at most `2**53`, and every total must be positive. Bounded integer summation checks totals before float rounding can hide an excess. Invalid data stop the run; they are not an unavailable biological endpoint.

The contrast is the equal-source-patient-sample mean score difference. It is not a gene-level fold change, a ratio of arithmetic gene abundance, a TF/pathway activation measure or an absolute RNA amount.

Ordered contrasts:

1. **PSC minus AIH**, primary; focal ETS2 gRNA1-down.
2. PSC minus ASC, fixed secondary.
3. ASC minus AIH, fixed secondary.

These contrasts have rank two, but all three remain in the fixed family. Cell types and genes are not independent biological replicates.

## Inference, numerical limits and full family

Use ordinary unequal-variance Welch inference on source-sample scores, sample variances (`ddof=1`), the actual Welch–Satterthwaite degrees of freedom, a two-sided t P value and a pointwise 95% t interval. Do not substitute a normal critical value or select a test from observed distribution checks.

Detect exactly constant vectors before calculating their means and variances. One zero-variance arm is valid ordinary Welch when total standard error is positive; flag it and retain the actual df. If both arms have zero variance, retain a finite descriptive difference but leave P, df and CI unavailable. Preserve real small variation; use no epsilon floor or variance filter. Scaled variance contributions avoid avoidable df underflow. Nonconstant variance underflow, nonfinite inference or unrepresentable positive tail probabilities have explicit unavailable statuses, not fabricated certainty or observed P=0/P=1.

Invalid count/score inputs fail data integrity. Inference-only numeric failure does not delete a finite descriptive effect or its complete hypothesis row. A naturally valid P=1 at t=0 is distinct from unavailable inference.

The sole discovery family is **all 30 Core program/contrast hypotheses**, including unavailable hypotheses. Report BH and BY, without selecting the more favorable adjustment. For unavailable tests retain raw P and displayed q as null; use P=1 **only inside the complete correction vector**, with explicit substitution flags and private numeric placeholders. Strict has the full 30-row diagnostic grid but no discovery q values. BH depends on valid tests and its dependence conditions; BY handles arbitrary dependence of valid marginal tests. Neither fixes confounding, small-sample calibration or wrong independent-unit assumptions.

## Fixed bootstrap and omission diagnostics

- Generate **10,000** within-diagnosis patient-sample bootstrap resamples with replacement, **PCG64 seed 2026091602**.
- Generate groups in order PSC, ASC, AIH, using original source-column order within each group. Reuse the same draws across every program, tier and contrast.
- Scores, source totals and memberships stay fixed. Use linear 2.5/97.5 percentiles for pointwise sensitivity intervals. No bootstrap P values, new correction family or reselected primary endpoint.
- Keep constant draws and flag degenerate intervals. A degenerate empirical interval does not establish zero population uncertainty. These intervals do not include annotation, clinical-confounding or unknown-person-dependence uncertainty.
- Omit each of the 64 source units for each of the 60 rows: **3,840 private omission records**. Recompute Welch summaries from the remaining fixed scores only. Do not renormalize, remap, change gates or remove a unit from the primary result.
- Out-of-contrast units retain explicit unchanged private rows. Public stability summaries use only the involved units: **47 for PSC−AIH, 34 for PSC−ASC and 47 for ASC−AIH**. All ranges, defined-effect counts, changes and reversal counts use this relevant subset, not 64.
- Reversal means strictly opposite nonzero signs. Missing effects are not evidence of zero reversals. Omission P values remain diagnostics, with no corrected discovery claims or independent-replication interpretation.

## Execution and verification contract

The accepted implementation is `scripts/analyze_psc_liver_programs.py`, with the source-blinded [method design](../reports/psc-liver-program-method-design.md) and [independent review](../reports/psc-liver-program-independent-review.md). The separate native-output verification adapter must be ready, tested and pinned before the first real execution. It must not relabel real results as synthetic fixtures or alter the already accepted independent math.

The root freeze will bind this protocol, exact implementation/tests, independent verifier/adapter/tests/design, accepted review receipts, native versions, source/metadata/membership pins and compiled specification digest. The native method requires Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. Root keeps `PYTHONHASHSEED=0` and BLAS/OpenMP thread variables at 1. No package installation, network acquisition or external model inference is required.

Execution requires the designated `config/psc-liver-program-freeze.json`, its exact SHA-256, state `frozen`, an explicit root role/reference/time and real-execution permission, accepted matching review, fixed method/source pins, and the acknowledgement `ROOT_AUTHORIZED_REAL_EXPRESSION_EXECUTION`. The approved primary output is a fresh ignored `work/psc-liver-program-plan/runs/primary-v1`. Mechanical replay changes only the output argument to a fresh `work/psc-liver-program-plan/runs/replay-v1`; it uses the same scientific plan and is not a sensitivity analysis.

Reject reused destinations, aliases, traversal, symlinks and source overwrites. Root validates all sealed artifacts before and after execution. The runner authenticates raw/count-cache hashes after authorization and before numeric loading. A pre-execution refusal creates no result directory. Later failure leaves an explicit failed-run receipt and unavailable status inventory; a directory without a valid completion receipt is not a completed analysis. Workflow JSON and file checks are not a security sandbox against a malicious local actor.

Independent real-output verification follows the [pre-effect adapter design](../reports/psc-liver-program-execution-verification-design.md). It authenticates actual inventories and hashes, reconstructs all 2,630,144 raw source count cells and axes, and checks 64 totals, all 1,280 score values, 60 endpoint rows, the complete Core correction, all 640,000 shared bootstrap indices and 600,000 differences, intervals, all 3,840 omissions and public 47/34/47 summaries. Counts, totals, axes, membership, inventory, statuses and RNG indices must agree exactly. Positive P/q, variance, SE and df use relative tolerance 1e-10 without an absolute floor. Other finite floating fields use relative 1e-10 and absolute 1e-12. Exact-zero, availability and probability/interval domain checks precede tolerances. The [sealed verification contract](../config/psc-liver-program-verification-contract.json) is a byte-identical public configuration copy of the ignored contract expected by the pinned adapter; it contains no individual scientific values. Reconstruction must restore that exact copy at `work/psc-liver-program-independent-review/execution-adapter/production-contract.json` alongside the separately required ignored source/acceptance artifacts. It does not make a clone alone a complete input bundle.

Native log-CPM intermediates are not saved; do not claim a direct comparison of every such native intermediate. Distribution functions and NumPy's generator remain shared software dependencies, not independently reimplemented mathematical libraries. The real audit needs its own explicit root acknowledgment and externally supplied exact adapter, plan and completion hashes. Exact audit scope, tolerances and any unavailable checks must be reported.

A separate clean mechanical replay compares actual complete outputs, including private artifacts. Any byte or numerical difference must be explained rather than hidden by overwriting receipts or rehashing altered scientific inputs. Publication is gated on independent verification and faithful reporting, not on effect direction or significance.

## Privacy, reporting and finite stop rule

Keep source-unit IDs, joins, counts, normalized expression, scores, library totals, individual omissions, bootstrap arrays and reconstruction caches in ignored directories. Publish aggregate-whitelisted group/program/contrast summaries, complete status grids and small verification/replay receipts only. Preserve all missing, opposing, null, unstable and numerically unavailable results. Do not disclose individual source values in plots or tables.

Report all ten programs in both tiers and all three contrasts. Separate raw, adjusted, sensitivity and descriptive quantities. Pointwise intervals are not simultaneous. Annotation agreement and numerical replay are not independent biological replication. No result licenses a clinical recommendation or an ETS2-specific causal claim.

After one frozen primary analysis, its fixed verification/replay and complete reporting, close this liver branch. Failed inputs or weak/null/positive/fragile results do not authorize a different program, source, gene exclusion, normalization, covariate guess or method. Reopening a held source or proposing a new experiment needs a separate question and authorization. The unfulfilled epithelial IL17-response-transfer question remains visible regardless of this analysis's result.

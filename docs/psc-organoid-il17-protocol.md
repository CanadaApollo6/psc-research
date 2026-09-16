# PSC cholangiocyte organoids: prospective IL-17A interaction protocol

**Status: frozen before any new effects at 2026-09-16T00:41:14.393412+00:00.** The [executable freeze](../config/psc-organoid-il17-execution-freeze.json) pins the qualified source bytes, population, feature universe, methods, code and software. Native structural validation must pass before inference. The 84 source/analysis tests and 24 independent-verifier tests pass. Published source narratives are already known. This is public-data reanalysis, not untouched or externally registered validation.

## Question and biological units

For GSE239283, estimate the disease-associated difference in modeled IL-17A log2 fold responses. The primary NB interaction is a **log2 ratio-of-ratios**, not a difference of raw counts:

`log2[(normalized mean_IL17 / normalized mean_CNT)_PSC / (normalized mean_IL17 / normalized mean_CNT)_nonPSC]`.

The paired NB model uses donor-specific baselines and a common response within each disease group. Its fitted interaction need not equal an unweighted mean of donor log changes. The prespecified Welch sensitivity instead estimates the equal-donor difference in mean `log2(CPM+1)` changes. These are related but different estimands; neither is a cell-level treatment effect.

There are four PSC and four non-PSC source donors. Each supplies one CNT and one IL17 culture library: **8 pairs, 16 libraries**, not 16 donors or 16 pairs. The 46,343 released cells are observations within these libraries, not independent biological replicates. The published culture exposure is 100 ng/mL IL-17A for 24 hours at passage 3.

Primary population: **all published retained organoid cells**, without a new treatment- or effect-selected cluster. This is a source-population definition, not a new proof that each cell has the same lineage/state. Author-defined clusters may be secondary only under a separate fixed specification. The non-PSC group comprises procedure controls, not healthy volunteers. Sampling route, age, clinical history and culture selection may confound disease-associated response differences.

## Input-quantity and identity gate

Use original raw RNA UMI/count measurements only for the proposed count likelihood. Check every file, sparse coordinate, feature/barcode axis, library and donor-treatment join. Integer validity alone is insufficient: SCT-corrected counts can also be integers. Compare merged per-cell quantities with the original RNA and SCT metadata, and preserve any discrepancy.

If needed, qualify original per-library processed RNA matrices and select the **same published retained cells** by exact source library plus barcode. Record raw features not present in the merged export, missing retained cells, duplicates and source-label transformations. Do not guess a historical gene filter merely to force a metadata total match. Define full-library sensitivity denominators as the complete qualified released **raw** feature universe, not an asserted reconstruction of every historical reference feature.

### Pre-effect raw-source scope decision — September 16, 2026

The original 16 Cell Ranger processed-library exports contain 36,601 source Ensembl features and 48,164 cells. Exact source-library/barcode joins recover all 46,343 published retained cells; 1,821 original-export cells are outside this fixed population. The original retained-cell counts total 1,128,543,232, versus 1,128,472,969 in historical RNA metadata. They also contain 68,460 more positive gene–cell entries. The original raw counts are therefore **not an exact reconstruction of the historically filtered RNA assay**. The difference alone does not identify the historical filter/export command, and no guessed filter will be applied.

Subject to final source-provenance and independent aggregation checks, use the **complete original exported raw feature universe for the fixed published cells**. This is an explicit new released-raw-universe analysis. It is scientifically distinct from feeding the merged SCT-assay-consistent matrix to an NB model, and from claiming recovery of every historical processing step. Full-source raw totals define the CPM sensitivity denominator. Keep the merged SCT comparison and historical RNA discrepancy in the evidence record.

Preserve stable feature IDs and as-published labels. The 36,601 IDs have 36,591 unique source gene labels; repeated labels are not a reason to collapse their counts. Do not merge genes by guessed symbols, collapse ambiguous IDs or choose representatives by P value. Exact source bytes, quantity/provenance evidence, selection map, method/code/software hashes and freeze time will be pinned in a separate execution plan. No real effects are permitted while that plan is draft.

## Proposed primary inference

- Sum raw counts within each of the 16 donor-treatment libraries. The count matrix has libraries, not cells, as columns/observations.
- Fixed eligibility: at least 10 counts in at least 4 of the 16 libraries. Keep every source row and exclusion reason. Use the same eligible family for all specified sensitivities.
- PyDESeq2 0.5.4 paired NB model: intercept, seven donor indicators, IL17 treatment, and `PSC × IL17`. Omit a redundant disease main effect, which is collinear with donor indicators. Check rank and leverage.
- The **interaction** is primary. Separately retain the treatment effect in non-PSC donors and the combined treatment effect in PSC donors, with the correct coefficient covariance. Their separate within-group significance is not evidence of an interaction.
- Require ratio normalization with a nonempty all-positive gene basis. No iterative normalization fallback. Request parametric dispersion; its documented built-in mean fallback is allowed and recorded. No count replacement, independent filtering or LFC shrinkage.
- Donor-plus-condition rows are not replicated identical design rows. Do not silently apply a Cook exclusion rule with no eligible replicates. Set automatic Cook filtering off; retain all-sample Cook and convergence diagnostics. No effect-based donor/sample exclusion.
- **Pre-effect Wald boundary safeguard:** a within-group treatment response requires positive pooled counts in both that group's CNT and IL17 arms. The interaction requires positive pooled counts in all four disease-by-treatment arms. An entirely zero arm can put a log-response MLE on a boundary; an entirely zero disease group can leave its response unidentified. A finite package/ridge result does not resolve this. Preserve all original package Wald fields as diagnostics, but mark affected contrast inference unavailable. A single all-zero donor pair is a diagnostic, not an automatic whole-group exclusion when the common response remains identifiable from other donors. This rule was specified using mathematical/synthetic cases before real effects.
- Use BH over the entire fixed eligible gene family, retaining missing P values and using 1 only as the correction input. The Wald safeguard does not shrink this family. Report BY sensitivity for the primary interaction. Clearly label small-donor-number model/calibration assumptions; model output is not proof of reliable asymptotic inference.

## Prespecified sensitivities and descriptive checks

1. Transform each library to `log2(1 + CPM)` using the complete unfiltered qualified raw released-feature totals. Compute one IL17-minus-CNT change for each donor. Compare the four-versus-four changes with Welch's test and its actual degrees of freedom. Correct the fixed eligible family; retain unavailable tests.
2. Leave one **donor pair** out at a time and report interaction-effect/sign diagnostics from the paired changes. These overlapping fits are not new replications and do not redefine eligibility.
3. Enumerate all **70** four-versus-four allocations of the eight paired changes. Report the two-sided absolute-mean-difference allocation-tail fraction as a **descriptive exchangeability calibration**, not a new FDR discovery criterion or disease-randomization test. Include the observed allocation and its complement; no Monte Carlo +1 formula. The minimum two-sided fraction is 2/70. Predeclare numeric tie handling in code/tests. Disease and donor history were not randomized, so exchangeability is an assumption.
4. Retain all donor-level changes in ignored reproducible work outputs and show aggregate/unlabelled donor plots where useful. Do not publish personal clinical metadata or source-to-person joins. Preserve missing, opposite and unstable directions.

Any focused program comparison requires separately pinned gene sets/mapping and an explicit testing family before its scores. No gene or program is selected from the new effects to create a confirmatory claim. Counts and model outputs are relative expression measurements, not absolute transcription rates or patient outcomes.

## Interpretation and completion

A treatment contrast can support an ex vivo response under the culture-design assumptions. The disease-by-treatment difference remains partly observational. It cannot show that IL-17 causes PSC or that an IL-17 intervention benefits patients. More cells do not repair the eight-donor sample size.

Complete all eligible-gene contrasts and sensitivities, coverage/quantity accounting, independent numerical checks and a clean-input replay before reporting a result. Stop the NB branch rather than feeding it an unqualified transformed assay; a different descriptive estimand would require a new explicit plan. Retain the published-response prior and do not claim discovery merely because a familiar pathway is reproduced.

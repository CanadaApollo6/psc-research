# Frozen MacroMap program-context protocol

## Status and question

This is a separate pre-effect protocol for `macromap-fixed-program-context-v1`. Its exact freeze time, source/preparation/code hashes and numerical contract are in [the execution freeze](../config/macromap-program-context-execution-freeze.json). The [method specification](../reports/macromap-program-design.md) and [independent review](../reports/macromap-independent-review.md) are incorporated by their frozen hashes. Preparation and synthetic testing were completed before real normalization, expression matching, scoring or model fitting. Earlier publications and the original ETS2 results were already known; this is not a claim of publication-level blinding.

**Question:** does the unchanged ETS2 gRNA1-down program, excluding ETS2, add held-out stimulus-context information beyond inflammation, interferon-gamma, oxidative-stress and apoptotic-signaling programs? MacroMap contains healthy-derived iPSC macrophages, not a PSC diagnostic or treatment cohort. The analysis proceeds independently of the completed organoid outcome.

## Explicit source and quantity decisions

- Use all **58,243 released count features ×4,698 samples** and the pinned original-study GENCODE v27 reference. All 273,625,614 original tokens passed the exact nonnegative-integer check, with positive library totals below 2^53. No raw reads or eQTL archive are needed.
- Authorize a **new full-released-universe** `log2(1 + 10^6 * count / released-library-total)` quantity. The 45 reference-only `_PAR_Y` features stay absent, not zero or collapsed. This does not reconstruct an unobserved historical universe or the original filtered 14,060-gene pipeline, and it does not change an earlier historical CPM gate.
- Accept the complete, hash-pinned original `source_value` exports, including old unmapped members, with agreement across modality exports and source-list audits. The original ETS2 ZIP remains absent. This is verified-export reuse, not a new archive read or an atlas-filtered membership subset.
- Use fresh v3 preparation with summary SHA-256 `eae76811c5b46388965f4f2c27624e650748b9d7adc216f75bff05b389e1f507`. It and all nine realized design artifacts replay identically. An earlier unfrozen preparation summary changed during development; the independent review preserves its initial bytes separately. No real source or scientific effect was changed by that event.
- There are 205 mapped source lines and four separate `NOMATCH` lines across 57 runs. Published unrelated-donor provenance is not an independent kinship check. The primary prediction population has **185 mapped both-time lines: 114 Mod_SmartSeq2 and 71 NEB**. Lines, not repeated stimulation rows, are the analysis units.

## Fixed analysis

Retain the complete nine-program response panel for nine strict stimuli at both 6 h and 24 h. The complete registry also retains PIC at each time as `unavailable_mock_identity`; basal control does not identify its lipofection effect. Pair treatment with its own line's same-time basal control. Do not substitute precursor samples or infer unavailable measurements.

Whole-run five-fold holdouts are fixed within protocol, shared across times and conditions, using the stated hash order and seed 20260916. Each outer fold learns abundance bins and matched sets from training-run, same-time basal controls only, once per line. Use all released features for normalization/binning, the fixed eight-program-plus-ETS2 background exclusion, 20 bins and 1,000 no-replacement matched sets. Preserve mapping/pool failures. Both training and held-out predictor matrices use that same outer fold's reference; a global table of each training line's own-fold response is not substituted.

The primary model is separate nine-class symmetric ridge multinomial regression per protocol/time/outer fold. Compare four fixed baseline predictors with the same four plus `ets2_g1_dn`. Weighted training loss gives each line equal weight and each available stimulus within that line equal weight. Learn predictor means and population scales on training paired changes only. Use lambda 0.1, the frozen optimizer/tolerances, no tuning and all nine original classes. Constant training features use scale one and retain their flags. No class, fold, time or failed fit may be dropped to rescue a result.

The single primary estimate is the baseline-minus-extended held-out **natural-log loss**: average available stimuli within each line/time, average the two times once with weights 1/2, then average the 185 lines. Equivalently, pool protocol estimates with fixed shares **114/185 and 71/185**. Report all four protocol/time cells, folds and failure statuses before the pooled number. Missing expected predictions make the endpoint unavailable, not a partial favorable estimate.

## Uncertainty and fixed sensitivities

Use 5,000 within-protocol run-cluster bootstrap draws, PCG64 seed 2026091601, lexical run order and linear 2.5th/97.5th quantiles. Keep every line, condition and time of a sampled run together, retain original protocol shares, and do not count times or folds as independent replicates. The primary and four-cell gain intervals hold OOF fits, scaling, references and folds fixed. They are **conditional fixed-prediction intervals**, not full-procedure/refitting or unconditional generalization confidence intervals. No unconditional superiority P value is proposed.

Protocol-specific response means/medians have pointwise run-bootstrap intervals. **Pooled response means and medians are descriptive only: no pooled response confidence intervals.** Preserve leave-line/run-out ranges, unadjusted and background components, signs and full ranges. No response-P screen or simultaneous full-panel coverage is claimed.

The fixed sensitivities are: (1) remove the mapped comparator union from the ETS2 predictor, retaining its parent gate, baseline, folds, lambda and cohort; (2) separately admit the four individually labeled `NOMATCH` lines, rebuilding training-only references/scales with the same original run assignment and metadata-fixed shares. Neither replaces the primary. Leave-run-out prediction diagnostics keep fixed predictions and original protocol shares; they are not refits.

A positive gain supports only limited incremental context information relative to these comparators. A weak/opposing estimate, protocol disagreement or loss after overlap removal limits utility. An interval crossing zero does not prove equivalence. No result proves ETS2 activity, mediation, PSC specificity, genetic allele direction or clinical benefit. Finish the complete fixed panel and stop rather than change programs, stimuli, times or regularization.

## Execution, verification and publication

Use the pinned `.venv` with Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. Set `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1` and `PYTHONHASHSEED=0`. No new dependency, GPU, service access or paid resource is authorized.

The guarded runner validates the external plan hash and refuses execution without `--execute`. It checks source/code/runtime/preparation pins, creates a fresh literal ignored destination, reserves the cache namespace exclusively through no-replacement publication, checks the raw-source pin before any matching/fit and rechecks pins/cache before completion. A failure gets a retained failed-execution receipt. Do not edit a frozen method to make an observed failure disappear; a necessary change needs an explicit amendment.

All normalized matrices, source/donor joins, line responses, per-fit arrays and model probabilities stay ignored under `work/macromap-design/`. Only aggregate outputs may be promoted to Git. Keep all expected 80 paired-model audit cells (160 small model fits across scopes/variants), even when unavailable.

The frozen independent verifier checks all fit cells, expected predictions, aggregate points and predictive bootstrap intervals. Its predeclared source-to-score audit covers fold 0 for both protocols/times/scopes, all nine programs and four hash-selected training plus four held-out pairs per context. Eight hash-selected source samples have full-universe normalization reconstruction. All response-panel points/statuses are checked; numerical response-interval replay is fixed to CIL at 6 h, both protocols, all programs/quantities/scopes. These are bounded computational checks, not biological replication or exhaustive independent raw-to-score reconstruction.

A clean replay may use a new plan that changes **only the ignored output directory** and records its link to the original freeze separately. Every scientific field, freeze timestamp, input/preparation/code hash, method and software setting must remain identical. It must not overwrite the primary run. Compare numerical payloads and report run-time/path receipt differences explicitly; do not expect timestamped execution receipts to be byte-identical.

# PSC organoid IL-17A method design

**Status:** implemented and tested on synthetic data only. This note is not the execution freeze and contains no real organoid effects. The root analysis plan must pin the final qualified inputs, code, software and UTC freeze before inference.

## Scientific scope and source gate

GSE239283 has four PSC donors (`PSC_4`, `PSC_6`, `PSC_8`, `PSC_9`) and four non-PSC procedure-control donors (`nonPSC_2`, `nonPSC_3`, `nonPSC_4`, `nonPSC_5`). Each contributes one `CNT` and one `IL17` culture. IL-17A exposure is 100 ng/mL for 24 h. The independent response units are **eight donors**, not 16 libraries or 46,343 cells. The controls are not healthy volunteers. Published IL-17 responses were known before this analysis; this is not untouched validation.

The primary population contains **all 46,343 published retained organoid cells**. This is a source-defined population, not an independent claim of identical epithelial states. There is no treatment-selected cluster analysis or automatic donor/cell exclusion.

The [completed source qualification](psc-organoid-input-qualification.md) establishes that every merged-matrix cell's total and positive-feature count matches `nCount_SCT` and `nFeature_SCT`. The 21,562-feature integer matrix is therefore SCT-assay-consistent, **not usable as unnormalized RNA input to the NB model**. This does not recover its exact historical export command.

The original per-library processed Cell Ranger count exports provide **36,601 unique Ensembl feature IDs** and exact matches for all 46,343 retained cells. The root accepted that full original released raw-count universe, not the 21,562-feature SCT universe. No extra source cells enter the analysis. Its retained-cell total exceeds the historical `nCount_RNA` metadata by 70,263 counts; this unresolved feature/export difference remains in the qualification record. No feature list was guessed to force historical margin agreement. The final plan must pin these original inputs and 36,601 expected features before inference. The script leaves feature count configurable only through that source-qualified freeze.

No Ensembl/symbol inference, alias substitution, version stripping, duplicate collapse or significance-based identifier choice occurs here. All original feature columns are retained. Source feature IDs must be unique before exact joins. Unequal source universes and source assay identity belong in the acquisition qualification, not in a guessed repair by the model code. "Complete" always refers to the specifically pinned released input universe, not an unverified historical full library.

## Primary paired count model

1. Sum qualified raw counts for every retained cell within its source donor × culture library. The acquisition code owns the exact cell, feature and library joins.
2. Retain a feature if it has **at least 10 counts in at least four of the 16 libraries**. Keep this family fixed across all methods and leaveouts. Retain all other features in public tables with `excluded_fixed_count_filter` and no invented effect or P value.
3. Use **PyDESeq2 0.5.4**. The explicit numeric design is an intercept, seven donor indicators, `treatment_IL17`, and `disease_PSC:treatment_IL17`. A disease main effect is redundant with donor indicators and is not included. The complete design has rank 10, residual df 6 and unweighted leverage 0.625 at each library. The code recomputes these diagnostics rather than assuming them.
4. Use `size_factors_fit_type="ratio"`. The normalization basis is the fixed eligible features with positive counts in every library. Check that this basis is nonempty before calling the library. Independently reproduce `exp(median(log(count/geometric_mean)))`, including the even-feature convention. An iterative/poscounts fallback is prohibited. This normalization assumes a suitable reference distribution of gene changes; a large common compositional response can violate that assumption.
5. Fit a parametric dispersion trend. The **built-in mean-trend fallback** is allowed and explicitly reported. Set `refit_cooks=False`, `independent_filter=False`, `cooks_filter=False`. Do not shrink LFCs or replace counts.
6. The primary two-sided NB Wald contrast is the interaction. Its positive sign means a larger model-based IL-17 log2 fold response in PSC than in non-PSC. Also export the non-PSC treatment coefficient and the PSC treatment contrast (treatment + interaction). The latter uses both variances **and twice their covariance**. All three contrasts share one fit and are not independent replications.
7. Use BH over the whole fixed eligible family. Keep source/model P values unavailable when unavailable; replace them with 1 **only in the correction vector**. The resulting q value for such a row is 1, and its test status stays unavailable. Report BY over the same family as a dependence sensitivity. Neither correction repairs small-sample test miscalibration. Each within-group contrast has its own secondary family; it does not enter the primary interaction discovery rule.

### Necessary count support for a finite log-response contrast

A narrow synthetic check was added before the real-data freeze. With every PSC count zero for an otherwise eligible gene, PyDESeq2 still returned an interaction of −0.104 and P=0.965. With the entire PSC IL17 arm zero, it returned −14.924 and P=6.23×10⁻¹⁹. These are **synthetic package outputs**, not organoid findings. Ridge penalties, mean clamps and numerical convergence do not create an identified finite response: a wholly unexpressed disease group has an unidentified response; one entirely zero arm places its unconstrained log-response at an infinite boundary.

The root approved this narrow prespecified support gate:

- A within-group NB treatment response requires positive pooled raw counts in both that group's CNT and IL17 arms.
- The interaction requires positive pooled counts in all four disease × treatment arms.
- On a failed gate, main `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and confidence limits are unavailable. The feature remains eligible, with correction-only P=1 and q=1 in the unchanged family. The other disease's supported within-group contrast can remain available.
- Native package values remain separately labeled `pydeseq2_unmasked_log2FoldChange`, `pydeseq2_unmasked_lfcSE`, `pydeseq2_unmasked_stat`, `pydeseq2_unmasked_pvalue`, and `pydeseq2_unmasked_finite_family_padj`. They are diagnostics, not rescued primary inference.
- `count_supported_contrast`, `zero_required_count_arms`, all four pooled count totals, and each disease's number of zero-total donor pairs are retained. An individual zero-total donor pair is flagged, not deleted and not automatically equated with a wholly unsupported group response when other donors support both arms.

This is a necessary count-support rule, not a guarantee of general identifiability, convergence or small-sample calibration. The log-CPM sensitivity can remain defined on such rows because its pseudocount-based estimand differs; it does not rescue an unavailable finite NB log fold response. Four new synthetic tests cover both group/arm boundaries, preservation of the other group's contrast, individual zero pairs, raw diagnostics and unchanged correction families.

The design has no replicated identical full-design rows. Standard Cook eligibility at three identical rows would therefore cover zero libraries. **No automatic Cook exclusion is claimed or performed.** Save every library × feature Cook distance, the diagnostic F(10,6) 0.99 cutoff, fit coefficients, dispersion estimates and convergence flags. High Cook values and nonconvergence are limitations to report, not post-result reasons to delete donors or features.

The NB interaction is a **model-based log2 ratio-of-ratios**. It is not a raw-count difference, a cell-level effect, or the equal-donor mean difference of log-CPM changes used in the sensitivity. Four versus four donors and six residual degrees of freedom give **limited NB Wald calibration**. No synthetic test, software replay or many-cell count repairs that limitation.

## Donor-level sensitivity analyses

### Paired log-CPM changes and Welch uncertainty

Calculate `log2(CPM+1)` using totals over **all unfiltered released features** in each qualified input library. Only then subset the fixed eligible genes. Form one `IL17 − CNT` change per donor. Compare four PSC with four non-PSC changes using Welch's t test:

- Effect: `mean(change_PSC) − mean(change_nonPSC)`.
- Variance: `s_PSC²/4 + s_nonPSC²/4`.
- Satterthwaite df: `(v_PSC + v_nonPSC)² / (v_PSC²/3 + v_nonPSC²/3)`.
- Two-sided t-tail P and t-based 95% interval, not an asymptotic normal interval.

Finite paired observations are required; no silent case deletion, imputation or donor-specific feature filters are allowed. A zero or numerically degenerate variance produces an unavailable test, not a manufactured P=0. BH/BY sensitivity columns retain the fixed family but do not replace the primary NB rule. Every donor change remains in ignored work outputs for local review.

### Whole-donor leaveouts

Omit each donor's **entire pair**, then recompute the mean-change difference on the same features. Each omission leaves either three PSC/four non-PSC or four PSC/three non-PSC donors. Report the full effect, leaveout minimum/maximum and nonzero sign agreement/opposition counts. Preserve every omitted-donor effect in ignored work. These are influence diagnostics, not confidence intervals, resampling at the cell level or independent replications. There is no newly selected count filter or new NB fit per omission.

### Exact 70-label allocation diagnostic

The parent approved this sensitivity before implementation. Enumerate all `choose(8,4)=70` assignments of four PSC labels to the eight observed paired changes. Use the absolute mean-change difference. Preserve all 70 effects and donor-label assignments in ignored work, and the complete feature-level tail fractions publicly.

- Include the observed allocation and its complement; no Monte Carlo `+1` is used.
- An allocation is as or more extreme if `abs(allocated) >= max(0, abs(observed) − tolerance)`.
- The fixed tolerance is `1e-12 + 1e-12 * abs(observed)`.
- The minimum possible two-sided fraction is **2/70**, and an exactly zero effect has fraction 1.
- Label the result `descriptive_tail_fraction`, not a new discovery P or FDR criterion.

The diagnostic requires observational exchangeability of donor changes across disease labels. Clinical and age differences can violate this assumption. Disease labels were not randomized. Exhaustive enumeration is computationally exact conditional on that assumption; it does not make this a randomized-disease causal test. It has little resolution for sparse genome-wide discoveries and is retained only as a transparent small-panel calibration check.

## Covariates and interpretation

Do not add covariates without eight exact, independently qualified donor values and a separately justified frozen response model. Baseline donor attributes are collinear with donor fixed effects. Donor blocking does **not** remove age-by-treatment, sampling-by-treatment, clinical-history or culture-selection response confounding. A significant response within one group and a nonsignificant response within another does not demonstrate an interaction. Different DEG list lengths do not demonstrate it either.

Ex vivo treatment can support a culture response mechanism under the design assumptions. The between-disease response difference remains partly observational. It does not establish that IL-17 causes PSC, that a genetic variant acts through the response, or that an IL-17 inhibitor benefits patients. No benchmark/target/program score, transfer endpoint or tissue-validation criterion is added by this module.

## Source interface and execution contract

Script: `scripts/analyze_psc_organoid_il17.py`.

Input tables are tab-separated by default; an explicit `separator` pin can select comma instead. Gzip is accepted. Physical duplicate headers, ragged rows, missing IDs, duplicate IDs, mismatched sample/feature sets, invalid counts and zero-total libraries fail before fitting.

- `counts`: `source_feature_id` then exactly 16 `library_id` columns. Source-feature IDs remain unchanged.
- `features`: `source_feature_id` plus all original annotation columns. Order may differ; identity must join exactly. Unknown columns remain unchanged, including literal `NA` labels.
- `libraries`: `library_id`, `donor_id`, `disease` (`PSC`/`nonPSC`), `treatment` (`CNT`/`IL17`), `n_cells`, `total_released_counts`, plus preserved source annotations. Library totals must exactly match the complete pseudobulk matrix.

The JSON plan requires top-level `status="frozen"` and a past UTC `frozen_at_utc`. Its `organoid` object contains:

- `series="GSE239283"`, `population="all_published_retained_organoid_cells"`, the exact `donor_disease` dictionary, `expected_library_count=16`, `expected_cell_count=46343`, and the qualified `expected_feature_count`.
- `feature_namespace`, `feature_universe_description`, and optional `feature_id_column` (default `source_feature_id`).
- `counts`, `features`, `libraries`, `count_quantity_qualification`: artifact objects with `path`, SHA-256, and optional byte length/format fields.
- `source_artifacts`: original full-file objects with path, SHA-256, URL and UTC retrieval time. Partial/header-only files are not usable input pins.
- `code_artifacts`: SHA-pinned implementation/protocol objects; the executing script is mandatory. Pin tests, method/protocol and acquisition code too in the final execution plan.
- `filter`, `model`, `sensitivity`: the exact required entries exported as `FILTER_SETTINGS`, `MODEL_SETTINGS`, `SENSITIVITY_SETTINGS`. Additional provenance fields are allowed but never alter computations.
- `software`: exact Python, NumPy, pandas, SciPy, PyDESeq2, anndata, formulaic and joblib versions, plus any further lock versions. `current_software()` returns the required installed pins. A full dependency lock should also be a pinned artifact.
- `n_cpus`: explicit positive integer.
- `count_quantity="raw_RNA_counts"`, `raw_count_gate_passed=true`, `ready_for_inference=true` for inference. The quantity qualification must be reviewed scientific evidence, not merely an integer check.

```text
.venv/bin/python scripts/analyze_psc_organoid_il17.py \
  --plan config/psc-organoid-plan.json --plan-sha256 <frozen-plan-sha256> \
  --output-dir data/derived/psc-organoid-il17 \
  --work-dir work/psc-organoid-il17
```

`--validate-only` validates the pinned structure, donor pairing, totals, design, eligibility and ratio basis, but never computes log-CPM changes, NB effects or allocation results. It can run with a blocked quantity gate to document the problem, but still requires a frozen validation plan. It is not an override allowing SCT inference. **Both output and work directories must be new, including for validation and replay.** Existing directories are refused rather than reused, cleaned or partially overwritten. Every file is created exclusively, which also rejects existing symlink/hardlink destinations. These guards prevent stale completed effects, input mutation and private-output redirection. Source textual counts must match nonnegative integer tokens exactly before numerical conversion; fractional/scientific notation is refused rather than rounded through float64. No script reads a network source, installs a dependency or changes an input.

## Outputs and privacy

Every public gene table has the original feature identifiers/annotations and eligibility information. Excluded and unavailable rows are retained. Numerical NA is written as an empty field; source literal `NA` strings are unchanged. Read identifier columns as literal strings, not automatic missing tokens.

| Public filename | Coverage and main fields |
|---|---|
| `organoid-feature-universe.csv.gz` | Every source feature; `eligible`, count-threshold/nonzero library counts, all-zero flag, mean count. |
| `organoid-primary.csv.gz` | Every source feature; interaction `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `q_bh`, `q_by`, 95% Wald limits, count-support status, unmasked native diagnostics, fit/Cook flags, treatment–interaction covariance, explicit `status`. |
| `organoid-within-group.csv.gz` | Two rows per source feature; `contrast` is `nonPSC_treatment` or `PSC_treatment`; same NB effect/uncertainty fields. |
| `organoid-welch.csv.gz` | Every source feature; two mean changes, `effect_log2_cpm`, `se`, `t`, `df`, `pvalue`, sensitivity q values, 95% t limits, four/four donor Ns and status. |
| `organoid-leave-one-donor-out.csv.gz` | Every source feature; full mean-change difference, min/max omitted-donor effects, all eight sign diagnostics and status. No donor labels. |
| `organoid-allocations.csv.gz` | Every source feature; observed difference, 70 allocations, extreme count, descriptive fraction and status. |
| `organoid-qc.json` | Source/code/software pins; feature/donor/library/cell counts; design, normalization, fit flags and scientific limits. |
| `organoid-output-manifest.json` | Deterministic public-output hashes and byte sizes. |

Validation-only writes `organoid-feature-universe.csv.gz` and `organoid-validation.json`, not effects. The ignored work directory contains the frozen plan copy, selected libraries and design, normalization, all Cook distances, natural-log donor/model coefficients, all unrenamed native per-feature estimates/flags, fit warnings/log, every donor change, every leaveout effect, and all allocation labels/effects. Individual library/donor expression values never enter public outputs. Gzip mtime is zero for deterministic replay.

## Synthetic verification and readiness

Run `.venv/bin/python -m unittest discover -s tests -p test_psc_organoid_il17.py -v`.

Tests cover exact sample/feature joins; duplicate physical columns and source IDs; unpaired, mislabeled or missing donor fields; literal `NA` annotation preservation; counts/integer/library guards; rank/leverage and positive PSC interaction direction; unchanged designs under large cell-count changes; fixed count filtering; even-basis ratio math and prohibited fallback; unfiltered CPM totals; full-family BH/BY against statsmodels; Welch results/df/intervals against SciPy; missing/degenerate donor cases; whole-donor leaveouts; all 70 allocations, ties, complements and minimum resolution; full covariance in the PSC treatment contrast; synthetic NB response and retained flags; source/code/software/freeze gates; no-effect validation-only mode; complete public feature coverage, ignored donor outputs and exact synthetic replay.

All 38 organoid-method tests passed in the project environment, including two byte-identical isolated synthetic fits and the four prefreeze count-boundary tests. An independent code review found no arithmetic defect in the NB design/contrasts, correction, Welch, allocation or leaveout methods. It did identify three implementation guards that were repaired before any real effects: destination aliases, stale output-directory reuse, and float-rounded textual counts. Five added tests reproduce and reject those cases. The independent reviewer rechecked the repaired implementation and all 34 tests and found no substantive remaining issue within that code-review scope. The full repository run completed 434 tests with the same seven historical errors: four missing ETS2/GSE84161 cached-input errors and three missing GSE84161 R-runtime errors. The complete run is recorded locally in `work/psc-organoid-methods-full-suite.log`; it is not a clean full-suite pass. These checks confirm implemented arithmetic and contracts, not real source assay identity, small-sample NB calibration or biological independence. **Final source pins and the root execution freeze remain required.**

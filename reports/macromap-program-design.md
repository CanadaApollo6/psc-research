# MacroMap fixed-program context: prospective method proposal

Prepared 16 September 2026. **Source-only preparation is complete. This is not the execution freeze.** The parent owns the final protocol, inference and findings. No real MacroMap program score, paired effect, expression-matching reference, classifier fit or prediction has been computed by this preparation. Source quantities, identifiers and metadata were inspected. Existing source publications and the prior ETS2 results are already known.

This work is independent of the organoid result. A weak or opposing organoid result must not cancel it. It addresses the [first analysis in the frontier review](psc-analysis-frontier-review.md#1-macromap-the-strongest-immediately-executable-next-analysis), using the [completed local qualification](research-path-review.md#completed-bounded-macromap-qualification). No network request, dependency installation, outreach, enrollment or spending was needed.

## 1. Question and interpretation

Does the fixed ETS2 gRNA1-down program, excluding ETS2, carry held-out **stimulus-context information** beyond the four previously fixed broad response programs? Finish two distinct outputs:

1. The entire paired response panel for the unchanged nine programs, all nine qualified stimuli and both times, including failed/unavailable rows.
2. One primary predictive endpoint: reduction in held-out, line-weighted multiclass log loss from adding ETS2 to the four-comparator baseline. Six-hour and 24-hour models remain separate. They enter the final endpoint once, with equal weight.

MacroMap is iPSC-derived macrophage stimulation, not new PSC subjects, a clinical diagnostic study, direct ETS2 perturbation, or proof of ETS2 mediation. Positive incremental information would support this particular context score relative to these particular comparators. Failure to improve, protocol disagreement or loss of the signal after overlap removal would limit utility. An interval crossing zero is imprecision, not proof of exact redundancy. No clinical benefit follows from either outcome.

The original study selected successfully stimulated samples using expression QC and published stimulus-response biology. Responsiveness itself is not a new finding. The paper describes unrelated starting HipSci donors; this is not a new kinship verification. The bounded Stankey-source check found no MacroMap reuse, not universal donor/source disjointness.

## 2. Source reuse and completed gates

Use the existing caches under ignored `data/raw/macromap-qualification/`. The new script hard-pins the raw file, annotation, Supplementary Data 11, source receipts, metadata, header, feature-universe audit and the historical program inputs. Retrieval was on 15 September 2026. Important unchanged sources:

| Input | Source / identity | SHA-256 |
|---|---|---|
| Released raw counts | [Zenodo 11563707](https://zenodo.org/api/records/11563707/files/MacroMap_raw_expression.txt.gz/content), 174,122,135 bytes | `eb9dd32b6f25ae9912e47dcc5e1f5170a98af919a42d48a4bca3a680fc5353d9` |
| Original-study GENCODE v27 | [Pinned tree `191fc67f84d631e5b4fc9a97993c142b348cf2c9`](https://raw.githubusercontent.com/andersonlab/macromap_eqtl/191fc67f84d631e5b4fc9a97993c142b348cf2c9/Data/gencode.v27.annotation_chr_pos_strand_geneid_gene_name.gtf) | `db11e855d959656bc3a25de65460ae12368ca13b1ae1a5493973f6ea433bdf41` |
| Supplementary Data 11 | [Original RNA-seq metadata worksheet](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-61670-9/MediaObjects/41467_2025_61670_MOESM13_ESM.xlsx) | `97805429aaf777d8d13f625375a0dee124f130a2499fd6237bec279295bd89d8` |
| Expression paper | [Panousis et al., DOI 10.1038/s41467-025-61670-9](https://doi.org/10.1038/s41467-025-61670-9), PMC12391345 | XML `d3f69f885db807ff964f9158c7b7f581bea649225d6c76a37e4d3ef34f38cb8f` |
| Original program plan | `config/ets2-program-plan.json` | `234560f51bc68af2b9f3258af222894e68c9e992c89ba7d2eae396673fb8c9c9` |
| Effective validation amendment | `config/ets2-validation-amended-rules.json` | `cbbaa6e626dad9f494367cf1a8c45c0db270f4c330f9bcd61a72280fe1b8ef6b` |
| Complete original source-value export | `data/derived/ets2-program-membership.csv` | `b6e468f000802e6f1ec6baa53722e8a7d3e0f98d6df4aa2bb49ef886efc199a3` |

The expression methods declare GRCh38, GENCODE v27, STAR 2.5.3a and featureCounts 1.5.3. The earlier companion splicing-paper link remains a distinct source. Neither an allele direction nor a genetic association is estimated here.

### Feature universe and the historical-gate boundary

The previous full-file audit verified 58,243 unique gene rows × 4,698 samples, all 273,625,614 counts finite, nonnegative and numerically integral, 2,397 all-zero genes, positive library totals and gzip EOF/CRC. The preparation CLI rechecks byte/hash pins and the actual gzip header, not count-valued rows. It independently joins all released versioned IDs to the pinned study annotation.

After the exact-token guard repair, a separately authorized local source-only scan confirmed that all **273,625,614 original count tokens are exactly finite, nonnegative integers below 2^53**, and all full released-universe library totals are positive and below that limit. Source hashes agree before/after the scan and gzip EOF/CRC passes. Only whole-file validity flags and dimensions are saved in ignored `work/macromap-design/exact-count-token-audit.json`; no normalized values, target-gene result, matching reference, program score or prediction was produced. The source bytes were not changed.

The reference has 45 additional `_PAR_Y` records. They remain absent, not imputed or collapsed. Count totals will use **all 58,243 released features**, including zeros and features outside every signature. No signature, abundance, protein-coding or differential-expression filter is allowed before totals or reference binning.

This is a **new released-universe normalization estimand**. It does not claim recovery of every historically counted feature or reconstruction of the author's filtered 14,060-gene workflow. The historical full-universe eligibility rules in the earlier validation protocol remain unchanged. The new MacroMap-specific protocol must expressly authorize this qualified released-universe quantity, rather than retroactively claiming that an unobserved historical denominator has been recovered.

### Recovering the unchanged source memberships

The original `data/raw/ets2-benchmark/ets2-public-release.zip` is absent locally. No replacement archive was downloaded. The complete, hash-pinned historical membership export retains each original `source_value`, including entries that did not map in the liver atlas and ETS2 exclusions. The two independent modality copies agree for all eight original lists. The four Ensembl-defined lists also agree with the complete pinned `ets2-source-gene-sets.csv`. The old audit pins the membership export and the original source members. This is reuse of a verified export, **not a new rereading of the absent ZIP**.

Use the entire source lists from this export, not their previously mapped feature subset. Match Ensembl IDs after removal of terminal numeric versions only. Otherwise require a unique exact symbol in GENCODE v27. Never use aliases, smallest P values or expression to choose mappings. Stable-ID collisions stop; symbol ambiguity makes that member unavailable, while both source count rows remain in the denominator. Exclude `ENSG00000157557` / `ETS2` from every program.

| Program | Intended after ETS2 exclusion | Mapped | Gate |
|---|---:|---:|---|
| ETS2 gRNA1-down, primary | 927 | 926 | pass |
| ETS2 gRNA1-up | 668 | 668 | pass |
| ETS2 gRNA2-down | 875 | 874 | pass |
| Chr21 enhancer-deletion-down | 141 | 141 | pass |
| Inflammation | 815 | 807 | pass |
| Interferon-γ response | 139 | 136 | pass |
| Oxidative stress | 436 | 431 | pass |
| Apoptotic signaling | 588 | 580 | pass |
| Primary minus mapped comparator union | parent 927 | 759 | parent gate passes; ≥10 remain |

The original ≥80% and ≥10-member gates are unchanged. The overlap sensitivity removes 167 mapped primary members. Its gate uses the parent coverage, not a newly chosen denominator. Related guide/enhancer programs are companions, not independent replications. They are not additional predictive features.

### Exact sample and pairing interfaces

The new standard-library XLSX reader reconstructs Supplementary Data 11's `RNA-seq` sheet and agrees with all 4,698 rows and 21 fields of the earlier CSV. Every matrix column equals `SampleID + "_" + RunID`; every source line is obtained by removing the exact `"_" + Stimulus_Hours` suffix from `SampleID`. The independent earlier sample join agrees. No anonymous suffix is treated as a new donor identity.

All static metadata are constant within a source line, including library protocol, sequencing run, HipSci label, sex and recorded culture/differentiation fields. One sample occurs per line/condition. Multiple cultures, technical splits, or conflicting mappings would stop the current procedure rather than trigger a new collapse rule. One line maps to only one run/protocol; every run is protocol-specific; mapped HipSci labels are one-to-one with source lines.

There are 205 mapped lines: 124 `Mod_SmartSeq2` and 81 `NEB`. Four separate `NOMATCH` lines are excluded from the primary analysis and never collapsed into one person. All 209 lines span 57 runs: 34 and 23 by protocol.

The strict 18-condition response panel uses every available within-line treatment / same-time basal `Ctrl` pair. The full 20-condition registry keeps PIC at both times explicitly `unavailable_mock_identity`; basal control does not isolate lipofection from PIC. Precursors are not stimulation controls. CIL and LIL10 are mixtures; MBP is myelin basic protein.

| Stimulus | Mod_SmartSeq2 6h / 24h | NEB 6h / 24h | Total 6h / 24h |
|---|---:|---:|---:|
| CIL | 100 / 98 | 69 / 72 | 169 / 170 |
| IFNB | 115 / 115 | 67 / 77 | 182 / 192 |
| IFNG | 112 / 115 | 70 / 76 | 182 / 191 |
| IL4 | 117 / 114 | 72 / 76 | 189 / 190 |
| LIL10 | 114 / 111 | 70 / 71 | 184 / 182 |
| MBP | 100 / 106 | 71 / 64 | 171 / 170 |
| P3C | 114 / 117 | 68 / 76 | 182 / 193 |
| R848 | 112 / 116 | 70 / 75 | 182 / 191 |
| sLPS | 114 / 111 | 65 / 76 | 179 / 187 |

These are source-design counts, not measured effects or independent condition-level replications.

## 3. Proposed fixed split and scoring procedure

### Whole-run five-fold holdout

Use five folds **within each protocol**, shared across times and every condition of a line. Include the complete source run inventory in the assignment so the identity sensitivity never reassigns primary runs. For each protocol, sort runs by `(sha256(f"20260916:{protocol}:{run}").hexdigest(), run)` and assign sorted ordinal modulo five, numbered 0–4. Do not balance folds using expression, stimulus effects or predictive performance. `RunID` remains its original string. No alternate split seed is tried.

The realized run assignment, joins and pairs remain ignored; the final freeze must hash them. Every line's test response is computed only in its assigned fold. Training and test lines/runs are disjoint. The source-only checks found all 20 protocol × time × fold design gates satisfied. No test-fold stimulus is absent. For the predictive cohort, test folds contain 20–26 lines / 6–7 runs in Mod_SmartSeq2, and 12–19 lines / 4–5 runs in NEB. Minimum training class coverage is 72 lines / 23 runs and 45 lines / 16 runs, respectively.

### Normalization, control matching and paired responses

For sample s and released gene g, let

`L[g,s] = log2(1 + 10^6 * count[g,s] / sum_all_released_genes count[g,s])`.

This is per-sample normalization. No normalization factor is estimated from another sample, and no additional outcome-based QC exclusion is allowed. It is not TPM, TMM, voom or an attempt to remove every composition effect.

For each protocol × time × held-out-fold stratum:

1. Select **training-run matched basal controls only**, from mapped lines with at least one eligible strict pair at that time. Each shared control enters once. The descriptive response reference may include training lines lacking the other time; predictive model fitting uses the both-time cohort below. Require ≥3 reference lines from ≥3 training runs.
2. Average each gene's transformed control abundance with equal line weight. No treatment or held-out control determines bins, matched sets or their QC summaries.
3. Sort all released genes by training-control abundance, then stable feature ID, then original versioned ID. Apply `numpy.array_split` into 20 equal-count bins. Measured zeros remain eligible. Exclude the union of all eight original mapped programs plus ETS2 from the background pool, including programs that would fail a mapping gate.
4. Preserve 1,000 matched sets, without replacement within a set and ascending bins per draw. Seed PCG64 from the first 16 SHA-256 hexadecimal characters of `20260912:macromap:protocol={protocol}:time={time}:fold={fold}:{program_id}`. Candidate order is the frozen abundance/ID order. Freeze program names exactly as in the source plan. Record selection hashes and matching errors. No bin widening or fallback on pool failure.
5. Score every training/test sample needed in this stratum as mean mapped-program `L` minus mean of its 1,000 matched-set means. The same stratum/program draws apply to all lines and all nine stimuli. Keep unadjusted mean, background mean, matched score and detection counts separately. No gene is removed because it is zero.
6. For each eligible line/stimulus/time, take treatment score minus its own same-time control score. A line's reported response always uses its held-out-fold scoring reference. Do not later replace it with an all-controls fit selected for a cleaner result.

The score formula, bin count, draw count, no-replacement sampling, seed family and mapping gates reproduce the earlier ETS2 conventions. **The new estimand is a run-cross-fitted, within-line stimulus response in bulk iPSC macrophages**, not the old atlas cell-type contrast. Fold/time/protocol-specific control references are a new prospective transfer procedure. Negative adjusted scores or differences do not by themselves demonstrate biological inhibition. Matched gene sets are abundance references, not independent donors or disease-null samples. Optional conditional gene-set percentiles would not be population P values and are not a new endpoint here.

The reference mean over 1,000 sets can be computed exactly, up to floating-point summation, by counting each gene's selection frequency and using those linear weights. Synthetic tests compare this with explicit random-set means. Keep the realized selection hash; do not replace the random-set average by a different deterministic pool mean.

## 4. Complete response panel

For each protocol, program and all 20 registry conditions, retain paired-line count, contributing-run count, mapping/matching status, median, mean, full range and positive/negative/zero counts. PIC has no effect estimate. Preserve unadjusted and matched-background changes as separate summaries. Report detection and matching discrepancies without using them to choose arms.

Require ≥3 paired lines for the original descriptive summaries. Retain leave-one-line-out median ranges when ≥3 pairs remain, and leave-one-run-out ranges only when ≥3 pairs remain. Neither range is a confidence interval. Companion program or matching failures do not remove other programs from the full table.

Report protocols first. A pooled descriptive row may then recompute the median/mean over all eligible lines for that exact condition; do not average two protocol medians or count protocol results as separate donor replications. Different donor pools prevent interpreting protocol differences as an identified library-preparation effect. Separate times remain separate response-panel endpoints; only the predictive endpoint below is combined over time.

Use pointwise run-clustered bootstrap intervals on the frozen cross-fitted line responses: 5,000 PCG64 draws, seed `2026091601`, lexical run order, sampling runs with replacement within protocol and keeping every line in a selected run together. Require ≥3 observed runs for an interval. **Final prefreeze clarification, 16 September 2026:** both pooled response **mean and median are descriptive only**, with `mean_interval_status=not_computed_prespecified` and `median_interval_status=not_computed_prespecified`. No pooled-response confidence interval or weighted-quantile convention is computed. This supersedes the earlier prefreeze request for a contrast-weighted pooled-mean interval; no real effects had been examined. Protocol-specific pointwise mean and median intervals remain. The separate primary **pooled predictive-gain interval remains mandatory**, using its fixed 114/185 and 71/185 line shares and equal time weighting. These are distinct outputs and uncertainty rules. Full-panel simultaneous inference is not claimed. No P-value or “significant stimulus” screen is planned. The script supplies the single-protocol response helper; the parent must preserve this stratification for any pooled bootstrap implementation.

These intervals condition on learned references and fixed fold assignments. They do not quantify the full variation of retraining the scoring procedure. Fold-specific descriptive summaries and leave-run-out ranges expose gross split/run dependence without pretending to resolve that limitation.

## 5. One primary held-out incremental prediction endpoint

### Predictive target population and missing measurements

Use lines with ≥1 strict paired condition **at each time**. This metadata-only rule yields **185 mapped lines: 114 Mod_SmartSeq2 and 71 NEB**. It excludes 20 lines lacking a paired time, not weak program responses. The response panel still uses every contrast-specific pair. The predictive cohort spans 34 eligible Mod_SmartSeq2 runs and 22 eligible NEB runs; the original assignment includes all 57 source runs.

Within line and time, use all available strict stimuli and weight them equally. Do not impute an absent stimulus. Each time's nine-class risk is conditional on this released, available-condition pattern; it is not performance on a complete unobserved factorial experiment. Keep its frozen missingness table. At least one evaluation observation from each of nine classes is required in each protocol/time cell. Test-fold class absence would be listed without moving the fold; a whole-cell class absence makes the primary endpoint unavailable.

### Exactly one model class, no tuning

Fit separate models for each protocol × time × outer fold. There is no pooled-protocol fit, feature selection, nonlinear model, stimulus search or regularization search.

- **Baseline columns:** `inflammation`, `interferon_gamma`, `oxidative_stress`, `apoptotic_signaling` paired changes, in that order.
- **Extended columns:** the same four plus `ets2_g1_dn` paired change.
- **Class order:** `CIL, IFNB, IFNG, IL4, LIL10, MBP, P3C, R848, sLPS`.
- **Training weight:** for line i with k_i available stimuli at the fitted time, each row weight is `1 / (N_train_lines * k_i)`. Weights sum to one, and every line has equal total weight. No class-frequency balancing is added.
- **Predictor standardization:** learn weighted column means and population standard deviations from the **training-fold paired changes only**, using those weights. Test values never enter this step. SD ≤1e−12 uses scale one and is flagged; a constant training feature has zero fitted slope. The baseline uses the same first-four-column transforms as the extended fit. This standardization is distinct from control-only gene-reference matching.

For classes c, predictors z, intercepts a and slopes B, minimize

`sum_rows w_r * [-log softmax(a + z_r B)[y_r]] + (0.1 / 2) * sum_features,classes B[j,c]^2`.

This is a **weighted mean**, not summed unnormalized row loss. All class slopes are penalized symmetrically. Intercepts are unpenalized. Use all nine softmax columns, zero initialization, and center each parameter row over classes to remove the common-logit shift. Both models use the fixed λ=0.1. No penalty change is allowed after an unsuccessful fit.

Use SciPy L-BFGS-B with analytic gradient, `maxiter=2000`, `ftol=1e-13`, `gtol=1e-8`, `maxls=50`. Require optimizer success, finite objective/coefficients/gradient and maximum absolute final gradient ≤1e−6. Evaluate natural-log probabilities via log-sum-exp, not clipping zero probabilities. Exact prediction ties use the first frozen class label. The primary metric is log loss, not accuracy; optional confusion tables and accuracy are descriptive only.

Before each model, require ≥3 distinct training lines and ≥3 distinct training runs **for every class**. If a training class, any required mapping/matching feature, or either paired model fit is unavailable, retain the explicit failure. Do not drop difficult rows/classes/folds, use a uniform-probability replacement, shrink the class set, or refit a different model. Any missing expected prediction makes the primary endpoint unavailable. Successful parts can be reported as incomplete diagnostics, not substituted into the complete primary estimate.

### Exact aggregation, once across time

For model m and line i at time t, define

`loss[m,i,t] = mean_available_stimuli(-log P_out_of_fold[m](true stimulus))`.

Let `gain[i,t] = loss[baseline,i,t] - loss[extended,i,t]`; positive means ETS2 improved held-out log loss. First report all four protocol × time mean losses/gains, and each fold's counts/status/losses. Then define each line's two-time gain once:

`gain[i] = (gain[i,6] + gain[i,24]) / 2`.

The single primary endpoint is

`D = sum_185_lines gain[i] / 185`.

Equivalently, average times equally within each protocol and pool the two protocol estimates with **fixed line shares 114/185 and 71/185**. Equal protocol weights would be a different estimand. Do not pool treatment rows directly, average folds equally despite unequal line counts, or count the two times as independent replications. Report the 6h/24h cells before D.

### Conditional run-cluster uncertainty

Use 5,000 PCG64 bootstrap draws with seed `2026091601`. Within each protocol, resample its eligible run IDs with replacement, in lexical run order. Each draw moves all conditions, times and lines of a run together. Recompute line-weighted protocol/time means with run multiplicities, then use the frozen 114/185 and 71/185 protocol weights and the same once-only 1/2 time weights. The 2.5th and 97.5th percentiles use `numpy.quantile(method="linear")`. Report primary and four-cell intervals, contributing runs, and fold/run diagnostics. Do not count 18 repeated conditions as 18 independent outcomes per line.

This bootstrap **holds fitted models, matched references, scale estimates and fold assignments fixed**. Its interval is conditional uncertainty in the observed held-out loss contrast, not full procedure/retraining uncertainty or an unconditional classifier-generalization confidence interval. Overlapping training sets create dependence not removed by treating folds as replicates. No formal unconditional superiority P value is proposed. A narrow conditional interval alone cannot establish a robust transferable mechanism.

## 6. Fixed sensitivities and stopping rule

1. Replace only the ETS2 predictor/program with `ets2_g1_dn_without_comparators`, using its unchanged parent mapping gate and separately seeded matching sets. Keep the same baseline, folds, eligible lines, times, λ, classes and metrics. It is a sensitivity, not a second primary endpoint or a choice of whichever performs best.
2. Keep a separately labeled identity sensitivity that admits all four individually named `NOMATCH` source lines when they have the required pairs. Reuse the original run-fold assignment, rebuild training-only references/scales using that sensitivity's training lines, and use its own metadata-fixed line shares. Its identities are less qualified; it cannot override the mapped-line primary result.
3. Retain all folds, protocol/time cells and leave-run-out influence ranges even when they disagree. Do not add signatures, drop unfavorable stimuli, select time points, search aliases, or change λ to rescue a weak result.

Stop after the complete frozen response panel, primary endpoint, these sensitivities and checks. Matching-pool failure or classifier failure is a result status, not permission to start another search.

## 7. Implementation, resource plan and data boundaries

Owned files:

- `scripts/analyze_macromap_program_context.py`: source-only CLI, fixed method helpers, and the complete post-freeze engine. Its default CLI never calls the numerical engine.
- `scripts/run_macromap_program_context.py`: separate parent-freeze-gated executable entrypoint; validation-only unless `--execute` is supplied.
- `tests/test_macromap_program_context.py`: synthetic method/source/precision/path/concurrency tests.
- `tests/test_run_macromap_program_context.py`: synthetic end-to-end engine, freeze/changed-source guards, audit reconstruction, failed-matching preservation, and the final no-pooled-response-interval rule.
- `reports/macromap-program-design.md`: this method/readiness proposal.
- Ignored `work/macromap-design/`: realized sample joins, pairing/folds, mapping details, metadata-only predictive cohort, source-only aggregate readiness and test logs. Real counts, scores and predictions must stay ignored if later produced.

Run source-only preparation with:

```text
.venv/bin/python scripts/analyze_macromap_program_context.py --output work/macromap-design/preparation-runner-v3
.venv/bin/python -m unittest discover -s tests -p test_macromap_program_context.py -v
.venv/bin/python -m unittest discover -s tests -p test_run_macromap_program_context.py -v
```

The CLI verifies source pins, reads only the raw count header, compares the original worksheet and cached joins, remaps full unchanged source lists, constructs the 20-condition registry and fixed folds, then writes ignored preparation artifacts. It emits `ready_for_inference=false`. It cannot launch an actual scoring/model run. Real execution is the parent's separate post-freeze step.

`load_prepared_design(directory, expected_summary_sha256=...)` is the intended numerical-run handoff. The parent supplies its frozen readiness-summary hash; the loader verifies all artifact hashes and restores integer times/indices and boolean identity flags from CSV. It returns `samples`, `pairs`, `predictive_pairs`, tuple-keyed `folds`, `mapping_gates`, `program_indices`, original `feature_ids`/`sample_ids`, and the eight-program-plus-ETS2 `pool_exclusion_indices`. This avoids interpreting the literal CSV string `False` as truthy. Loading a design does not authorize inference.

Every preparation or execution uses a **fresh run directory** below the literal `work/macromap-design/` anchor. Reusing a run directory is an error, not an overwrite operation. Source/destination symlinks, including dangling links and aliased path components, are rejected. Literal `..` components are rejected before path normalization, so a symlink followed by parent traversal cannot bypass a helper's guard. Existing preparation snapshots and source files remain unchanged. Numerical caches reserve an exclusive namespace through cache and sidecar publication; memmaps use an exclusively created file handle, and finished files are published without replacing an existing name. A successful `.lock` marker prevents namespace reuse. Synthetic concurrent writers cannot share or keep changing a completed cache inode.

### Separate parent execution freeze and executable runner

`execution_plan_template(prepared_directory, output_directory)` returns a draft dictionary with `execute_real_expression=false`. The parent must review it, set a timezone-aware freeze timestamp, explicitly authorize real expression, and externally pin the final JSON bytes. The final plan binds both method and runner code hashes, Python/NumPy/SciPy versions, the prepared-summary hash and the entire finite method contract. It must use a new ignored output directory. No readiness timestamp or prior preparation hash grants permission by itself.

The separate entrypoint validates only unless `--execute` is explicit:

```text
.venv/bin/python scripts/run_macromap_program_context.py --plan <parent-frozen-plan.json> --expected-plan-sha256 <sha256>
.venv/bin/python scripts/run_macromap_program_context.py --plan <parent-frozen-plan.json> --expected-plan-sha256 <sha256> --execute
```

The real-data path checks all source/preparation/code/software pins, creates one full released-universe normalized cache, and compares that cache's source hash with the frozen raw-source pin **before any matching or fit**. Raw-source identity/size/modification metadata must remain unchanged during normalization. It then executes both fixed identity scopes, all fold/protocol/time references, every response program and both fixed predictive variants. Both training and held-out predictor matrices use the **same current outer fold's training-control reference**; a global table of own-fold out-of-fold training responses is never substituted. Before a completion record, the plan, source, preparation and code/software pins are rechecked, along with the cache hash. An error writes a private failed-execution receipt, not a complete primary result.

All expected fit identities remain in `fit-manifest.json`, including missing-reference, matching, class or fit failures. Missing expected held-out keys make the primary endpoint unavailable. The runner never silently drops failed folds, substitutes scores/probabilities, changes λ or searches a model. Only aggregate tables may later be promoted by the parent; the runner itself publishes everything below ignored `work/`.

### Manageable CPU/memory strategy after freeze

- Stream the existing concatenated gzip with the standard gzip reader, not a first-member-only decoder. Check exact source feature/sample order and exact lexical count validity before float conversion, then accumulate full released-universe library totals in float64. A nonnegative integer/zero-fraction syntax fast path handles ordinary count exports; other numerical syntax uses exact `Decimal` validation. Fractional near-integers and negative/positive underflow tokens must fail rather than round to integers or zero.
- Write one gene × sample float64 memmap, **2,189,004,912 bytes** (~2.04 GiB). Float64 preserves integral source values below 2^53; reject larger values/totals rather than lose integer precision. The raw gzip remains the original count source.
- In 256-row blocks, transform the memmap to log2(CPM+1). Exclusive, no-replacement publication and a completed sidecar identify finished normalized data. Do not hold a full pandas matrix or several full transformed copies. Tests verify chunk invariance, alias rejection, concurrent namespace isolation and failed-cache cleanup.
- Read each fold's training controls in bounded blocks if necessary; a complete reference vector is only 58,243 doubles. Original sample IDs and original versioned feature IDs are retained. A reference control slice is roughly 50 MB, not a whole experiment copy.
- Convert sampled reference sets into linear frequency weights. Nine program/background columns require less than 10 MB per stratum; score sample blocks or multiply a gene block into a compact sample × program result. Avoid an observations × 1,000 × genes tensor. Preserve hashes of selections and matching QC.
- Each identity scope has 20 fold/time/protocol cells and two predictor variants. Each variant audit is self-contained and includes its baseline and extended fit; the unchanged baseline is thus recomputed identically for the nonoverlap audit. Across both identity scopes this is 160 very small fits, not 160 independent replications or a model search. Bootstrap fixed predictions from run-level sufficient sums/counts rather than refitting thousands of models or copying count matrices.
- Allow roughly 4–8 CPU cores, 8–16 GB RAM and <5 GB of new cache/temporary disk for the main run; run core model code with a pinned BLAS-thread setting to avoid oversubscription. A 0.5–3-hour envelope is conservative planning, not a measured real-score runtime. No GPU, eQTL archive or raw-read reprocessing is needed.

### Agreed private independent-verification bundle

The bundle schema is `macromap-private-audit-v1`:

- Per stratum, `score-reference-arrays.npz` contains the full training-control reference abundance vector, gene-bin assignments, excluded pool indices, flattened target indices and program offsets, all nine target/background weight vectors, matching availability, training-control indices/IDs, selected sample indices/IDs, raw/background/adjusted sample scores, detection fractions and feature/sample axis hashes. Unavailable background weights/scores are NaN, not zero. `score-reference-audit.json` records protocol, time, outer fold, identity scope, program statuses, reference counts, seeds, selection hashes, bin pool/target counts and baseline-matching discrepancies.
- Per fit, `fit-arrays.npz` contains `train_x/test_x` (unscaled paired changes, n×5), standardized matrices, `train_y/test_y` (0–8), `train_lines/test_lines`, `train_runs/test_runs`, sample IDs/indices, `train_pair_ids/test_pair_ids`, `predictor_names`, `scaler_mean/scaler_scale`, constant-feature flags, line weights, `baseline_coefficients` (5×9 including intercept), `extended_coefficients` (6×9), `baseline_logp/extended_logp`, and per-row losses/gain. Failed fits omit unavailable coefficient/probability arrays; they are not filled with fake predictions.
- `fit-manifest.json` retains all 80 expected scope/variant/protocol/time/fold identities and statuses and links their NPZ and stratum audits. Each `fit-audit.json` contains optimizer convergence, objective/gradient, class order, λ and the training-only preprocessing rule. Fold diagnostics distinguish design-gate status from actual fit status.
- Each scope retains the complete 20-condition pairing registry, expected predictive keys, line responses, response panel, loss rows, endpoint results and leave-run-out influence diagnostics. The root execution record includes frozen plan/preparation/source/code/software pins and the normalized cache's hash/shape/dtype. A final file-hash manifest pins the entire private bundle.

For bounded independent source-to-score verification, fix outer fold 0 for both protocols and both times. Select four training and four held-out predictive pair keys by SHA-256 of `macromap-independent-v1:` plus the canonical key; check all nine programs. A canonical pair key is compact JSON `[protocol,run,line,integer_time,stimulus]`, with no separator spaces. Separately select eight of all 4,698 source sample IDs using the same prefix/hash ordering and reconstruct their full 58,243-gene CPM denominators from the raw file. These are verification subsets, not new model comparisons or endpoints. The independent reviewer owns those computations and any required missing-output status.

### Suggested final schemas

Public aggregate outputs should preserve explicit status fields:

- **Mapping:** `program_id, intended_genes, mapped_genes, missing_genes, ambiguous_genes, parent_coverage, overlap_removed, gate_status`.
- **Design:** `identity_scope, protocol, time, stimulus, available_pairs, runs, status`; plus protocol/fold training and test counts, missing-class and fit statuses. No line IDs.
- **Response:** `identity_scope, protocol_or_pool, time, stimulus, program_id, quantity, n_lines, n_runs, median, mean, min, max, positive, negative, zero, interval_status, interval_low, interval_high, leave_line_out_range, leave_run_out_range, status`.
- **Prediction:** `identity_scope, program_variant, protocol_or_pool, time_or_equal_time, fold_or_all, n_lines, n_runs, n_pairs, baseline_log_loss, extended_log_loss, gain, interval_status, interval_low, interval_high, status`.
- **Fit audit:** `protocol, time, fold, model, lambda, optimizer_status, iterations, gradient_inf, objective, reference_line_count, reference_run_count, matching_status, scaling_status`, with source/code/runtime hashes.

Ignored tables should keep the source-preserving `sample-map.csv`, `pair-registry.csv`, `predictive-pairs.csv`, `run-folds.csv`, complete program-member audit, selected reference/control IDs, bin/selection hashes, normalized cache, line response rows and all nine log probabilities for every expected held-out row. Aggregate routines compare the full expected prediction keys to realized keys; a missing output is not silently dropped.

## 8. Readiness and remaining boundary

**Completed:** source pins; original worksheet replay; sample/run/protocol and time-pair joins; full unchanged source-list recovery; all mapping gates; metadata-only prediction cohort; deterministic folds and class coverage; synthetic score/multinomial/cluster checks.

Prefreeze review found and repaired exact-token, output-alias and concurrent-cache guard defects using synthetic probes. No real source was corrupted. The fixes reject rounded fractional/underflow counts, source/destination symlinks, literal parent traversal, reused run directories and aliased artifacts. Cache namespace reservation remains exclusive through cache and sidecar publication. Raw-source changes during normalization and source/code/plan/cache changes before completion stop execution rather than publish a result. All repairs precede any real normalization, matching or fit.

All **71 current native MacroMap tests pass**: 40 method/guard tests, eight executable-runner tests and 23 independently maintained verification tests. Their end-to-end inputs are synthetic. Both successful and matching-failed 80-fit manifests are tested, with explicit status/expected-key retention and independently reconstructed scaler/log-probability records. Independent guard replay passes **100 assertions** against the final method and runner hashes. The independent source-only audit also reproduces 205 mapped lines, 185 predictive lines, 57 runs, all 20 fold gates, nine programs and 167 overlap removals. Logs are in ignored `work/macromap-independent-review/final-targeted-tests.log` and `source-only-final-check.json`.

Full-repository logs remain under ignored `work/macromap-design/`; seven known errors concern absent legacy ignored caches/runtime rather than MacroMap. A transient independent-test wrapper still expecting pooled-bootstrap arguments was repaired by its owner when the final descriptive-only rule removed those arguments. The current 71-test MacroMap pass does not claim an all-repository pass.

Fresh `preparation-runner-v3/` and `preparation-runner-replay-v3/` reproduce all ten source-only artifacts byte for byte. All nine realized gene/sample/pair/fold design artifacts match the preceding preparation exactly; only preparation-code provenance changes. Prior preparation snapshots remain untouched. Native preparation uses Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. The draft `execution-plan-template-v3.json` remains explicitly unauthorized; the parent must freeze and pin it separately after review.

**Not yet observed:** real matching pool sufficiency; real expression/reference/scaling values; optimizer convergence on real inputs; response effects; predictive losses; clinical outcomes. These are assessed only after the parent's distinct execution freeze. Record failures without changing the method.

The absent original ETS2 ZIP is an explicit provenance boundary, not a missing membership input under the approved complete-export route. No author response or further source search is required for this proposed analysis. The remaining action is parent review and freezing of the exact protocol, realized metadata artifacts, input/code hashes, runtime versions and status rules. Preserve all earlier scientific files and results.

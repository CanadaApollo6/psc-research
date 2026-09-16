# MacroMap results: bounded independent interpretation review

Reviewed 16 September 2026. This is a **new post-execution review**, not an edit to the frozen [prefreeze review](macromap-independent-review.md), protocol or code. It reads the completed verifier receipt and primary aggregate outputs. It does not rerun verification, scoring, models or bootstrap calculations. No new model was fitted to explain the result.

## 1. Evidence and disposition

The inspected aggregates match the primary completion manifest. The completed frozen verifier reports `independently_verified_with_declared_source_subset`.

| Record | SHA-256 |
|---|---|
| Root execution freeze | `7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df` |
| Primary `execution-complete.json` | `e85445627230c7bdae025df48c80e01c34952da0c476237bfe099f351cee1408` |
| `postfreeze-verification-v1.json` | `3cfaf66e53481883f04916935e7c1945fe664caf32b5cf29f962fb192d665e96` |

Primary outputs are under ignored `work/macromap-design/primary-execution-v1/`. The verifier receipt is under `work/macromap-independent-review/`. The root reports execution after pre-effect commit `88c7464`.

**Disposition:** the inspected outputs support a positive held-out log-loss improvement for the prespecified fixed model comparison. They do **not** establish uniquely new conditional ETS2 information, ETS2-specific activity or causality, PSC validation, or clinical utility. No technical discrepancy was found in these inspected outputs. The audit limits below remain material.

## 2. Primary result, without substituting a sensitivity

The primary population is 185 mapped lines: 114 `Mod_SmartSeq2` and 71 `NEB`. Losses average available stimuli within line/time, then average the two times once and use the original protocol line shares. The unit below is natural-log loss, not classification accuracy or a clinical effect size.

| Primary cell | Baseline loss | Extended loss | Gain | Conditional 95% interval for gain |
|---|---:|---:|---:|---|
| Mod_SmartSeq2, 6h | 1.552563 | 1.402409 | 0.150154 | [0.138836, 0.161265] |
| Mod_SmartSeq2, 24h | 1.653192 | 1.556822 | 0.096370 | [0.075380, 0.110284] |
| NEB, 6h | 1.682576 | 1.503107 | 0.179469 | [0.161579, 0.194691] |
| NEB, 24h | 1.824523 | 1.717555 | 0.106968 | [0.094502, 0.122046] |
| **Primary pooled, equal-time** | **1.660703** | **1.529782** | **0.130921** | **[0.124191, 0.136547]** |

All four primary cell gains and all 20 primary held-out fold gains are positive. Fold gains range from 0.069599 to 0.200020. The 56 fixed-prediction leave-run-out evaluation summaries range from 0.129730 to 0.133208. These are descriptive consistency checks, not independent replications. Leave-run-out evaluation does not remove that run from the already learned references or training fits.

The intervals use 5,000 within-protocol run-bootstrap draws. Their status is `fixed_prediction_conditional_only`. They hold OOF models, references, scales and folds fixed. They do not estimate full retraining uncertainty or justify an unconditional superiority P value. Overlapping training sets remain dependent.

### Prespecified sensitivities

| Scope and predictor | Lines | Gain | Conditional 95% interval |
|---|---:|---:|---|
| **Mapped, original ETS2 program — primary** | **185** | **0.130921** | **[0.124191, 0.136547]** |
| Mapped, comparator-overlap removed | 185 | 0.147878 | [0.141669, 0.153906] |
| All source lines, original program | 189 | 0.131721 | [0.125157, 0.137276] |
| All source lines, overlap removed | 189 | 0.148769 | [0.142750, 0.154782] |

The all-source sensitivity admits four separately named, less-qualified `NOMATCH` lines; its line shares are 114/189 and 75/189. It does not establish their identities. References and models are rebuilt under that fixed sensitivity, rather than merely appending four test rows.

Positive gain persists after removing the 167 shared target members. This does not make the sensitivity the primary result. The larger point estimate is not permission to choose the more favorable program, and no test of the difference between these gains is claimed.

## 3. Why this is a fixed-model gain, not proof of unique information

Adding a feature at a fixed ridge lambda changes both the available representation and the geometry of its penalty. If the added feature is correlated with existing predictors, the extended model can change effective shrinkage along directions already represented by the baseline.

An exact duplicate is a useful mathematical example, **not an assertion about the observed ETS2 score**. Splitting an existing coefficient equally across two identical columns preserves predictions but halves that direction's squared-coefficient penalty. Held-out loss can therefore change without adding new input information. Standardizing each column and retaining the same numerical lambda do not remove this issue.

The actual ETS2 predictor may also add useful predictive variation. This review has not fitted alternatives to separate these explanations and does not claim that regularization explains the observed gain. The result demonstrates improvement of this **particular frozen ridge comparison**. It is not an estimate or proof that conditional mutual information `I(stimulus; ETS2_score | comparator_scores)` is positive, nor evidence of uniquely ETS2-specific biology.

Removing shared gene members does not remove expression correlation. The nonoverlap score also has its own matched background and training scaling. Persistence of its gain cannot by itself isolate a distinct ETS2 mechanism or eliminate the penalty-geometry explanation.

**Suitable result wording:** “Adding the fixed ETS2 gRNA1-down program improved held-out, line-weighted log loss relative to the prespecified four-program ridge baseline in this released macrophage stimulus-context dataset.” Avoid replacing this with “independent ETS2 information,” “ETS2-specific activation,” or “PSC validation.”

## 4. Response panel and retained flags

Each identity scope retains all 1,512 aggregate rows:

- 1,458 `available_descriptive` summaries, including raw, background and adjusted quantities;
- 54 PIC program/protocol/time rows with `unavailable_mock_identity`, not estimated effects;
- 972 protocol-specific pointwise conditional interval rows;
- 486 pooled descriptive rows with `not_computed_prespecified` intervals.

**Neither pooled response means nor pooled response medians have a confidence interval.** This does not remove the mandatory pooled predictive-gain interval. No response P-value discovery screen or simultaneous full-panel coverage is claimed.

Descriptively, the main adjusted-program medians are positive for all nine 24h stimuli in both protocols. Six-hour signs vary; IL4 and sLPS have opposite median signs between protocols. These are descriptive checks of the complete panel, not a selection of “significant” stimuli. They do not support a universal activation story. A negative adjusted response means a decrease relative to its matched expression background, not demonstrated ETS2 inhibition. A positive response is likewise not a direct measure of ETS2 activity.

Protocol differences are not identified causal library-preparation effects: their donor pools differ. Time-specific matching references also differ, so these summaries are not a formal time-interaction test. The original publication's stimulation-based expression QC and known responsiveness limit novelty and the target population. Prediction remains conditional on the released available-condition pattern, not an unobserved complete stimulus grid.

## 5. What passed, and what was not checked exhaustively

The completed frozen verifier reports:

- **80/80 fit cells verified**, covering both identity scopes and both program variants. All four endpoint summaries, cell intervals and fixed-prediction run-influence outputs passed.
- Maximum saved-model gradient norm `3.003239766319332e-7`, below the frozen `1e-6` criterion. The maximum coefficient difference from the independent refits was `1.4228528661952566e-6`, within the frozen comparison tolerance.
- **Eight source samples × 58,243 genes** checked against independently reconstructed normalization. Maximum absolute difference was **`1.7763568394002505e-15`**, not literal zero. This is floating-point-level agreement.
- **Eight fold-0 source-to-score contexts** passed: both protocols, times and identity scopes, all nine programs, and the prespecified training/test pair subsets. No program failed in those audited contexts.
- All response-panel points/statuses passed. Independently replayed response intervals cover **CIL6h only**: 54 protocol/program/quantity cells in each scope. Other response interval values were not numerically replayed.

Raw normalization was not independently reconstructed for every sample. Numerical reference/bin/draw/score reconstruction was not repeated for all folds. Other folds received identity, paired-X, scaler, model and loss checks against the saved outer-fold state. A failed primary model would not be replaced by the independent solver; no such failed fit cell occurred here.

The root separately reports an exact-endpoint clean numerical replay and 266/267 matching manifest artifact hashes, with only execution-start metadata differing as prespecified. That is a repeatability check, distinct from independent mathematics. This reviewer did not repeat that execution or independently rehash the whole expression cache; final full-file replay receipts belong to the root audit.

The existing source boundaries remain: normalization uses the qualified released feature universe, and historical program membership is recovered from the pinned complete export rather than a new reading of the absent original ZIP. This review does not change those source qualifications. The root's reported seven baseline full-repository test errors are not silently relabeled as an all-repository pass.

## 6. Biological and reporting boundary

MacroMap is healthy-derived iPSC macrophage stimulation, not a PSC cohort or an ETS2 intervention in this analysis. Related guide/enhancer programs are not independent biological replications. The result supports a bounded claim about this fixed score/model comparison in stimulus context. It does not prove ETS2 dependence, mediation, genetic allele direction, a PSC mechanism or a treatment effect.

Keep the primary result, all sensitivities, both protocols, both times and the complete response panel visible. Keep these results distinct from PSC organoid disease-response results. No clinical recommendation follows.

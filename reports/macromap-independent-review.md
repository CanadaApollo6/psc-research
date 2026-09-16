# Independent MacroMap prefreeze review

Reviewed 16 September 2026. **Ready for a separate root execution freeze, with the limits below. This review is not execution authorization.** No real MacroMap normalization, matching, program scoring, model fit or prediction was performed by this review. Its numerical checks used synthetic data. Source-only checks used the existing public-data caches.

This reviews the [fixed-program design](macromap-program-design.md), its guarded runner and the private audit schema. It does not change the programs, prediction endpoint, model, tuning rules or cohort. The completed organoid implementation and its freeze were not changed.

## 1. Reviewed artifacts and results

| Artifact | SHA-256 |
|---|---|
| `scripts/analyze_macromap_program_context.py` | `cabd1d9c12a6cd0f7f86a74dedf40526fff897b644208da41368abcf571647ef` |
| `scripts/run_macromap_program_context.py` | `dc3f976ee05ee681e98b3cdd6edf1b15e9030d9cd25f47455739486f2adf8a3c` |
| `tests/test_macromap_program_context.py` | `13e284cedbe32ea231155cdd91cef4ea7c623d183cc6dfc5bd6eec55071d3844` |
| `tests/test_run_macromap_program_context.py` | `048456a3c78b38d043bb4283ee86a634c05c7ddb1963dccd5321ed7e514aca12` |
| `scripts/verify_macromap_program_context.py` | `31b69fca6bcbe39e9204d3816203b39c5cb96762ef2e56ec9562dede77c71a7a` |
| `tests/test_verify_macromap_program_context.py` | `80c2492231cb8585a6f0bd431abcd4cb5b8518532427d7a83a45643bd49cb895` |
| `work/macromap-design/preparation-runner-v3/source-only-readiness.json` | `eae76811c5b46388965f4f2c27624e650748b9d7adc216f75bff05b389e1f507` |

Completed native checks:

- **71 targeted tests passed:** 40 method/guard tests, eight runner tests and 23 independent-verifier tests.
- **100 separate independent guard assertions passed** against the exact method/runner hashes above.
- The independent source-only audit passed against the final preparation hash. It reconstructed 58,243 features, 4,698 samples, 209 source lines, 205 mapped lines, 185 predictive lines, 57 runs, 20 fold gates, nine programs and 167 removed overlap features.
- The synthetic integration exercised both identity scopes and all 80 expected fit cells. It checked source normalization, outer-fold references, paired training/test scores, model state, expected losses, full response-panel points and complete prediction uncertainty. A missing manifest cell failed verification.

The commands used the project interpreter, `.venv/bin/python`. Final preparation records Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. Logs are ignored under `work/macromap-independent-review/`: `final-targeted-tests.log`, `source-only-final-check.json` and `guard-checks/final-review-cabd1d9.md`. This is a targeted pass, not a claim that every repository test passes. The methods worker separately records known full-suite failures caused by absent legacy caches/runtime.

## 2. What was independently checked

The verifier does **not** import the primary scoring or model helpers. Tests import them only as subjects under test.

### Source identities and fixed design

The source-only audit checks byte/hash pins, the cached released header, the original metadata CSV, static line metadata, exact treatment/control identities, same-time pairing and all expected registry keys. It reconstructs the both-time prediction cohort and whole-run hash folds. Separate `NOMATCH` source lines remain separate; they are not collapsed into a donor.

It independently joins the released feature universe to the original GENCODE v27 annotation and remaps the complete historical program-membership export. The two modality copies agree. The reference-only 45 `_PAR_Y` records remain absent. ETS2 exclusion, parent mapping gates and the comparator-overlap subtraction are preserved. The original ETS2 ZIP is absent; this is reuse of its pinned complete export, not a new reading of that ZIP. The primary preparation replays the original worksheet; this verifier uses the pinned original CSV as its metadata oracle.

The initial review used preparation hash `22ef449dc33f7b7b86dabb807c6ef5d7de5412d9863d3b0a1b7e12c63deda478`. An unfrozen preparation run replaced that summary at its default path. The hash check correctly rejected it. A byte-identical initial summary and its nine still-hash-identical artifacts were preserved under `work/macromap-independent-review/reviewed-preparation-initial/`. The final review uses the separate fresh v3 directory above, not a silent substitution of pins.

### Scoring and model mathematics

Independent checks cover:

- All released features in the CPM denominator, including measured zeros and genes outside signatures.
- Training-run, same-time controls once per line. Treatments and held-out controls do not set the reference. Reference-only lines may lack the other time and need not enter the classifier cohort.
- Stable abundance/ID bin ordering, fixed PCG64 seeds, per-bin draws without replacement, union exclusions, exact selection hashes and no pool fallback.
- Explicit matched-set means versus compressed gene-frequency weights, detection fractions and treatment-minus-own-control direction.
- Equal-line training weights and weighted population standardization, including constant columns.
- The nine-class **weighted mean** cross-entropy objective, symmetric slope penalty λ=0.1, unpenalized intercepts and centered class coefficients. The independent implementation checks analytic gradients/Hessians and uses `trust-exact` in a sum-zero class basis; the primary implementation uses L-BFGS-B.
- Independent coefficients, objectives and held-out log probabilities; missing classes, nonfinite inputs, failed fits, leaked line/run identities and altered scalers are rejected or explicitly unavailable.
- Equal stimulus weights within line/time, one equal-time average and the original protocol line shares. Unequal treatment-row, fold and run sizes do not become equal-weight replicates.
- Conditional run bootstrap by explicitly expanding sampled runs, rather than reusing the primary sufficient-sum/count implementation. Leave-run influence keeps the original protocol shares.

A key integration check requires **both training and test X to use the same outer fold's training-control reference**. A global table of each line's own-fold response is not reused to train another fold. Such reuse could expose an outer test run's controls to training preprocessing.

## 3. Issues found before freeze and their disposition

All substantive findings were sent to the root and methods author. The reviewer did not edit the primary implementation.

| Finding | Repair and independent result |
|---|---|
| Float-first parsing accepted `9.9999999999999999` as 10 and `-1e-400` as zero. | Lexical/exact Decimal checks now precede float conversion. Near-integer fractions, underflow and invalid counts fail; valid integral syntax remains supported. Synthetic regressions pass. |
| Preparation could trust a symlinked ignored anchor or overwrite an existing artifact's symlink target. | Literal private containment, nonaliased paths, fresh run directories and exclusive writes. Reproduction no longer creates or overwrites the target. |
| A dangling `.partial` symlink could redirect the cache payload. | Exclusive file-handle creation and alias checks cover the full cache namespace. Stale payloads, sidecars and reservations are rejected. |
| Concurrent writers could share a memmap inode. A failed second writer could alter the first writer's completed cache while its sidecar still said complete. | An exclusive reservation spans cache/sidecar publication; the successful marker remains. Publication cannot replace another file. In the forced race, the loser never maps a file and the winner's hash remains valid. |
| `abspath` could erase `link/..` during a check while a writer opened the original path. | Literal `..` is rejected before normalization. Text, CSV, NPZ and publication regressions pass. |
| Fold-gate status could overwrite actual failed-fit status in diagnostics. | Design-gate status and actual fit status are now distinct. The complete manifest retains failed cells. |
| A normalized source hash was recorded without being compared again before scoring. | The runner checks source identity/stability, compares the normalized-source pin before matching/fitting and rechecks plan, code, runtime, sources, preparation and cache before complete publication. Source-change and fail-fast tests pass. |

No substantive issue remains in this scoped review. This is **not** a general filesystem sandbox. Low-level helpers accept ordinary paths for synthetic use. The supported real-data route is the separately frozen CLI, which enforces fresh ignored output containment. No unchanged legacy source value was corrected or overwritten by this review.

## 4. Final response-interval clarification

Before any real scores, the root finalized the narrower rule:

- Protocol-specific response means and medians retain pointwise run-bootstrap intervals.
- **Pooled response mean and median are descriptive only. Neither has a pooled confidence interval.** Their interval status is `not_computed_prespecified`.
- The mandatory pooled **predictive-gain** interval remains unchanged, with fixed protocol line shares and equal time weights.

This final clarification supersedes an intermediate request for a fixed-share pooled-response mean interval. No weighted-median convention, response-P discovery screen or outcome-dependent endpoint selection is introduced.

All intervals condition on frozen OOF references, scaling, models and folds. They do not capture full procedure/retraining uncertainty. Overlapping training sets are not independent fold replicates. Pointwise response intervals do not provide simultaneous full-panel coverage.

## 5. Frozen read-only verification scope

After a separate root freeze and completed execution, the verifier can read the `macromap-private-audit-v1` bundle. It authenticates the externally supplied plan, verifier and completion-receipt hashes; checks code/software/preparation/source/cache pins and chronology; rejects incomplete inventories; and does not alter an input bundle.

The scope is fixed before effects:

1. **Raw normalization:** eight source sample IDs, chosen by ascending SHA-256 of `macromap-independent-v1:` plus the original ID. Reconstruct their counts and totals over every released feature and compare the normalized cache. This does not independently validate every count token.
2. **Source-to-score:** outer fold 0, both protocols, both times and both identity scopes. For each context, choose up to four training and four held-out predictive pairs by that prefix plus compact JSON `[protocol,run,line,integer_time,stimulus]`. Reconstruct the full control reference, bins, all nine programs' draws/weights, sample scores and paired X. This is eight contexts, not every numerical reference in the run.
3. **All model cells:** require all 80 expected scope/variant/protocol/time/fold entries. Check source axes, training/test identities, controls, saved paired X against that outer-fold score state, scalers, independent fits, probabilities and losses. Failed cells remain failed; the independent solver never replaces a failed primary fit.
4. **Prediction outputs:** independently reconstruct the primary points, all four protocol/time cells, conditional intervals and leave-run influence for every scope/variant. Missing expected predictions make an endpoint unavailable.
5. **Response outputs:** check all panel keys, points, ranges, signs, source statuses and leave-out ranges against all saved line responses. Independently replay interval values for **CIL at 6h**, both protocols, all nine programs and all three quantities, in both scopes. Other response interval values are not numerically replayed; their shared algorithm is covered by the synthetic tests. No pooled response interval is accepted.

The response-interval subset is a computation check, not a filter on the reported panel. All stimuli, times, programs and unavailable rows remain required. Source-to-score reconstruction outside the declared subset and full-pipeline bootstrap refitting are not claimed.

Example postfreeze command (hash placeholders must come from the root's records):

```text
.venv/bin/python scripts/verify_macromap_program_context.py \
  --preparation work/macromap-design/preparation-runner-v3 \
  --preparation-sha256 eae76811c5b46388965f4f2c27624e650748b9d7adc216f75bff05b389e1f507 \
  --execution-directory <frozen-private-output-directory> \
  --plan <root-freeze.json> --plan-sha256 <frozen-plan-sha256> \
  --execution-complete-sha256 <completed-receipt-sha256> \
  --verifier-sha256 31b69fca6bcbe39e9204d3816203b39c5cb96762ef2e56ec9562dede77c71a7a \
  --output work/macromap-independent-review/postfreeze-verification.json
```

Without execution arguments, the CLI performs source-only checks. The real-output mode has not been invoked. A new root freeze is still required; changing a hash to bypass a failure is not verification.

## 6. Scientific boundary

These are healthy-derived iPSC macrophage stimulation data. The study's prior expression QC and published responsiveness are already known. The proposed endpoint measures incremental stimulus-context information beyond four fixed comparators. It is not PSC validation, direct ETS2 perturbation here, proof of ETS2 activity/dependence or mediation, a genetic allele-direction result, or clinical evidence.

Positive or negative findings must preserve protocol/time and sensitivity results and the full response panel. A conditional interval crossing zero is imprecision, not proof of exact redundancy. No result authorizes a clinical recommendation or outcome-driven changes to the frozen method.

# Independent PSC organoid method review and verifier

**Status: pre-freeze, synthetic verification complete. 16 September 2026 UTC.**

This review covers the [protocol](../docs/psc-organoid-il17-protocol.md), [method implementation](../scripts/analyze_psc_organoid_il17.py), [method note](psc-organoid-method-design.md), and [input qualification](psc-organoid-input-qualification.md). It does **not** report real organoid effects. No real-data verification receipt has been created. The previous frontier review is unchanged.

## Decision and scope

The independent verifier is ready for a frozen run:

- `scripts/verify_psc_organoid_il17.py`
- `tests/test_verify_psc_organoid_il17.py`

The verifier imports neither the analysis module nor PyDESeq2 helpers. It rebuilds the core calculations from source count/metadata tables and exported native fit arrays. NumPy and SciPy are shared numerical dependencies, not independent software ecosystems. The integration tests call the analysis module **only to produce synthetic outputs**; expected mathematics come from the verifier and separate library/manual references.

A real run still requires the root's reviewed raw-count decision, final code/software/source pins, UTC freeze, primary fit and independent verification. A numerical pass is not evidence that the observational comparison is unconfounded, that small-sample Wald P values are calibrated, or that an IL-17 intervention benefits patients.

## Source and biological-unit review

The fixed panel has eight donors, four PSC and four non-PSC procedure controls, with two cultures per donor. The 16 libraries and 46,343 retained cells are not independent donor replicates. The population token is `all_published_retained_organoid_cells`; it does not assert independently proven lineage purity.

The source qualification distinguishes two different quantities:

- The 21,562-feature merged integer matrix matches SCT margins. It must not be supplied as original RNA counts.
- The selected original Cell Ranger RNA exports contain 36,601 features for the same 46,343 retained cells. Their full-feature total is **1,128,543,232**, versus **1,128,472,969** in historical RNA metadata. The **70,263** difference is preserved. No guessed feature filter should force agreement.

The verifier checks frozen file bytes/hashes, full feature identity and annotation preservation, exact integer tokens, exact library totals, donor pairing, population, software versions and method settings. It requires the reviewed raw-RNA gate. It does **not** infer assay identity from integers, repeat sparse cell-to-library aggregation, or recover the historical reference bundle/export command. Those remain source-qualification responsibilities. The frozen root decision must explicitly distinguish admissibility of the original full raw universe from reproduction of the historical RNA metadata universe.

## Independent mathematical checks

| Component | Independent reconstruction |
|---|---|
| Eligibility | At least 10 counts in at least four of 16 libraries; all source rows retained. |
| Normalization | Eligible genes positive in all 16 libraries; sorted log ratios and their middle element(s), followed by exponentiation. An even basis uses the mean of the two middle **log** ratios, not the arithmetic median of raw ratios. No fallback. |
| Design | Explicit intercept, seven sorted-donor indicators, treatment and PSC×treatment. Rank 10, six residual df, leverage 0.625; a disease main effect is redundant. Cell counts never become model replication weights. |
| Native NB reconstruction | Natural-log coefficients, final native dispersion and independently computed size factors produce fitted means. Solve linear systems for `H M H`, with `M = X′WX`, `W = mu/(1+dispersion*mu)`, and `H = (M + 1e-6 I)^−1`. |
| NB contrasts | Interaction, non-PSC treatment and PSC treatment. The PSC variance includes both variances and **twice their covariance**. Rebuild log2 coefficients/SE, signed Wald statistics, two-sided normal tails and 95% normal intervals. |
| Families/missingness | Rebuild BH and BY over the whole eligible family. Missing or count-boundary-masked P values remain unavailable; only correction inputs become 1. Native package diagnostics remain separate. No ordinary absolute tolerance excuses incorrect tiny P/q values. |
| Donor changes | `log2(CPM+1)` with totals over all unfiltered raw features, then one exact IL17-minus-CNT change per donor. |
| Welch | Equal-donor four-versus-four comparison, sample variances, actual Satterthwaite df, t tails and t-based intervals. Positive-variance df is 3–6. Zero/numerically degenerate variance remains unavailable; no case deletion. |
| Leaveouts | All eight complete-donor omissions on the same paired-change features, all estimates, minima/maxima and sign counts. No NB refit and no new filter. |
| Allocations | Independently enumerate all 70 four-label bit masks, in the documented order. Rebuild **unstudentized** mean differences, observed labels, complements and tails using `1e-12 + 1e-12*abs(observed)`. Zero effect has fraction 1; minimum is 2/70. No FDR on these descriptive fractions. |
| Diagnostics | Compare native flags, dispersion trend label/fallback flag, full Cook-distance summaries and diagnostic F(10,6) cutoff. No automatic Cook exclusions. The allowed trend labels are parametric and built-in mean fallback. |
| Integrity/privacy | Check public manifests, original annotations, fixed row coverage and private array axes. Public outputs are disjoint from raw/work inputs. Verification is read-only; an optional receipt is created exclusively and cannot overwrite inputs/results. |

The NB interaction is a log ratio of treatment ratios. The Welch contrast is a mean difference of donor log-CPM changes. They are different estimands. Separate within-group significance does not test their interaction. Donor blocking removes baseline levels, not disease-associated response confounding by age, source, clinical history or culture selection.

## Pre-effect boundary problem found and corrected

Full design rank does not ensure a finite, identified response estimate for every gene. The root requested two targeted synthetic cases before freezing:

1. **Both arms are zero in one disease group.** That group's response has no likelihood information. Its response and the interaction are unidentified.
2. **One complete disease×treatment arm is zero.** The corresponding unpenalized log response is on an infinite-estimate boundary. Ordinary interior Wald inference is not justified merely because software returns finite numbers.

A separate direct PyDESeq2 0.5.4 simulation used 96 genes, seed 72865, rank 10 and the same design. With both non-PSC arms zero, the package returned finite interaction P≈0.9893 and non-PSC P≈0.9964. With the PSC CNT arm zero, it returned interaction log2FC≈11.3965, SE≈1.7163 and **P≈3.13×10⁻¹¹**; the PSC-response P was ≈1.23×10⁻¹⁰. Both examples reported LFC convergence. These are **synthetic software demonstrations**, not biological results. Numerical regularization/clamping and a convergence flag do not supply missing likelihood information.

The method worker added, before any real effects, the following narrow count-support rule. The verifier derives it independently:

- A within-group response needs positive pooled counts in both that group's CNT and IL17 arms.
- The interaction needs all four pooled arms positive.
- Affected main coefficient, SE, statistic, P and interval fields are unavailable with status `unavailable_boundary_or_unidentified_group_log_response`.
- The fixed gene family is unchanged. Those missing P values contribute 1 only to correction; their q values are 1.
- The raw package values remain in `pydeseq2_unmasked_*`; required zero arms and pooled counts are explicit.
- The other group's response can remain supported. It is not masked merely because the opposite group has no counts.

An individual zero-total donor pair is flagged, not removed. Its donor intercept has a nuisance boundary, but other donors can still identify the common group response. It is not equivalent to an absent entire group/arm and does not trigger a new gene filter. This narrow rule is **not a proof that every retained gene satisfies every regularity assumption**. Small-sample estimation and remaining boundary/calibration limits remain explicit.

Synthetic tests now verify both masks, the unaffected group, unmodified raw package diagnostics, the individual-zero-pair flag and the unchanged multiplicity denominator. No real-data pattern motivated the correction.

## Conditional fixed-dispersion optimization

The verifier selects at most 16 genes **from counts and exact IDs before reading fitted effects**:

1. Keep fixed-eligible genes with at least 10 counts in every library.
2. Hash the UTF-8 concatenation of `PSC-organoid-conditional-NB-v1:` and the exact unchanged source ID.
3. Select the first 16 by ascending SHA-256, with exact ID as the deterministic tie-breaker.

For those genes, an independent implementation optimizes the fixed-dispersion NB likelihood with the pinned 1e-6 coefficient ridge. It uses an independent count-based least-squares start, analytic score and observed Hessian, and SciPy `trust-exact`. Synthetic finite differences and SciPy NB log probabilities check the objective derivatives and likelihood differences.

This comparison is restricted to interior native fits with all fitted means above 0.5 and coefficient magnitudes below 30. Other selected fits are recorded as unavailable or outside scope, never replaced by more favorable genes. The optimizer must reach Newton decrement ≤1e-5. Agreement requires twice the objective improvement ≤0.001 and each of the three contrast differences ≤0.02 of its native SE; objective numerical slack is 1e-6. These are prospective numerical tolerances, not biological significance thresholds.

**This does not independently estimate gene-wise dispersions, fit the trend, estimate the prior or repeat MAP shrinkage.** The complete Wald reconstruction is conditional on supplied native dispersions and coefficients; the small optimizer panel adds conditional coefficient-fit evidence. Cook arrays are checked for consistency of summaries, not independently regenerated from the whole dispersion pipeline.

## Tests completed and remaining execution

```text
.venv/bin/python -m unittest tests.test_verify_psc_organoid_il17 -v
.venv/bin/python -m unittest discover -s tests -v
```

- All **24 independent-verifier tests pass**. A synthetic 112-feature run retains 110 eligible features, verifies the count-boundary masks, reconstructs every public/native calculation and passes all 16 prespecified conditional fits. Tests also detect a changed private donor value and public file corruption.
- The full repository run executes **460 tests** and has the same **seven historical errors**: four missing ETS2/GSE84161 cached-input errors and three missing GSE84161 R-runtime errors. There are no verifier test failures. This is not a clean full-suite pass.
- The reference Welch test emits an expected precision-loss warning for its deliberately constant-group case. The explicit unavailable-variance and one-constant-group rules are tested separately.
- No catalog parser, prior report, frozen prior output or source file was edited for this task. No network request, dependency change, Git operation or researcher contact was needed.

After the final freeze and primary run:

```text
.venv/bin/python scripts/verify_psc_organoid_il17.py \
  --plan config/psc-organoid-plan.json \
  --plan-sha256 <frozen-plan-sha256> \
  --output-dir <existing-primary-public-directory> \
  --work-dir <existing-primary-private-work-directory>
```

The plan must pin both the primary script and verifier before effects. Pin the tests, review/protocol, acquisition code and dependency lock too. `--structure-only` checks inputs, eligibility, design and pre-effect selection without reading/rebuilding effect outputs; it still requires a frozen plan and reviewed raw gate. `--verification-json <new-path>` writes a receipt **only when requested**. Without that flag, the verifier prints only an aggregate completion summary.

Real-data numerical agreement, real fit warnings, conditional-panel coverage and clean-input replay remain unverified at this pre-freeze stage. The final analysis must preserve negative, unavailable, opposing and donor-sensitive results. Numerical agreement is not independent biological replication, an equivalence test, a causal disease experiment or clinical evidence.

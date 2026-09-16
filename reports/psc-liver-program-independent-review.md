# Independent review: pediatric liver fixed programs

## Decision and boundary

**Accepted for prospective method/code preparation. Not an execution approval.**
The final reviewed implementation passes the source-only, synthetic numerical and
workflow-guard checks below. Root must create and approve a separate execution
freeze before opening real expression payloads. No real liver count row,
normalization, score, effect or hypothesis test was read or computed by this review.

The accepted question is deliberately narrower: **ten unsigned programs × three
contrasts = 30 Core hypotheses**. These are the nine original/nonoverlap programs
plus the official ECM set. The epithelial IL17A-alone response remains **HOLD**
outside this panel. This does not answer the blocked epithelial response-transfer
question. It is not independent mechanistic confirmation of MacroMap. All 27
specifically protected metadata/code/report artifacts, including MacroMap, remain
unchanged. There was no network, install, researcher contact or Git operation.

The machine-readable acceptance is
`work/psc-liver-program-independent-review/accepted-review.json`. Its
`real_execution_reviewed` and `real_execution_authorized` fields are false.

## Fixed method and interpretation

- Use all **41,096 original rows** in each library total. Core, Strict and program
  subsets cannot replace that universe. Validate exact nonnegative integer counts
  and safe totals before float conversion. A zero library or invalid input stops
  the run; do not delete, round, replace or repair source values.
- For each mapped gene, compute `log1p(1e6 * count / full_library_total) / log(2)`.
  The score is the positive equal-member mean after the fixed annotation gate.
  Keep measured zeros and tiny values. Missing members are not measured zeros.
  There is no abundance matching, variance scaling, learned weighting or
  expression-based membership selection. Original companion labels are not inverted.
- The focal endpoint is ETS2 gRNA1-down. **PSC minus AIH** is primary; PSC minus ASC
  and ASC minus AIH are fixed secondary contrasts. Welch tests the difference of
  patient-sample score means, with sample variances, actual Satterthwaite df,
  two-sided P and pointwise 95% t intervals.
- Root approved ordinary Welch when exactly one arm has zero variance and the
  other has positive variance. Flag the constant arm. Detect exact constancy
  before reduction rounding. Both-zero variance retains only the descriptive
  difference. Nonrepresentable variance, nonfinite inference or tail underflow
  also retains a finite descriptive difference with explicit unavailable
  inferential fields. Never manufacture observed P=0/1 or epsilon variance.
- BH/BY use the **complete 30-test Core family**, including unavailable tests.
  Only correction inputs replace missing P with 1; public raw P and q remain null
  on unavailable tests. The full numerical correction vectors remain private.
  Strict retains every status row but has **no discovery q values** and cannot
  replace a failed or less favorable Core result.
- Bootstrap **10,000** patient-sample resamples with **PCG64 seed 2026091602**.
  Group order is **PSC, ASC, AIH**, retaining original within-group source-column
  order. Share each draw across programs, interfaces and contrasts. Use linear
  2.5/97.5 percentiles; no bootstrap P family. Retain and flag degenerate intervals.
- Keep all omissions privately. Public stability metrics count only involved
  units: **47/34/47**, not all 64. Unchanged third-diagnosis rows cannot inflate
  robustness counts. Scores and memberships stay fixed; no discovery reselection.
  Omission ranges are not confidence intervals.

These are unadjusted, observational **whole-biopsy transcript-abundance
composites**. The article states 64 patients, with 17 PSC, 17 ASC and 30 AIH source
columns, but an independent person/visit crosswalk and sample-linked clinical,
stage and batch keys are absent. Root accepted that explicit independence
assumption; the columns are not independently verified unique donors. Repeated
persons would require a different analysis. ASC stays separate from PSC.

A score is not absolute RNA concentration, a summed-program fold change,
cell-intrinsic pathway/TF activity, protein abundance, fibrosis stage, causality or
a clinical outcome. Relative normalization, the pseudocount and cell composition
matter. Different programs are not interchangeable activity scales. Shared genes,
patient samples and contrasts are not independent replications. Nonoverlap removes
specified genes; it does not prove ETS2 specificity.

Welch and percentile bootstrap intervals are **pointwise, not simultaneous or
causal**. Degenerate bootstrap draws do not establish zero population uncertainty.
BH depends on suitable dependence assumptions; BY still requires valid constituent
P values. Neither correction repairs confounding or unverified independence.
Current Ensembl116 identity does not reconstruct the historical GRCh38 counting
GTF or assignment settings. GSE159676 remains HOLD; NoPSC atlas remains PARK.

## Source-only membership result

The official ECM membership is MSigDB **M610 / R-HSA-1474244**, version
**2026.1.Hs / Reactome 95**. Its denominator is **321 unique official symbols**;
340 original Ensembl IDs are provenance, not additional members or weights.
It is an **unsigned ECM-organization pathway**, not a directional fibrosis or
IL17A-response signature. **Literal ETS2 is absent**; no new target exclusion was
applied.

| Interface | ECM mapped | Original denominator | Annotation gate |
|---|---:|---:|---|
| Core exact symbol–Ensembl bijection | 301 | 321 | Pass |
| Strict additional reciprocal Entrez | 299 | 321 | Pass |

Twenty symbols have multiple current Ensembl links under both interfaces. Strict
adds one missing Entrez relation and one nonreciprocal Entrez exclusion. All 642
source-member/interface statuses and indexes are retained. No alias, alternative
source-ID bridge, favorable subset denominator or missing zero was used. These are
identity results, not measured expression or disease effects.

Every original ordered membership and all 20 compiled program/interface records
agree with the independent source comparison. The nonoverlap set still subtracts
only the original four comparator programs, **not ECM**. Its gate is the original
parent coverage gate plus ten retained genes, not a new 80% gate on the reduced set.
The canonical compiled SHA-256 is
`17a291f47fbb0ec0e815ac2f592ba4d0f657f3f5e3cbe4689dae5ecc44c90698`.

## Exact audit scope and results

The initial scope preceded testing. Separate source and method amendments record
root's later ten-program/30-test decision and fixed bootstrap/omission rules,
always before real effects. Source-only checks reused the closed identity audit;
they did not rerun its historical raw-source reconstruction.

| Final check | Result |
|---|---:|
| Independent verifier unit tests | **61 passed** |
| Author method/source/guard tests, rerun natively | **109 passed** |
| Complete final artificial numerical bundle | **584,546 scalar comparisons passed** |
| Boundary and full-size approved-RNG checks | **62 checks passed** |
| Independent source compilation comparison | **11,372 checks passed** |
| Negative/source-only entry and alias guards | **42 checks passed** |
| Temporary synthetic output/privacy writer | **51 checks passed** |
| Specifically protected pre-existing artifacts | **27/27 unchanged** |

The full artificial numerical bundle uses a 40-row × seven-unit matrix, ten
unsigned programs, both interfaces, all 60 contrast/status rows and the full
30-test Core correction family. It includes constant, measured-zero, unavailable
coverage and Strict-only unavailable endpoints. Every normalization value, score,
Welch field, shared bootstrap contrast, percentile interval, available unit
omission and public relevant-unit summary is compared. Denominator, score, q,
draw and missing-row mutations are detected. A second fresh native run gives a
**byte-identical numerical verification receipt**.

A separate 64-unit artificial design checks **640,000 bootstrap indexes** under
the approved seed/order, **30,000 constant bootstrap differences**, and public
omission denominators 47/34/47. Boundary checks cover exact equal/different
constants, one-zero-arm Welch, positive variances down to about `1e-300`, missing
inference, tail underflow, and unsafe integer/float boundaries. The extreme
`1e-200` score-scale fixture has unrepresentable variance: the author reports
explicit numerical unavailability, retains the finite difference and does not
label the groups exact-zero. It is not an observed biological null.

The independent verifier imports **no author helper**. It uses separate Python
`fsum`, scaled-SD and adjustment calculations. The separately pinned adapter
imports the author only to produce artificial outputs. SciPy distribution
primitives and NumPy PCG64 remain shared dependencies; this is not an independent
special-function or RNG-library implementation. Numeric tolerances are relative
`1e-10`, absolute `1e-12`; identities, row inventories and statuses are exact.
The tail/constant boundary probes use explicit additional assertions.

Entry checks refuse draft/unapproved execution, changed pins, matrix-as-plan or
matrix-as-source substitutions, unqualified source adapters, symlinks, known
payload hardlink aliases, reused directories and public destinations. Refusals
precede real payload reads. The writer test used **temporary synthetic inputs and
injected temporary authorization only**. It verifies owner-only run directories,
exact private plan bytes, complete artifact hashes, private omission inventories
and no individual IDs in the public aggregate. It does not certify a positive
real-data authorization. Workflow markers are not cryptographic proof of identity.

No real numerical replay or biological result is certified here. The verifier CLI
has only `--source-only` and `--synthetic-only`; actual run-output auditing needs a
separately authorized, pinned execution receipt and an explicit native-output
adapter. The full repository suite was not rerun by this review, and no clean
whole-repository test claim is made.

## Repairs resolved before acceptance

The first author draft invented effect `4.44e-16` and P about `8.70e-6` for
`[2.1]*30` versus `[2.1]*17`. Exact-constant reduction now returns effect zero with
unavailable inference. Scaled Welch ratios fix representable tiny-variance df;
numerical inability is distinct from exact zero. Stable log transforms and direct
sign comparisons replace cancellation/product-underflow risks. The final code
also implements root's in-contrast-only omission summaries and unavailable
inference rule; the earlier five-versus-seven synthetic failure is resolved.

Source/plan whitelists, qualified-manifest semantic binding, exact integer totals,
private-path ancestry, metadata alias guards and truthful failed-stage receipts
were added before acceptance. The reviewer did not patch author-owned files.
The independent verifier also corrected its own float-boundary guard and initial
schema/tolerance assertions. Earlier logs, pins, failed probes and the iterative
review text remain in the ignored review directory; they are development evidence,
not unresolved current blockers or real RNA findings.

## Final pins and reproduction

| Artifact | SHA-256 |
|---|---|
| `scripts/analyze_psc_liver_programs.py` | `f7d5e49d2309c2e21a2ea20d0a7d1b8bda8643e140b0c781db3edafd34fd5ec7` |
| `tests/test_psc_liver_programs.py` | `0a0fa5f0c3f454514d22057d63c970a809a548f5a8ef9cdf9106dc1fa8c086a5` |
| `scripts/verify_psc_liver_programs.py` | `82403f2cb63ae1eeb7e4f4a1232b21ae6772515d0e1b3ead1ca7397737a62391` |
| `tests/test_verify_psc_liver_programs.py` | `64fbcfd02c3a7d179d98fdc8f13b37d2160a274db02c7f3ca7cfc871ea5b0edc` |

Native environment: Python **3.12.14**, NumPy **2.5.3**, SciPy **1.18.1**.
Complete source, adapter, audit and preservation pins are in the acceptance JSON.

```text
.venv/bin/python -m unittest discover -s tests -p test_verify_psc_liver_programs.py -v
.venv/bin/python -m unittest discover -s tests -p test_psc_liver_programs.py -v
.venv/bin/python scripts/verify_psc_liver_programs.py --source-only
.venv/bin/python scripts/verify_psc_liver_programs.py --synthetic-only
.venv/bin/python work/psc-liver-program-independent-review/compare_synthetic.py
.venv/bin/python work/psc-liver-program-independent-review/boundary_probe_final.py
.venv/bin/python work/psc-liver-program-independent-review/guard_probe.py
.venv/bin/python work/psc-liver-program-independent-review/output_guard_probe.py
.venv/bin/python work/psc-liver-program-independent-review/compare_source_compilation.py
```

**Remaining gate:** root-owned, hash-bound execution approval. No current file or
successful synthetic check supplies that permission, observed RNA evidence,
independent biological replication or a clinical recommendation.

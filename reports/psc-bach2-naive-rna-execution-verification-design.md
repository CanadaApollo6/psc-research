# BACH2 real-output audit: pre-effect design

**Preparation and synthetic validation only. This document grants no expression access.**

## Fixed scope

The separate adapter is `scripts/verify_psc_bach2_naive_rna_execution.py`. It audits only the accepted-for-preparation endpoint in `psc-bach2-followup-design-review.md` and the root method decision. It does not relax or replace the existing metadata/synthetic verifier.

- One exact source triple: `ENSG00000112182` / `BACH2` / `Gene Expression`.
- Source feature row **12,314, one-based**. This is not a genomic coordinate.
- Exact states `C0: CD4+ TN RTE` and `C1: CD4+ TN mature`: 9,320 + 7,339 = 16,659 cells within the fixed 55,460 retained cells.
- Eight source-labeled male PSC donors; four `SNP` and four `noSNP`. Sample identity remains in every barcode join.
- For each donor, independently stream BACH2 units B and **all 36,601 RNA-row units R in the same selected cells**, excluding all 38 ADT rows. No prior technical margin supplies R.
- `Y = log2(1 + 1000000 B/R)`. Equal-donor mean SNP minus mean noSNP. One nominal actual-df Welch interval/P. All eight whole-donor leaveouts are point estimates only.

Zeros remain. All eight donors are required. Bad membership, axes, source integrity, R or source units make the entire endpoint unavailable. One constant arm with positive SE is valid. Both constant arms give a descriptive point only. Native numeric unavailability does not become a forced P=0 or a perfectly certain interval.

There is one scientific primary, one independent verification and a mechanical replay. These are not independent biological samples. No new model, gene, state, ADT endpoint, clinical join, normalization, test family or biological claim is added.

## Default and future authorization

Default mode authenticates **metadata only**. It requires the compiled plan's external SHA-256. It reads only the three allowed compiled-plan files, fixed source annotations and individual integrity/acquisition receipts, code, the ignore policy and native runtime files. It does not open archives, matrices, numeric result JSON, closed cohort quantity JSON or technical margins. It does not write result files.

Example, after the final preparation plan is available:

```text
.venv/bin/python scripts/verify_psc_bach2_naive_rna_execution.py \
  --plan-directory work/psc-bach2-naive-rna-plan/candidate-v4 \
  --plan-sha256 <externally-verified-final-plan-sha256>
```

The accepted author, independent verifier, adapter, tests, this protocol and test receipt must be pinned by root **before effects**. An evolving preparation candidate is not final authority.

Future real audit needs **both** `--real-audit` and `--acknowledge-root-authorized-real-audit`. It also requires all of these external arguments:

- Exact plan directory and plan SHA-256.
- Completed author run directory under `work/psc-bach2-naive-rna-plan/runs/<simple-name>`.
- Separate author authorization path/SHA-256.
- Separate root audit seal path/SHA-256.
- Exact last-written `completion.json` SHA-256.
- A fresh audit destination under `work/psc-bach2-naive-rna-execution-review/<simple-name>`.

The two root JSONs must be distinct, simple-name files directly below `work/psc-bach2-naive-rna-freeze/`. Their role and size limits are checked **before opening**. A matrix cannot be passed as an authorization, seal or plan. The author authorization retains its native exact schema; it cannot contain arbitrary new fields. Root will create separate primary/replay author authorizations outside the author's worktree. Their destination/auth hashes differ; their scientific plan does not.

The adapter never writes a usable real author authorization or root seal. `author_authorization()` and `seal_bindings()` define expected schemas for comparison; the synthetic tests use them only in isolated temporary roots with fabricated source pins. Neither schema construction nor metadata success grants access.

### Exact root audit seal contract

New plan/result/root-contract JSON uses sorted keys, two-space indentation, finite JSON and a final newline. Duplicate keys and noncanonical contract bytes fail. Original cached sources retain their original bytes; the candidate's original-file SHA and canonical-object SHA are pinned separately. The root seal has exactly these fields:

- `schema_version: 1`, `status: ROOT_AUTHORIZED_ONE_FIXED_BACH2_REAL_OUTPUT_AUDIT`, `authorizer: root`, `real_matrix_audit_authorized: true`, `one_fixed_endpoint_only: true`.
- `plan_directory`, `run_directory`, `audit_output_directory`, `authorization_path`.
- `plan_sha256`, `selection_sha256`, `source_manifest_sha256`, `implementation_sha256`, `runtime_sha256`, `public_schema_sha256`, `authorization_sha256`, `candidate_sha256`.
- `completion_sha256`, `protocol_sha256`, `adapter_sha256`, `adapter_tests_sha256`, `design_sha256`, `verifier_sha256`, `verifier_tests_sha256`.

All values must match the independently authenticated current files and explicit arguments. This post-completion audit seal must use the already accepted pre-effect implementation/protocol pins; it is not permission to change methods after observing results. Root's outer freeze also binds aggregate test receipts. No seal is created by this preparation task.

## Authentication and read order

1. Check acknowledgment, required external hashes and dedicated input/output roles.
2. Authenticate root seal and author authorization bytes. Check fresh ignored output policy.
3. Authenticate the exact compiled-plan inventory: `plan.json`, `selection.private.json`, `public-schema.json`. Reject missing/extra members. Independently rederive source donor mapping, full axes, exact prefixes, retained keys, selected keys and state membership.
4. Verify the fixed annotation/receipt allowlist, all source/selection/schema/candidate hashes, implementation pins and native runtime snapshot. Matrix pins are compared to the individual closed integrity receipts without reading their bodies. No broad qualification quantity JSON is opened.
5. Check the root seal and native author authorization exactly. Authenticate the complete author inventory: `authorization.used.json`, `completion.json`, `measurements.private.json`, `public-aggregate.json`.
6. Authenticate the externally pinned completion record and all **three actual opaque payload hashes and byte lengths**, before parsing either numeric JSON. Completion binds the exact sorted inventory, output directory and complete native provenance. It must not predate any payload. No missing, extra, pending, failed or stale run is accepted.
7. Exclusively create the fresh audit output. This consumes the destination even if later checks fail. Only now open the eight real decoded MatrixMarket streams, if future root authorization exists.
8. Independently parse every coordinate to EOF. Require exact header/exporter/dimensions, all admissible positive exact integer tokens, strict column-major uniqueness, nnz, bytes and SHA-256. Bounds are 50 million coordinates per matrix, 100 decimal digits per positive value and 512 bytes per line. Check all input coordinates, including ADT/unselected entries, but accumulate only selected-cell RNA R, BACH2 B and BACH2-positive selected-cell counts. No ADT or unselected-gene marginal endpoint is accumulated.
9. Require every stream's source integrity and valid R before computing Y. Compare every private donor record, all eight private leaveouts and all allowed public fields to the independent reconstruction.
10. Recheck metadata, code, root seals, complete inventory and all payload pins before the aggregate PASS receipt. Any error leaves a generic unavailable receipt, never a reduced-donor estimate or rescue.

Canonical relative paths, link-free components, regular single-link input files, distinct inventory inodes, fixed namespaces and exclusive fresh output directories reject symlinks, hard-link aliases, traversal, reused destinations and source overwrites. These are protections for this **bounded pinned-file workflow**, not universal security or protection against an adversary controlling the process/kernel/filesystem.

## Numerical criteria fixed before values

Root approved these acceptance criteria before value access:

| Quantity | Required agreement |
|---|---|
| Source units/counts, joins, axes, hashes, schemas, inventory, status and availability | Exact |
| Zero/nonzero, sign of Y/means/Delta/LOO/SD/SE, exact group constancy | Exact; no epsilon |
| Nonzero Y, group means, fractions, Delta, LOO, sample SD, SE, df and P | Relative error at most **1e-10**, with **no absolute floor** |
| Each interval endpoint | `(bound - independent Delta) / independent SE` agrees with the corresponding signed independent critical value to relative **1e-9** |
| Probability | Strictly positive and at most one if reported; exact Delta=0 with available inference requires P=1 |

An absolute tolerance must not hide an extremely small nonzero point, variance, SE or probability. Reported interval bounds must strictly straddle the independently reconstructed center when inference is available. Degrees of freedom must be within [3,6]. The higher arithmetic precision is a calculation check, not hundreds-digit scientific precision.

The independent reference evaluates exact rational `Q = 1 + 1000000 B/R` using **4,000-digit direct Decimal natural logarithms**, not the author's near-one atanh series. Group means use products of Q. Delta uses the exact ratio of group products; 3/4 leaveouts use the exact powered product identity. Exact rational equality determines cancellation and sign before logarithms. Sample variances use centered exact Fraction moments of independent high-precision Y, not float64 donor scores or the production parser/calculator.

The source bounds imply R < 5×10^107 and Q ≤ 1,000,001. Distinct donor Y can differ far below float64 resolution. For a 3/4 leaveout, a conservative numerator/denominator product bound is under about 2,664 decimal digits, including the Q numerator scaling. Four thousand direct-ln digits cover this cancellation domain with a large margin. No finite absolute epsilon is used to declare a constant arm or a zero effect.

For P, the independent reference evaluates the Student-t incomplete-beta identity with a 100-digit Decimal hypergeometric series. Symmetry/complement puts its series argument at most 1/2. The series has a 512-iteration hard bound and a relative 1e-80 term stop. Its positive tail is retained even below the float64 subnormal range. The log-beta normalizer uses float64 `math.lgamma`; the P acceptance tolerance is 1e-10, not the arithmetic working precision. Near P=1, complement formation uses the statistic directly to avoid cancellation.

The independent critical value inverts `scipy.special.betainc` with `scipy.optimize.brentq`; it does not use `scipy.stats`. Native `scipy.special.stdtr`/`stdtrit` only corroborate the author's finite/underflow availability status. They do not supply the reference P or CI. A native underflow can justify an unavailable P only while the independent mathematical tail remains positive. A native unavailable quantile leaves point-only output, not an alternate interval.

**Shared dependencies:** Python Decimal/Fraction/math, the pinned Python runtime and native NumPy/SciPy environment are shared with the author. Special-function availability checks share underlying native functions. Independence is in source selection/reconstruction, parsing and calculations; it is not complete software-stack independence.

## Output and privacy

`measurements.private.json`, every donor quantity/identifier/join and all eight donor-indexed leaveouts stay private. The adapter keeps the independent reconstruction in memory. Its final `aggregate-verification.json` contains only aggregate comparison counts, PASS status, approved provenance and file/protocol pins. It does not repeat means, effects, P, intervals, donor values, identifiers or LOO values. CLI errors contain no raw exception/source text.

The public scientific payload is checked against a separately embedded recursive allowlist, including nested scalar types. Extra or nested private fields fail. Public group summaries, contrast, complete LOO count/range/sign summary, QC and provenance must exactly project the private outputs already checked numerically.

## Synthetic validation and handoff

Owned tests fabricate native-schema sources: 36,601 RNA plus 38 ADT rows, the target at row 12,314, eight donor/sample labels, 59,315 source cells, 55,460 retained and 16,659 selected cells. Temporary-root author execution produces actual native payloads and its last-written completion. This producer is a fixture only; the auditor never imports production numerical/source-parsing helpers.

Coverage includes ordinary inference; zero B; one/both constant arms; exact zero Delta with positive variance; extremely small nonzero variance; fourth-order cancelling Delta; positive probability below native tail range; native numeric unavailability; all eight leaveouts; invalid denominators; malformed/hash/axis/coordinate gates; authorization, namespace and fresh-output guards; stale/partial inventories; payload tampering; and private nested leaks. The default is tested with matrix/archive opens forbidden.

Final aggregate test receipts, logs and exact code/test/design/protocol/dependency pins are stored under `work/psc-bach2-naive-rna-execution-review/`. They identify a stable handoff, not real-data execution authority. No real BACH2 values, selected-donor RNA/ADT quantities, real effects, source-value runs, network, installs, Git actions or edits to others' owned files are part of this task. The prior source reviewer disclosed opening a closed cohort qualification JSON without extracting/using its technical totals. This adapter does not open that JSON, and no such totals enter the reference.

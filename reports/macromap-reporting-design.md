# MacroMap post-fit aggregate reporting design

This document describes presentation only. It does not amend the frozen science.
The root owns interpretation, publication authorization, and reproduction/navigation documents.
The frozen methods, source preparation, primary results, and original review are read-only.

## Inputs and authorization

`python scripts/report_macromap_program_context.py` validates aggregate inputs by default.
`--publish` is a separate, root-controlled action. It requires fresh destinations.
No publication output from the real primary is generated during reporter development.

Required external pins authenticate:

1. The frozen plan and the primary `execution-complete.json`.
2. The private PASSED independent receipt. Its exact successful status is
   `independently_verified_with_declared_source_subset`.
3. The root-inspected, byte-identical public verifier receipt at
   `reports/macromap-independent-verification.json`.
4. The root's public mechanical replay receipt at
   `reports/macromap-replay-verification.json`.
5. The root-approved `reports/macromap-independent-optimizer-diagnostics.json`
   aggregate diagnostic, including execution/model commitments checked against
   the sealed inventories without opening model arrays.
6. The prepared source-only summary and every aggregate JSON/CSV actually read,
   against their already sealed inventories.

The reporter does **not** rehash/read raw or normalized expression, model NPZ,
probabilities, individual paired scores, individual held-out losses, or source
sample/line tables. The complete private inventory was authenticated by the
frozen independent verifier and root's replay check. This reporter authenticates
only its approved aggregate payloads. It does not rerun a model, score, bootstrap,
statistical test, source preparation, or numerical verification procedure.

The root inspected the private verifier receipt before publishing its public
copy. The reporter nevertheless constructs a smaller explicit scalar-only
verification summary. It does not copy nested fit checks, selection hashes,
source checks or other arbitrary receipt objects. It retains the actual maximum
normalization-subset error, including a small **nonzero** value; it never replaces
that value with a claim of exact zero error. Mechanical replay is not biological
replication. Its only differing listed artifact is expected execution-start
metadata; all other listed artifact hashes and all endpoints agree.

## Fixed export registry

All eleven tables are deterministic, UTF-8, gzip-compressed CSV. Gzip mtime is zero
and the gzip header has no original filename. Empty cells mean null/unavailable,
not zero. Each table has an explicit column whitelist in `COLUMNS`.

| File | Rows | Scope |
| --- | ---: | --- |
| `responses.csv.gz` | 3,024 | Both identity scopes; both protocols and pooled descriptive rows; all programs, times, conditions and quantities |
| `predictions.csv.gz` | 20 | Both scopes × both fixed variants × four protocol/time cells and pooled endpoint |
| `leave-run-summary.csv.gz` | 12 | Both scopes/variants; omissions grouped by protocol and an all-omission summary |
| `model-folds.csv.gz` | 80 | All fixed outer folds and statuses, including failed states |
| `model-fits.csv.gz` | 160 | Both components of every paired fit, including unattempted/failed components |
| `matching-status.csv.gz` | 360 | Nine programs in each of 40 scoring contexts |
| `scale-status.csv.gz` | 80 | Approved independently checked constant-feature counts/names; no means/SDs |
| `independent-optimizer-qc.csv.gz` | 160 | Separate raw rerun solver flags, numeric criteria and aggregate diagnostics |
| `mapping.csv.gz` | 9 | Source mapping gates, parent coverage and fixed overlap-removal counts |
| `paired-design.csv.gz` | 80 | Both scope-specific paired-response registries, including PIC |
| `cohort-design.csv.gz` | 4 | Source/mapped and fixed predictive counts by scope/protocol |

Each response scope has 1,458 quantitative rows plus **54 PIC stubs**. PIC has no
qualified mock identity. Its source row has no quantity or effect estimate. The
export keeps this exact key convention, not three invented zero-valued rows.

Primary predictive N is **185 = 114 + 71**. Identity-inclusion sensitivity N is
**189 = 114 + 75**. Both times receive weight one half, exactly once. The protocol
weights remain their original line shares. The nonoverlap program and identity
inclusion are fixed sensitivities, not independent cohorts or a best-result
selection. Only mapped-primary/full-program/pooled gain is labeled the primary
endpoint. Protocol/time cells remain labeled cells, not extra primary endpoints.

Response cohorts can exceed the both-time prediction cohort. Their original
counts and contrast-specific protocol shares are validated independently. They
are not capped at 114/71 or replaced with predictive weights.

## Missingness, uncertainty and consistency checks

- Validate complete Cartesian keys and reject duplicate, missing or extra cells.
- Retain original known status codes. Unknown states stop reporting for review.
- Reject nonfinite aggregates, bool-as-count, malformed numeric values, and
  arbitrary extra response/endpoint metadata. Unavailable effect values remain
  null; a numeric zero is accepted only for an actually available result.
- Retain input-score status counts even when raw responses remain available but
  matched background or adjusted responses fail.
- Reject **all pooled response CIs**. Both pooled-response mean and median remain
  descriptive with explicit `not_computed_prespecified` interval statuses.
- Preserve protocol-specific pointwise, fixed-score run-bootstrap intervals.
- Preserve all available predictive CIs, including the mandatory primary pooled
  CI. These are **conditional fixed-OOF** intervals, not full-procedure or
  unconditional algorithm intervals. Known unavailable-uncertainty source
  states remain explicit; no interval is fabricated.
- Check original scope-specific n/weights, paired loss differences, and pooled
  means against the four original cell values. Check held-out line-weighted fold
  aggregates against those cell values. This is arithmetic validation of saved
  summaries, not model recomputation.
- Check manifest, fit audit, fold identities/statuses and original counts.
- Leave-run summaries keep counts, signs, gain ranges and failure-status counts.
  They never emit omitted-run IDs. Their gain always refers to the **pooled**
  predictive endpoint after an omission, not a new protocol-specific effect.
- Read mapping count strings with exact Decimal validation, including integral
  decimal forms such as `926.0`. Keep blank derived/parent fields null.

## Privacy and scaling limits

No individual line/donor/sample/run identity, index vector, source gene list,
matching realization, fitted coefficient, scaler vector or probability array is
published. Free-text optimizer/reference error messages are never copied.
Public tables are constructed from fixed labels, known statuses and validated
scalars, not recursive serialization or a broad denylist.

The original fit-audit JSON does not contain scaler values or constant-feature
flags. The root supplied the separate approved aggregate optimizer diagnostic
with SHA-256 `586168f8e9c8cf6e246aba7e235445297e22cf0d57f7536d2207bfa9028b1232`.
The reporter validates all 160 model keys and 80 model commitments against the
manifest/completion metadata. It does **not** open or hash any NPZ file.

Constant-feature names are restricted to the actual fixed model predictors:
four comparators for baseline and those four plus the nominated variant for
extended. The 80 scale rows use the extended five-feature summary and confirm
that baseline flags equal the corresponding four-feature subset. The approved
diagnostic independently reconstructed saved constant flags at the frozen
training SD threshold of `1e-12`. Means and scales remain unreported.

The diagnostic contains **120 true and 40 false raw solver-success flags**, with
all **160 frozen numerical criteria passing** and no constant-feature models.
False flags are preserved, not hidden or converted into numerical failures.
These are **reruns of the identical frozen independent helper**, not stored
original SciPy result objects. The helper did not return message/status codes;
no failure reason is inferred or invented. Its unchanged gradient limit is
`1e-7`, and coefficient rtol/atol remain `2e-5`. This post-fit capture neither
replaces the original verification nor changes a primary fit, optimizer,
threshold, endpoint or interpretation rule.

## Figures

Both PNG and SVG are generated with explicit deterministic metadata and a fixed
SVG hash salt. The runtime manifest records Python, NumPy, Matplotlib, Pillow and
zlib versions. Replay determinism is conditional on the same native rendering
runtime and fonts.

- **Prediction:** four protocol/time cells and the pooled endpoint, with both
  fixed variants. Identity-inclusion sensitivity is a separate panel, explicitly
  not an independent cohort. Whiskers show conditional fixed-OOF 95% CIs only.
  Failed results and unavailable intervals stay visible. Gain is baseline minus
  extended mean log loss; positive values mean lower held-out log loss.
- **Responses:** mapped-primary median **adjusted** responses for all nine
  programs × all 20 registered conditions, separately for both protocols.
  The two panels share a symmetric diverging color range. PIC and any unavailable
  effect are grey and labeled N/A, not rendered as zero or silently removed.
  Scores have arbitrary relative units. There is no response significance screen.

Complete response data for all quantities and both identity scopes remain in
the table, irrespective of the display subset. Neither figure establishes PSC
validation, ETS2 activity/dependence, mechanistic exclusivity, or clinical benefit.

## Destinations, manifests and replay

The intended first public destinations are:

- `data/derived/macromap-program-context/`
- `reports/figures/macromap-program-context-prediction.{png,svg}`
- `reports/figures/macromap-program-context-responses.{png,svg}`

Destination parents must already exist. The table directory and all four figure
paths must be fresh. Source/destination symlinks, including dangling links and
literal parent traversal, are rejected. The reporter never replaces an existing
file. It rejects overlaps with immutable bundles and raw inputs. Input hashes
are rechecked before writing and before the success manifest. A failed partial
publication is not success: it has no success manifest, must not be committed,
and requires a new destination. This is not a hostile concurrent-filesystem
sandbox or a transaction spanning all destination files.

`reporting-manifest.json` records input pins, each authenticated aggregate
payload, runtime, explicit table columns/row counts, output hashes and limitations.
It excludes destination-specific paths and wall-clock timestamps. A clean
reporting replay uses new table and figure-prefix paths. All table, SVG, PNG and
manifest bytes must match. The primary scientific artifacts are never rewritten.

## Validation boundary

All 25 native synthetic reporter tests pass. They cover full registry
projection, scalar whitelists, unavailable/zero distinction, status-preserving
failures, fixed weights, wrong hashes, malformed schemas, forbidden-read paths,
exact count strings, symlinks/parent traversal, fresh destinations, input changes,
byte-identical SVG/PNG/gzip/manifest replay, approved-QC hash/identity checks,
false solver flags versus numeric pass, four/five-feature mask consistency and
constant-feature-name privacy. Independent read-only schema
review identified presentation defects; repairs apply only to this new reporter.

Validation-only reading of authenticated real aggregate inputs also passes.
It reports the expected table cardinalities, including the optimizer-QC table,
without creating publication data or figures. The core authenticated reporting
path passed bounded independent static privacy/schema review. The optimizer
extension also passed its separate narrow metadata/projection review. The
combined native MacroMap suite passes **96 tests**, including the 25 reporter
tests. Synthetic-only prediction/response figures were visually inspected;
no primary biological figures were generated by this worker. These are
presentation checks, not new biological tests or numerical re-verification of
the frozen science.

The final CLI requires `--plan`, `--plan-sha256`,
`--execution-complete-sha256`, `--verification-receipt`,
`--verification-receipt-sha256`, `--public-verification-summary`,
`--public-verification-summary-sha256`, `--replay-verification`,
`--replay-verification-sha256`, `--optimizer-diagnostics`,
`--optimizer-diagnostics-sha256`, `--data-output`, and `--figure-prefix`.
Add `--publish` only for root-authorized final generation. Root handoff records
supply the exact pinned command and fresh replay destinations. No real primary
publication was generated by this worker.

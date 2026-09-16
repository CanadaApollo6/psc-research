# PSC liver: bounded aggregate-reporting design

## Scope and approval boundary

This is a reporting-only adapter for the already completed, independently verified
GSE303271 analysis. It does not import or execute the scientific runner, load raw
counts or private outputs, recalculate scores/statistics/corrections, or change the
frozen protocol. It makes no new scientific choice. Publication needs a separate
root acknowledgment. The reporter worker creates only ignored previews; root
reviews the actual figure, authenticates reporting replay, writes the biological
narrative, and decides whether to publish.

The sole allowed scientific input is the exact `aggregate-public-summary.json`
named in `work/psc-liver-program-root-freeze-v1/root-reporting-handoff.json`.
Its SHA-256 is `350ee5d1e447b1cab4621595611de99f7a787af40102187ca0377257c9a598d0`.
The reporter also reads only that handoff, the frozen plan, outer seal, named
independent audit/replay receipts, protocol and verification design, plus its own
new script/test/design and installed rendering metadata. It does not follow any
private source path mentioned by a receipt. All eight input files have literal
expected paths, sizes and SHA-256 pins in the reporter. Hashes are checked before
JSON parsing, again after rendering and before completion.

The external root freeze is SHA-256
`5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d`.
The outer seal is SHA-256
`59eb70cb6a539d84799300332c1ece54bf5384e9655d26884dacf06b62df7c53`.
The native completion SHA-256 is
`b3b4e5d35dc48c2fa52604dadc0f45dd6c10f2c7a45567e924935140d49ce524`.
The checkpoint is `c1ce7f25934c3797eed4660cfed414749c8b46bf`.
Native `root_freeze_written=false` means the native writer did not create this
external root freeze. It does not mean that no freeze exists.

## Schema, privacy and completeness

Validation precedes every projection or result-directory creation. JSON duplicate
keys and nonfinite literals are refused. Every aggregate object has exact keys;
every nested list, string, null, boolean and count has an explicit typed contract.
Boolean values cannot masquerade as integer counts. All variable numerical leaves
must be finite with valid domains. Nested Welch/bootstrap intervals and omission
ranges have order checks. Exact status-dependent availability is enforced.

This adapter deliberately accepts **one verified snapshot**, in which all 60
Welch and bootstrap results are available. An unobserved unavailable-status branch
fails closed and requires a separately reviewed input contract; it is never
converted to a valid P value, a zero, or a missing row. Existing null fields remain
JSON null or empty CSV cells. Empty reason/variance-group lists remain `[]`.
Strict q values, family/correction fields, and optional targets cannot be replaced
with zero or one. The private numerical `p_for_correction` field is checked but
not exported. No correction or inferential calculation occurs in reporting.

The fixed source order is Core then Strict, ten programs in protocol order, then
PSC minus AIH, PSC minus ASC, ASC minus AIH. Exactly 60 unique ordered rows are
required. Core retains all 30 hypotheses with both BH and BY; Strict retains all
30 rows as identity sensitivity only, with no discovery q values. Duplicate,
missing, added or reordered rows are refused. Reused group summaries must match
exactly across contrasts. Coverage preserves the nonoverlap parent gate: its
component coverage flag is false but its parent gate is true. ECM uses 321 official
symbols and mapped counts 301/299, not 340 Ensembl records.

All exported audit checks, tolerance fields and runtime strings have a recursively
typed fixed specification. Receipts are not copied wholesale. The manifest and
completion have closed schemas and independently checked actual file inventories.
The root handoff reports zero actual source-identity token hits in this aggregate.
The reporter does not open those IDs to repeat that check. Exact byte authentication
and recursive leaf restrictions are a **bounded, hash-specific privacy check**,
not a universal guarantee or a security sandbox against arbitrary local code or
filesystem control.

## Complete public projections

A preview mirrors these later public logical paths:

| File | Contents |
|---|---|
| `data/derived/psc-liver-programs/endpoints.csv` | All 60 effects, raw/log P, Core BH/BY, pointwise Welch CIs, n/means/variances, SE/df/t, statuses, family flags and complete flattened coverage |
| `data/derived/psc-liver-programs/group-summaries.csv` | All 60 tier/program/group n, means and sample variances (`ddof=1`); no SD recalculation |
| `data/derived/psc-liver-programs/bootstrap-intervals.csv` | All 60 pointwise percentile intervals, method, requested/valid counts, statuses and no-P flag |
| `data/derived/psc-liver-programs/omission-stability.csv` | All 60 public count/range/change/sign summaries; only 47/34/47 relevant-unit denominators; no individual omissions |
| `data/derived/psc-liver-programs/mapping-summary.csv` | All 20 fixed tier/program mappings, original/member denominators, parent/component gates and target metadata |
| `data/derived/psc-liver-programs/held-endpoints.json` | IL17A response-transfer HOLD, no score/test/empty surrogate; GSE159676 HOLD and NoPSC atlas PARK |
| `data/derived/psc-liver-programs/verification-provenance.json` | Explicit aggregate verification counts, tolerances, input/code pins and interpretation/privacy limits |
| `reports/figures/psc-liver-programs.png` and `.svg` | Three contrast panels, all ten programs, separate Core/Strict pointwise Welch CIs and a zero reference |
| `data/derived/psc-liver-programs/report-manifest.json` | Nine payload hashes/sizes, row counts, runtime/font pins and logical destinations |
| `data/derived/psc-liver-programs/completion-receipt.json` | Ten preceding file hashes/sizes; written atomically last; not self-listed |

The endpoint table keeps group variances as published rather than computing SDs.
CSV uses UTF-8, LF, fixed columns, round-trip numeric serialization and no rounding.
The five table row counts are 60/60/60/60/20. All 11 final files are checked against
the actual write inventory. No whole native private output is published.

The provenance distinguishes the independently verified 2,630,144 raw/cache cells,
64 totals, 1,280 scores, 60 endpoints, 30-test Core correction, 640,000 shared RNG
indices, 600,000 bootstrap differences and 3,840 private omissions from this
reporter's aggregate-only work. No saved native log-CPM comparison is claimed.
The prior native replay checked 24 actual hashes and all 12 native files including
completion. Root's 73-file seal-preservation attestation is labeled as such; this
reporter rehashes only its allowed eight inputs, not forbidden native sources.

## Figure and interpretation

The plot shows mean `log2(1 + CPM)` composite differences and **pointwise,
unadjusted** Welch 95% intervals. It is not a raw fold-change plot, a simultaneous
confidence display, an individual-value plot, or a TF-activity result. Strict is
visibly separate from the Core discovery family. All Core corrections remain in
the complete table. Bootstrap intervals and omission summaries are complete tables,
not extra selected figures. Omission ranges are not confidence intervals.

All ten scores use positive equal mapped-member weights; gRNA1-up is not inverted.
A positive unsigned contrast does not imply activation or induction. ECM is
unsigned constitutive organization, not a fibrosis-response score. The focal
comparison remains ETS2 gRNA1-down, PSC minus AIH; reporting does not select a new
endpoint. Whole-biopsy composition, clinical and technical confounding, no healthy
group, and article-assumed rather than independently established patient-sample
independence remain explicit. No equivalence, generic-inflammation, PSC causality,
cell-intrinsic regulation, clinical benefit or independent MacroMap validation is
claimed. The held IL17 response-transfer question remains unanswered.

## Determinism, write safety and native commands

Use `.venv/bin/python` and the installed Matplotlib only. Rendering uses Agg,
explicit packaged DejaVu Sans regular/bold files, fixed SVG hash salt, fixed size,
DPI and axis range, and no varying date metadata. Matplotlib, NumPy (rendering
library dependency), Pillow, FreeType, Python, platform and font hashes are recorded.
The only runtime cache is under ignored `work/psc-liver-program-reporting/.matplotlib`.
Byte determinism is established by a fresh render in the recorded runtime; it is
not promised across arbitrary library, font or platform changes.

All input and output paths must have canonical spelling. Traversal, aliases,
symlinks, hardlinked input files, nonregular inputs, reused destinations, existing
figure targets, source collisions, extra files/directories and incomplete bundles
are refused. Only fresh named children of `work/psc-liver-program-reporting/` are
valid previews. No validation refusal creates a result directory. Rendering is
completed in memory before a preview is staged. Actual files are checked before
an atomic no-replace completion installation. Owned partial files are cleaned on
write failure. These are workflow safeguards, not concurrent-attacker isolation.

Publication requires an explicit root acknowledgment, a complete reviewed preview,
its externally supplied completion SHA-256, matching input/code/runtime pins and
full inventory/projection validation. It copies the exact preview bytes; it does
not rerender the figure during publication. Root must first visually inspect the
actual PNG and authenticate the separate fresh reporting replay. An unrelated
completion hash is not a review approval. Publishing to an existing destination is
never a repair path. Publication across the new data directory and two figure files
is not a single filesystem transaction; completion is the last commit marker.

The exact code/test/design pins and fully substituted commands are supplied in
ignored `work/psc-liver-program-reporting/reporting-handoff.json` after validation.
The command forms are:

```text
.venv/bin/python -m unittest discover -s tests -p test_report_psc_liver_programs.py -v
.venv/bin/python scripts/report_psc_liver_programs.py validate --expected-reporter-sha256 REPORTER_SHA256
.venv/bin/python scripts/report_psc_liver_programs.py preview --expected-reporter-sha256 REPORTER_SHA256 --output-directory work/psc-liver-program-reporting/preview-v1
.venv/bin/python scripts/report_psc_liver_programs.py replay --expected-reporter-sha256 REPORTER_SHA256 --preview-directory work/psc-liver-program-reporting/preview-v1 --expected-preview-completion-sha256 PREVIEW_COMPLETION_SHA256 --output-directory work/psc-liver-program-reporting/replay-v1
.venv/bin/python scripts/report_psc_liver_programs.py publish --expected-reporter-sha256 REPORTER_SHA256 --preview-directory work/psc-liver-program-reporting/preview-v1 --expected-preview-completion-sha256 PREVIEW_COMPLETION_SHA256 --root-publication-ack ROOT_APPROVES_PSC_LIVER_AGGREGATE_PUBLICATION
```

The tests use only this authenticated aggregate and metadata, with mutations or
sandbox copies under the reporter's ignored work directory. They do not import the
scientific implementation or execute native analysis. They cover grids, source
hashes, every nested schema class, null/zero and Strict-family rules, all omission
denominators, duplicate/nonfinite JSON, unsafe paths, partial inventories, write
failures, exact projections, deterministic figure replay and sandbox-only publication.
No package installation, network access or scientific rerun is required.

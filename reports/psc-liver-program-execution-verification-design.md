# Pre-effect design: independent native liver-output verification

## Status and limits

This is a **post-fit agreement audit design**, prepared before any real RNA values
were opened for this task. It is not a new endpoint, alternative analysis,
execution approval, or biological result. A separate explicit root permission is
still required for both the native analysis and the later real audit.

The existing accepted author implementation, tests, independent math, math tests,
and method/review reports remain byte-for-byte unchanged. The new adapter never
imports author helpers. It may use only the SHA-authenticated independent math
primitives in `scripts/verify_psc_liver_programs.py`. It never calls
`verify_numeric_bundle`, changes that function's synthetic-only condition, or
labels real data synthetic.

Source mapping is reused by **external authentication**, not accepted on the
native output's assertion. The independently accepted review receipt binds the
canonical compiled-program SHA-256
`17a291f47fbb0ec0e815ac2f592ba4d0f657f3f5e3cbe4689dae5ecc44c90698`.
The private runtime compilation must match it exactly. The complete ECM identifier
audit must also match the already reviewed source-only artifact. Original source,
qualifier/preparation code, test, handoff, manifest, membership and cached official
ECM page/JSON pins are checked. No source is fetched again.

## Fixed agreement scope

Under future root authorization, the adapter checks:

1. All **41,096 × 64 = 2,630,144** original gzip count cells against the frozen
   NumPy cache. It parses UTF-8 tab-delimited CSV independently. The first nonblank
   row must have an empty first cell and the complete literal column axis. Every
   following row must match the closed original feature-label axis. Blank lines
   are formatting, not features. Gzip EOF/CRC must complete. Counts must be exact,
   finite, nonnegative decimal integers no greater than `2**53`. Exact full
   libraries must be positive and no greater than `2**53`. There is no rounding,
   count repair, source-row deletion, alias rescue or zero-fill.
2. All **64 exact full-source library totals**, all **20** tier/program score
   vectors (**1,280 values**), and all **60** ordered endpoint/status rows.
3. The complete **30-test Core** raw/correction/BH/BY vector, including exact
   missing-P placeholders and public null fields. Strict has no discovery q values.
4. All **640,000 shared PCG64 indexes**, all **600,000 bootstrap differences**,
   and every bootstrap interval/status. The fixed seed is `2026091602`, B is
   `10000`, group order is PSC/ASC/AIH, and within-group order is the original
   source-column order. An unavailable score vector is never replaced with zeros.
5. All **3,840 private omission rows**, including every Welch field and status.
   Outside-contrast statistics must be an exact copy of the saved native baseline,
   with exact delta zero and no sign reversal. All 60 public stability summaries
   use only the relevant **47/34/47** units.

The source-approved production compilation makes all 20 score vectors available.
Legitimate zero-variance or numerical-inference-unavailable endpoints retain the
fixed rows and missing-value policy. These are checked, not reclassified as
biological nulls. Artificial fixtures separately exercise unavailable score
vectors and missing native bootstrap keys.

**No native logCPM matrix exists.** The adapter computes its own reference
normalization to derive scores, but it does not fabricate a native normalization
matrix comparison. Its receipt explicitly says that comparison is unavailable
and is not claimed. Source/cache equality is not a replay of sequencing,
alignment, historical GRCh38 counting-GTF assignment or original clinical labels.

## Agreement criteria fixed before values

- Counts, totals, source axes, membership digest, schema, row inventory, statuses,
  RNG indexes and fixed plan fields are **exact**.
- Positive P/Q values, variances, standard errors and df use relative tolerance
  **`1e-10` with no absolute floor**. This also covers private correction vectors.
- Other finite floating fields use relative **`1e-10`**, absolute **`1e-12`**.
- Domain checks apply before tolerance: available P/Q must be in `(0,1]`,
  positive variances/SE/df cannot be zero, exact-zero variance must be zero,
  unavailable inferential fields must be null, unsigned scores cannot be negative,
  intervals must be ordered, and a degenerate interval must be exactly degenerate.
- Missing Core correction inputs and their private adjusted placeholders are
  exactly one; public raw P/q remain null. No observed P=0 is accepted.

The independent review exposed draft-adapter gaps in global float tolerance,
fixed-plan equality, missing source pins, numerical-failure translation and exact
outside-contrast copies. These were repaired using artificial fixtures only.
`pre-effect-agreement-criteria.json` records the tighter criteria. The frozen math
was not edited. Its unrepresentable-SD errors are translated only into explicit
unavailable-variance states with finite descriptive effects; no epsilon variance
or substitute valid inference is introduced.

## Authentication, inventory and privacy

The real entrypoint requires externally supplied hashes for the adapter, original
root-frozen plan and native completion receipt. It authenticates the accepted
review, author/code/tests, source lineage and Python/NumPy/SciPy versions before
numerical work. The plan must use the original schema; only state, root
authorization, independent-review acceptance and expected compiled digest may
change from the prepared template. Other plan fields are exact. Accepted plan
paths are:

- `config/psc-liver-program-plan.json`
- `config/psc-liver-program-freeze.json`
- `work/psc-liver-program-plan/proposed-plan.json`

A run must be an owner-private dedicated child of
`work/psc-liver-program-plan/runs/`. It must contain exactly the native **12 files**:
the completion receipt and its **11 listed artifacts**. The completion receipt
does not list itself. It must assert successful complete execution of the same
plan. Every actual artifact size and SHA must match the externally authenticated
completion receipt before numerical parsing. `RUNNING.json`, `FAILED.json`,
partials, extras, missing files, duplicate manifest entries, reused or unsafe
paths, symlinks and hardlinks are refused. NPZ member names, duplicate entries,
headers, shapes, dtype and bounded sizes are checked before array allocation;
pickle is disabled. JSON duplicate keys and nonfinite literals are refused.

The original native plan copy must match the frozen plan **byte for byte**.
Sources and the complete native inventory are reauthenticated after computation.
The adapter checks its loaded code hash again before sealing a receipt.

The output is a fresh ignored child directory under
`work/psc-liver-program-independent-review/execution-adapter/runs/`, mode `0700`,
with an atomically published `audit-receipt.json`, mode `0600`. Its fields are
aggregate verification counts, fixed criteria, software and artifact hashes.
There are **no gene labels, sample/person IDs, unit totals, scores, effects, P
arrays or omission values**. Errors use fixed codes, not values or tracebacks.
Failure does not issue a success receipt. These are workflow safeguards, not
cryptographic proof of who supplied a root marker or protection against arbitrary
code/filesystem control.

## Native commands

Metadata-only default, permitted now:

```text
.venv/bin/python scripts/verify_psc_liver_program_execution.py
.venv/bin/python -m unittest discover -s tests -p test_verify_psc_liver_program_execution.py -v
```

Future real audit, **not run by this preparation**:

```text
.venv/bin/python scripts/verify_psc_liver_program_execution.py --verify-real \
  --expected-adapter-sha256 ADAPTER_SHA256 \
  --plan config/psc-liver-program-freeze.json \
  --expected-plan-sha256 ROOT_FROZEN_PLAN_SHA256 \
  --execution-receipt work/psc-liver-program-plan/runs/RUN_ID/completion-receipt.json \
  --expected-execution-receipt-sha256 NATIVE_COMPLETION_SHA256 \
  --root-real-audit-ack ROOT_AUTHORIZES_COMPLETE_GSE303271_NATIVE_OUTPUT_AUDIT \
  --output-directory work/psc-liver-program-independent-review/execution-adapter/runs/AUDIT_ID
```

The real CLI cannot accept an artificial contract. Tests inject an explicitly
artificial contract into a private verification seam. They invoke the immutable
native writer only inside temporary roots, with generated counts and fake
identities. This does not certify positive real authorization. One full artificial
fixture has **41,096 rows and 64 columns**, including positive unmapped filler
rows, and all approved-size resampling/omission outputs. Smaller fixtures cover
unavailable and tampered states. A repeat verifies identical aggregate receipts.

Final implementation/test/contract/report pins, the focused test log and the
independent child's review are recorded in
`work/psc-liver-program-independent-review/execution-adapter/pre-effect-acceptance.json`.
No real run is certified by that pre-effect receipt.

## Interpretation remains unchanged

Ten unsigned programs and three contrasts remain fixed. The focal ETS2 gRNA1-down
PSC-minus-AIH comparison is primary. The epithelial IL17A-response question stays
HOLD outside this panel. GSE159676 stays HOLD and the NoPSC atlas stays PARK.
Whole-biopsy composition, relative-expression scaling, unknown clinical/stage/
batch covariates and paper-asserted—not independently verified—patient-sample
independence remain limitations. Welch and bootstrap intervals are pointwise,
not simultaneous or causal. Agreement cannot establish TF activity, fibrosis
stage, clinical benefit or independent mechanistic confirmation of MacroMap.

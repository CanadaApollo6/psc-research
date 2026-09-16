# Reproduce the fixed PSC liver program analysis

## Scope: two different reproduction tasks

The [completed report](../reports/psc-liver-programs.md) answers one closed,
observational question: how ten fixed unsigned transcript-abundance composites
compare among deposited PSC, ASC and AIH diagnostic liver biopsies. The focal
comparison is ETS2 gRNA1-down, PSC minus AIH. This guide distinguishes:

1. **Native scientific reproduction:** restore the exact frozen source, method,
   runtime and acceptance bundle; run the authorized primary, independent audit
   and mechanical replay into fresh ignored destinations.
2. **Aggregate presentation reproduction:** authenticate the already verified
   aggregate and receipts; reproduce the complete tables and figure without
   accessing private unit-level inputs or recalculating statistics.

A clone alone is **not a complete execution bundle**. Public reports, hashes and
configuration files do not supply missing raw counts, ignored caches, source joins,
acceptance bodies or runtime files. Missing or mismatched source, method or runtime
artifacts block reproduction. Do not edit a pin, substitute a newly fetched source,
regenerate an acceptance receipt, or silently rehash an environment/path difference
to make a check pass. Record the blocker. An environment adaptation would need a
separately authorized, explicitly different reproduction record.

Here **root means the coordinating research agent**, not the operating system
superuser. No `sudo` or elevated operating-system account is required or implied.
The acknowledgment strings below are explicit workflow gates, not access grants
or proof that an arbitrary caller has coordinator approval.

Read the [evidence log](evidence-log.md), [frozen protocol](psc-liver-program-protocol.md),
[source qualification](../reports/psc-liver-input-qualification.md),
[identifier qualification](../reports/psc-liver-program-inputs.md), and
[external-program qualification](../reports/psc-liver-external-program-qualification.md).
The [method design](../reports/psc-liver-program-method-design.md),
[independent verification design](../reports/psc-liver-program-execution-verification-design.md)
and [reporting design](../reports/psc-liver-program-reporting-design.md) define the
separate computational and presentation checks.

## 1. Restore and authenticate the native bundle

### Frozen public artifacts and ignored metadata

The pre-effect checkpoint was `c1ce7f25934c3797eed4660cfed414749c8b46bf`.
Use the exact immutable scientific files bound by the
[execution seal](../config/psc-liver-program-execution-seal.json), not a revised
method selected after results.

| Artifact | SHA-256 |
|---|---|
| [Frozen plan](../config/psc-liver-program-freeze.json) | `5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d` |
| [Execution seal](../config/psc-liver-program-execution-seal.json) | `59eb70cb6a539d84799300332c1ece54bf5384e9655d26884dacf06b62df7c53` |
| [Protocol](psc-liver-program-protocol.md) | `a36faad8bffee1e5f4f0789c8ce77827e1378477751c7e9e29ecd6bb1d3dcf90` |
| [Public verification contract](../config/psc-liver-program-verification-contract.json) | `af87402ee46bee56751f8acbc54bdf932d92db0924b434bf588666868bb53785` |
| Native runner, `scripts/analyze_psc_liver_programs.py` | `f7d5e49d2309c2e21a2ea20d0a7d1b8bda8643e140b0c781db3edafd34fd5ec7` |
| Native-output adapter, `scripts/verify_psc_liver_program_execution.py` | `34115ac11cf8fb9932f80846afb89d8cfcd04f2543476f1872a765aeeb5c55e1` |
| Independent math, `scripts/verify_psc_liver_programs.py` | `82403f2cb63ae1eeb7e4f4a1232b21ae6772515d0e1b3ead1ca7397737a62391` |
| Ignored `work/psc-liver-program-root-freeze-v1/private-seal-ledger.json` | `006ecb2e66020aea9cd6bdb80d0224e415cda43865a68cd6f24aa2b0e66a9eeb` |

The ignored ledger binds **73 files**. Restore its actual original body and every
listed body at its expected path, including the scientific tests, preparation
receipts, source-only preflight, identifier/membership artifacts, original review
acceptances, runtime capture artifacts and logs. Its pin is not a replacement for
those files. Authenticate all seals before and after native work. The accepted
review and preparation receipt are distinct from the later native run completion.

Important ignored restoration paths include:

- `work/psc-liver-program-independent-review/accepted-review.json`;
- `work/psc-liver-program-independent-review/execution-adapter/pre-effect-acceptance.json`,
  `pre-effect-agreement-criteria.json`, `review/candidate-review-acceptance.json`
  and the other adapter files listed in the ledger;
- `work/psc-liver-program-plan/completion-receipt.json`,
  `root-method-acceptance.json` and the complete original `root-source-preflight-v1/`
  bundle;
- `work/psc-liver-program-root-review-v1/root-pre-effect-acceptance.json`,
  `native-runtime.json`, `native-runtime.log` and the other sealed runtime/review
  artifacts;
- the complete qualified source metadata, membership, feature-index and source-join
  bodies specified by the contract's `source_metadata_pins` and the ledger.

The adapter requires
`work/psc-liver-program-independent-review/execution-adapter/production-contract.json`.
Restore this as an **exact byte copy** of the public verification contract linked
above, with SHA-256 `af87402ee46bee56751f8acbc54bdf932d92db0924b434bf588666868bb53785`.
This public copy permits recovery of that one contract, not fabrication of the
other ignored acceptance or source files. Do not substitute an artificial fixture
contract or relabel real inputs as synthetic.

### Original expression and source metadata

Two separately pinned expression payloads must be restored exactly. A saved cache
alone cannot support the independent raw-source audit.

| Original path | Bytes | SHA-256 |
|---|---:|---|
| `data/raw/psc-liver-qualification/GSE303271_raw_counts.txt.gz` | 3,021,812 | `966be15cf8cd893d0da3ec62d73ee0a56a34932a9d1c990c3d0d58b903a24046` |
| `work/psc-liver-qualification/GSE303271-submitted-float64.npy` | 21,041,280 | `4a19908237dfaeb2f6c713fd6c7f60743d9b31b556a643c9c0dfcab7159003d7` |

Restore the exact submitted feature order and source-sample join at
`work/psc-liver-qualification/GSE303271-features-as-submitted.tsv` and
`work/psc-liver-qualification/GSE303271-sample-source-join.json`, plus the complete
qualified feature audit, program coverage and mapped indices under
`work/psc-liver-program-inputs/`. Restore the original legacy membership exports,
Ensembl 116 reference, official ECM source caches, membership CSV and
`work/psc-liver-external-programs/method-schema-handoff.json` named by the contract.
Do not substitute an updated annotation, alias mapping, gene list or source export.
No source-patient identifiers or individual values belong in the public guide.

The required compiled-program **canonical** digest is
`17a291f47fbb0ec0e815ac2f592ba4d0f657f3f5e3cbe4689dae5ecc44c90698`.
It is not the serialized compiled JSON's file hash. Both identities must satisfy
their own frozen checks.

### Native runtime and environment

Use the restored `.venv/bin/python`, Python **3.12.14**, NumPy **2.5.3** and
SciPy **1.18.1**. The recorded platform is
`Linux-7.2.3-arch1-3-x86_64-with-glibc2.44`.
The executable SHA-256 is
`8f9781a98200d9ecda7e00464e4c64b1327abae788ae8e6979d5c859311410c7`.
The seal also pins both package `RECORD` hashes and the verified installed-file
inventory digests: **929 NumPy files** and **1,433 SciPy files**. Matching version
strings alone does not recover those exact distributions or the original runtime
metadata bodies. Preserve the original path-sensitive metadata. Do not amend it to
hide a changed interpreter location, environment path or package build.

Set these environment values before starting the documented native commands:

| Variable | Value |
|---|---|
| `PYTHONHASHSEED` | `0` |
| `PYTHONDONTWRITEBYTECODE` | `1` |
| `OMP_NUM_THREADS` | `1` |
| `OPENBLAS_NUM_THREADS` | `1` |
| `MKL_NUM_THREADS` | `1` |
| `NUMEXPR_NUM_THREADS` | `1` |

The plan must remain `state=frozen`, with the existing exact authorization,
accepted review and method/source pins. The protocol uses all **41,096 released
rows** in each library total, **64 source patient-samples** (17 PSC, 17 ASC, 30 AIH),
positive equal mapped-member means of `log2(1 + CPM)`, ordinary Welch inference,
and the complete 30-hypothesis Core BH/BY family. Strict is a separate 30-row
identity sensitivity, not a second discovery family. Bootstrap uses the fixed
10,000 shared PCG64 draws and seed `2026091602`; omission reporting uses only
47/34/47 in-contrast denominators. This summary is not permission to change any
other protocol choice.

## 2. Archived native primary, audit and replay commands

The following are the actual command records from ignored
`work/psc-liver-program-root-freeze-v1/root-execution-resume.json`. They are not
commands run by this guide. They require coordinator authorization and a fully
restored, authenticated bundle. **Do not rerun them in the completed workspace.**
Their result directories already exist, and the historical `>` redirections would
also truncate existing logs. Use them only in a separately prepared reproduction
workspace where all stated output directories and log destinations are fresh.
Never delete or overwrite the original completed runs to make a destination fresh.

Primary:

```text
.venv/bin/python -B scripts/analyze_psc_liver_programs.py --execute --plan config/psc-liver-program-freeze.json --expected-plan-sha256 5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d --acknowledge-real-execution ROOT_AUTHORIZED_REAL_EXPRESSION_EXECUTION --output-directory work/psc-liver-program-plan/runs/primary-v1 > work/psc-liver-program-root-freeze-v1/primary-execution.log 2>&1
```

Independent audit, only after authenticating the actual primary completion and
its complete payload inventory:

```text
.venv/bin/python -B scripts/verify_psc_liver_program_execution.py --verify-real --expected-adapter-sha256 34115ac11cf8fb9932f80846afb89d8cfcd04f2543476f1872a765aeeb5c55e1 --plan config/psc-liver-program-freeze.json --expected-plan-sha256 5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d --execution-receipt work/psc-liver-program-plan/runs/primary-v1/completion-receipt.json --expected-execution-receipt-sha256 b3b4e5d35dc48c2fa52604dadc0f45dd6c10f2c7a45567e924935140d49ce524 --root-real-audit-ack ROOT_AUTHORIZES_COMPLETE_GSE303271_NATIVE_OUTPUT_AUDIT --output-directory work/psc-liver-program-independent-review/execution-adapter/runs/primary-v1 > work/psc-liver-program-root-freeze-v1/primary-independent-audit.log 2>&1
```

Mechanical replay, with the identical scientific plan and only a fresh result
destination (and separate log destination):

```text
.venv/bin/python -B scripts/analyze_psc_liver_programs.py --execute --plan config/psc-liver-program-freeze.json --expected-plan-sha256 5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d --acknowledge-real-execution ROOT_AUTHORIZED_REAL_EXPRESSION_EXECUTION --output-directory work/psc-liver-program-plan/runs/replay-v1 > work/psc-liver-program-root-freeze-v1/replay-execution.log 2>&1
```

The primary and replay native completion SHA-256 is
`b3b4e5d35dc48c2fa52604dadc0f45dd6c10f2c7a45567e924935140d49ce524`.
A different result must be reported as a mismatch, not approved by replacing the
expected hash. Native `root_freeze_written=false` only means that the native writer
did not create the external coordinator freeze.

### Native inventory and verification scope

Each completed native run has **11 payload files plus its completion receipt:
12 files in total**, not 11. The exact payload names are:

```text
aggregate-public-summary.json
private-all-unit-omissions.json
private-bootstrap-differences.npz
private-bootstrap-draws.npz
private-compiled-programs.json
private-correction-vector.json
private-external-identifier-audit.json
private-frozen-plan-original.json
private-library-totals.json
private-scores.json
private-source-units.json
```

`completion-receipt.json` lists those 11 payloads and does not list itself. Private
outputs stay ignored. Missing files, extra files, failed-run markers, unsafe paths
or a missing valid completion are not a completed analysis.

The accepted independent audit compared all **2,630,144 raw/cache count cells**,
64 full library totals, 1,280 scores, 60 endpoint rows, the complete Core 30-test
correction, 640,000 shared RNG indices, 600,000 bootstrap differences and all
3,840 private omissions. It checked all 60 public omission summaries. Exact
inventory/status/index checks and field-specific tolerances are in the frozen
verification design. It did **not** compare a saved native log-CPM matrix: that
intermediate was not saved, so such a native-matrix comparison is unavailable and
is **not claimed**.

The independent audit receipt SHA-256 is
`9a5827803aa44a9124e59bb491d9b9f6d021b1c5fdd79d82cfbed73c1b856004`.
The separate native mechanical replay compared **all 12 corresponding files,
including completion**, byte for byte: 24 actual file hashes, all identical.
It was not a new scientific analysis or biological replication.

## 3. Public aggregate inventory

The public presentation has **11 files: nine CSV/JSON data files and two figures**.
This differs from the native 12-file private-run inventory.

| Under `data/derived/psc-liver-programs/` | Scope |
|---|---|
| `endpoints.csv` | 60 rows: 30 Core and 30 Strict; source effects, raw/log P, Core BH/BY, pointwise Welch intervals, variance/SE/df/t, status and mapping coverage |
| `group-summaries.csv` | 60 tier/program/group rows; source n, mean and **sample variance (`ddof=1`)**, not newly calculated SD |
| `bootstrap-intervals.csv` | All 60 pointwise percentile interval/status rows; no bootstrap P values |
| `omission-stability.csv` | All 60 aggregate range/count/change/sign summaries, with 47/34/47 relevant-unit denominators only |
| `mapping-summary.csv` | All 20 tier/program mapping summaries, including the nonoverlap parent gate and ECM's 321-symbol denominator |
| `held-endpoints.json` | IL17 response-transfer HOLD, GSE159676 HOLD and NoPSC PARK; no zero substitute |
| `verification-provenance.json` | Explicit aggregate verification counts, criteria, source/reporting pins and limitations |
| `report-manifest.json` | Nine payload hashes/sizes, five table row counts, rendering/runtime/font pins |
| `completion-receipt.json` | Ten preceding files, including the manifest; does not self-list; installed last |

The two figure files are [PNG](../reports/figures/psc-liver-programs.png) and
[SVG](../reports/figures/psc-liver-programs.svg). They show all ten programs in
three contrast panels, with Core and Strict pointwise Welch intervals separately.
No individual scores or source-unit values are plotted. The
[public manifest](../data/derived/psc-liver-programs/report-manifest.json) and
[completion](../data/derived/psc-liver-programs/completion-receipt.json) authenticate
all published bytes. The completion SHA-256 is
`ef924593491e9c5b5b578e177a68f24de8f6af74c39206fe41594e60fed01211`;
the manifest SHA-256 is
`f02cb5180f99335746254d5ba18f717ea5fc8949b39e716966d28e823a45d6a4`.

## 4. Reproduce presentation without rerunning science

Restore the **eight exact inputs** listed in the public reporting manifest:

- the approved native `aggregate-public-summary.json`;
- `work/psc-liver-program-root-freeze-v1/root-reporting-handoff.json`;
- the frozen plan and execution seal;
- `work/psc-liver-program-independent-review/execution-adapter/runs/primary-v1/audit-receipt.json`;
- `work/psc-liver-program-root-freeze-v1/root-mechanical-replay-verification.json`;
- the frozen protocol and native-output verification design.

The input aggregate SHA-256 is
`350ee5d1e447b1cab4621595611de99f7a787af40102187ca0377257c9a598d0`.
Restore the exact reporter, reporter tests and reporting design pinned by the
manifest. They are later presentation artifacts, not part of the earlier native
checkpoint. The reporter SHA-256 is
`2c502da9107deacab445df7f2e81d0496e0e479c4875e8bf4509df11bc892fba`.
The tests SHA-256 is
`72942df8c321db4ad6217e6268104dcc953f1dd034c892bbcf01b5a3f1f9ea99`;
the reporting-design SHA-256 is
`31a2a187728384f0c637f7fca64b987978eab5c63b4d3c94611c8791f6e6569c`.
No raw counts, source-patient IDs, private scores/totals/omissions or arrays are
opened by this reporting route.

Rendering uses the recorded Python/NumPy versions, Matplotlib **3.11.2**, Pillow
**12.3.0**, FreeType **2.14.3**, Agg and the two hash-pinned packaged DejaVu Sans
fonts. Fixed SVG salt and omitted varying date metadata make repeat rendering
deterministic in that runtime, not across arbitrary library/font/platform changes.
The approved reporting tests passed **30 tests**. Coordinator visual review found
the actual PNG readable, and a separate fresh coordinator replay verified all
**11 presentation files**, including completion, with **22 actual hashes**.

Validation is read-only with respect to result artifacts:

```text
.venv/bin/python -B scripts/report_psc_liver_programs.py validate --expected-reporter-sha256 2c502da9107deacab445df7f2e81d0496e0e479c4875e8bf4509df11bc892fba
```

Fresh ignored preview and replay examples follow. The dedicated parent
`work/psc-liver-program-reporting/` must exist; these named children must not.
The names below are examples for a new reproduction, not existing run repair.
The expected completion stays fixed; an environment-induced difference is a
reproduction mismatch, not permission to replace the expected digest.

```text
.venv/bin/python -B scripts/report_psc_liver_programs.py preview --expected-reporter-sha256 2c502da9107deacab445df7f2e81d0496e0e479c4875e8bf4509df11bc892fba --output-directory work/psc-liver-program-reporting/reproduction-preview-v1
.venv/bin/python -B scripts/report_psc_liver_programs.py replay --expected-reporter-sha256 2c502da9107deacab445df7f2e81d0496e0e479c4875e8bf4509df11bc892fba --preview-directory work/psc-liver-program-reporting/reproduction-preview-v1 --expected-preview-completion-sha256 ef924593491e9c5b5b578e177a68f24de8f6af74c39206fe41594e60fed01211 --output-directory work/psc-liver-program-reporting/reproduction-replay-v1
```

The reporter validates the whole recursive schema before projection. It accepts
only the exact authenticated completed snapshot, preserves null/empty fields and
all rows, exports no correction placeholders, and recalculates no P values,
corrections, scores or model estimates. Its privacy check is bounded to those
specific bytes; it is not a universal privacy guarantee.

### Archived publication record: do not overwrite existing outputs

The public data directory and both figure targets already exist in the completed
repository. **They must not be overwritten.** Publication is not needed to inspect
or replay a public result. Do not delete existing targets to bypass the guard.

The historical publication command below is an archive record only. It could be
used only in a separately prepared fresh checkout/staging tree where all public
targets are absent, after restoring the complete original ignored `preview-v2`
bundle, verifying the externally supplied completion hash, authenticating a fresh
reporting replay, and obtaining separate coordinator approval of the actual PNG.
A normal checkout of the completed milestone already contains the public targets
and is not such a destination. The historical log path must also be absent.
Publication copies the reviewed preview bytes; it does not rerender its figures.

```text
.venv/bin/python -B scripts/report_psc_liver_programs.py publish --expected-reporter-sha256 2c502da9107deacab445df7f2e81d0496e0e479c4875e8bf4509df11bc892fba --preview-directory work/psc-liver-program-reporting/preview-v2 --expected-preview-completion-sha256 ef924593491e9c5b5b578e177a68f24de8f6af74c39206fe41594e60fed01211 --root-publication-ack ROOT_APPROVES_PSC_LIVER_AGGREGATE_PUBLICATION > work/psc-liver-program-reporting/root-review-v1/root-publication.log 2>&1
```

## 5. Interpretation and finite stop

The effects are mean log-CPM-plus-one composite differences, not raw fold changes,
TF activation, induction or cell-intrinsic regulation. All programs retain positive
equal mapped-member weights; gRNA1-up is not inverted. ECM is unsigned organization
membership, not a directional fibrosis-response score. Welch and bootstrap
intervals are pointwise, not simultaneous; omission ranges are not confidence
intervals. Missing fields remain missing, not observed zero effects or P=1 tests.

Whole-biopsy cell composition, clinical/stage/medication and technical confounding
remain. The paper's source-patient independence assertion is a working assumption,
not an independently recovered person/visit crosswalk. There is no healthy group.
Computational agreement is not equivalence, proof of generic inflammation or PSC
causality, independent mechanistic validation of MacroMap, biological replication,
or a clinical recommendation.

The epithelial IL17A-alone response-transfer endpoint remains **HOLD outside this
family**, not an empty or zero-scored endpoint. GSE159676 remains HOLD and NoPSC
remains PARK. The fixed liver question is closed after its verified primary,
mechanical replay and complete reporting. This guide authorizes no additional
source search, sensitivity, endpoint, annotation change or scientific rerun.

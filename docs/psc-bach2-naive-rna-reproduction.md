# Fixed BACH2 RNA check: reproduction

**16 September 2026. Not execution permission.** Primary, full-eight-matrix audit and scientific replay passed; 102 metadata/code seals and 16 source pins stayed unchanged. “Root” means coordinator, **not OS superuser**; no `sudo`. Authority: [protocol](psc-bach2-naive-rna-protocol.md) and [freeze](../config/psc-bach2-naive-rna-execution-freeze.json).

## Exact restoration; a clone is incomplete

The freeze required checkpoint `a1b0a6e1d880f9b97c49bdbc0358aef7db408a47` committed/pushed before execution. Do not repeat historical Git operations.

Restore exact bytes/paths from an authorized copy, using the complete 102-file inventory and 16 payload pins in `work/psc-bach2-naive-rna-freeze/private-seal-ledger.json` (hash bound by the freeze):

- `work/psc-bach2-naive-rna-plan/candidate-v4/plan.json`, adjacent `selection.private.json` and `public-schema.json`. Plan/hash are pinned in the freeze. Its preparation-only status stays unchanged, not authorization.
- Metadata/SDRF, eight feature/barcode axes, individual integrity/acquisition receipts, cached reference/provenance, candidate/decision/acceptance records and frozen code/tests. Exact matrix/archive payloads also require separate source-access permission; receipt hashes cannot replace missing payloads.
- Native authorizations, audit seals, freeze/ledger and archived completed outputs/verification receipts. Keep `work/` ignored and path/link guards intact. Selection, joins and individual outputs stay private.
- The [public schema](../config/psc-bach2-naive-rna-public-schema.json) and [verification protocol](../config/psc-bach2-naive-rna-verification-protocol.json), **and their ignored originals at ledger-listed paths**. Public copies are byte-identical, not substitute input paths.

Restore Python **3.12.14**, NumPy **2.5.3**, SciPy **1.18.1**, executable bytes, import origins and all **2,362** pinned installed files. Absolute interpreter: `/home/riels/Projects/Personal/psc-research/.venv/bin/python`. Version labels alone are insufficient. Preserve `PYTHONDONTWRITEBYTECODE=1`, `PYTHONHASHSEED=0`, and set each of `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS` to `1`.

Missing inputs, relocation/runtime drift or hash changes **fail closed**. No silent recompile, rehash, substitution, resealing or relaxed checks. Stop for a prospective coordinator decision if exact restoration fails.

## Native authorization, audit seal and completion

Native authorization: simple `.json` filename under **`work/psc-bach2-naive-rna-freeze/`**, maximum **8,192 bytes**. Role/name/size/symlink checks precede body access. An external SHA authenticates exact scientific bindings and fresh destination.

Audit authorization is separate: same namespace, maximum **16,384 bytes**, exact frozen `seal_bindings()` schema. It also binds verifier/adapter code/tests, design/protocol, native authorization and run/audit paths. Only actual `completion_sha256` is filled after authenticating completion/inventory; other per-run bindings are fixed beforehand. Root must separately permit private-result preloads and source-matrix access. Do not invent seal fields.

Exactly four native files: `authorization.used.json`, `measurements.private.json`, `public-aggregate.json`, `completion.json`. Last-written completion binds actual three payload bytes/SHA, inventory, destination and eight provenance digests. Partial/failed/extra/missing/modified outputs fail. Failed destinations remain consumed.

A fresh destination requires matching **external native authorization plus a new audit seal and source/preload permission**, not just another `--output-directory`. Never reuse existing paths.

## Completed commands — ARCHIVAL: DO NOT RERUN

Verbatim from `work/psc-bach2-naive-rna-freeze/root-execution-resume.json`, at original repository root/frozen environment. Paths/logs are occupied.

Primary:
```text
.venv/bin/python -B scripts/analyze_psc_bach2_naive_rna.py --execute --plan-directory work/psc-bach2-naive-rna-plan/candidate-v4 --authorization work/psc-bach2-naive-rna-freeze/primary-v1.json --authorization-sha256 b9703443a6585517e39ab5dc57d6287d5ccad34b3280e306acd268bc9244c886 --output-directory work/psc-bach2-naive-rna-plan/runs/primary-v1 > work/psc-bach2-naive-rna-freeze/primary-execution.log 2>&1
```

Independent audit:
```text
.venv/bin/python -B scripts/verify_psc_bach2_naive_rna_execution.py --real-audit --acknowledge-root-authorized-real-audit --plan-directory work/psc-bach2-naive-rna-plan/candidate-v4 --plan-sha256 b6826aa5adb6a3ff26c6f81f5fb16017df5e9580487a96afc2ff491554257b06 --run-directory work/psc-bach2-naive-rna-plan/runs/primary-v1 --authorization work/psc-bach2-naive-rna-freeze/primary-v1.json --authorization-sha256 b9703443a6585517e39ab5dc57d6287d5ccad34b3280e306acd268bc9244c886 --root-seal work/psc-bach2-naive-rna-freeze/primary-v1-audit.json --root-seal-sha256 053220f4f1a97adea7bde81c670ab5b87d25c595ed9767297cba780e53c97ce1 --completion-sha256 850ff604b819fb0aae660fabc9a42248d511a5a334b07764afdb46d3dd961c4f --output-directory work/psc-bach2-naive-rna-execution-review/primary-v1-audit > work/psc-bach2-naive-rna-freeze/primary-independent-audit.log 2>&1
```

Scientific replay:
```text
.venv/bin/python -B scripts/analyze_psc_bach2_naive_rna.py --execute --plan-directory work/psc-bach2-naive-rna-plan/candidate-v4 --authorization work/psc-bach2-naive-rna-freeze/replay-v1.json --authorization-sha256 dd99adc6a35846c87048af773b64a87a31fe0f10c2a933099f60fcdd5462808f --output-directory work/psc-bach2-naive-rna-plan/runs/replay-v1 > work/psc-bach2-naive-rna-freeze/replay-execution.log 2>&1
```

## Field-specific scientific replay

`root-mechanical-replay-verification.json` records:

- `measurements.private.json`: **byte-identical**; no authorization/provenance subtree exists.
- `public-aggregate.json`: only `provenance.authorization_sha256` changes. All other fields are exact.
- `authorization.used.json`: only `output_directory` changes.
- Completion: only destination, authorization provenance hash, authorization payload bytes/hash and public payload hash change; its own bytes/hash follow. Private payload pins and all other fields match.

No whole-provenance exclusion. Root’s early assertion wrongly expected private authorization provenance to change; private files were already byte-exact. Correcting the assertion caused no rerun or scientific/source/code repair. These are computational checks, not independent biological replication.

## Presentation-only reproduction

The [reporter design](../reports/psc-bach2-naive-rna-reporting-design.md) and `work/psc-bach2-naive-rna-reporting/reporter-handoff.json` fix seven inputs: root reporting handoff, public aggregate/schema, outer freeze, primary completion, independent verification, mechanical replay verification. No source/private quantities or joins; no new statistics.

The following are the exact recorded handoff commands: validate; the root presentation replay; and verification of the existing preview. **The displayed `root-replay-v1` is now occupied and must not be rerun.** For any later authorized presentation replay, choose a new simple preview name and use its returned completion SHA. Do not overwrite or repin an existing bundle:

```text
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py validate
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py preview --output work/psc-bach2-naive-rna-reporting/root-replay-v1
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py verify --bundle work/psc-bach2-naive-rna-reporting/preview-v1 --completion-sha256 785b9ebd45de81daafb3ef5d694d7830427fe632de52796a8441d6c8a9b85226
```

The **12-file** bundle: under `data/derived/psc-bach2-naive-rna/`: `public-aggregate.json`, `public-schema.json`, `groups.csv`, `endpoint.csv`, `influence-qc.json`, `aggregate-verification.json`, `mechanical-replay-verification.json`, `reporting-info.json`, `manifest.json`, `completion.json`; plus `reports/figures/psc-bach2-naive-rna.png` and `.svg`. Manifest binds ten payloads/schemas; completion binds manifest/inventory. All eleven presentation payload/manifest files replay exactly; only completion's `origin_preview_directory` and hash differ.

After root inspects bundle and PNG/SVG, these conditional commands require **separate** approval. Expected SHA requires exact template bytes:

```text
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py approval-template --bundle work/psc-bach2-naive-rna-reporting/preview-v1 --completion-sha256 785b9ebd45de81daafb3ef5d694d7830427fe632de52796a8441d6c8a9b85226
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py publish --bundle work/psc-bach2-naive-rna-reporting/preview-v1 --completion-sha256 785b9ebd45de81daafb3ef5d694d7830427fe632de52796a8441d6c8a9b85226 --approval work/psc-bach2-naive-rna-freeze/reporting-publication-v1.json --approval-sha256 e0b794b0bebab41b03ccdba794daaed789b9c7ced524ab975894a309c17406e9
```

Only root may save template stdout exclusively after review; never overwrite approvals/destinations. Publish writes completion last. Publication paths: [manifest](../data/derived/psc-bach2-naive-rna/manifest.json), [figure](../reports/figures/psc-bach2-naive-rna.png), [root narrative](../reports/psc-bach2-naive-rna.md). The reporter does not write narrative; this guide performed no publication.

## Stop

Keep one BACH2 C0/C1 endpoint, eight donors/four per arm, 36,601 RNA rows and all frozen Welch/omission/audit rules. Source-format units are **not certified untouched UMIs**. Pooled states do not resolve composition, causal, protein or clinical claims. Individual plots/omissions stay private; aggregates are not formal differential privacy. After this endpoint, diagnostics, audit, replay and reporting, **close regardless of result**. No new genes/states/cohorts/models or post-effect repairs.

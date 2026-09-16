# Reproduce the frozen MacroMap context analysis

## Scope and prerequisites

The [result report](../reports/macromap-program-context.md) distinguishes fixed-model stimulus-context performance from PSC or causal validation. The [protocol](macromap-program-context-protocol.md), [method specification](../reports/macromap-program-design.md) and [independent review](../reports/macromap-independent-review.md) are immutable pre-effect artifacts. Freeze: **2026-09-16T02:22:59.550114+00:00**; plan SHA-256 **`7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df`**; pre-effect checkpoint **`88c7464`**.

Commands require the exact pinned source/cache bundle and project `.venv`. Raw counts, source/line joins, learned references, predictions and audit arrays are deliberately ignored and are not supplied by a Git clone. Source URLs and versions are recorded in the method report and frozen source pins. Downloading a current response is not proof of recovering historical receipt/metadata bytes. Do not rewrite pins to bypass absent or changed caches. A genuinely changed source needs an explicit newly qualified plan, not silent recovery under this freeze.

Use Python3.12.14, NumPy2.5.3 and SciPy1.18.1 in the unchanged project environment/lock. Use one BLAS/OpenMP thread and deterministic Python hashing:

```text
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
```

The main normalized cache is 2,189,004,912 bytes, plus private score/fit/audit outputs. No GPU, new service, raw-read reprocessing or eQTL archive is required. The measured original primary run took about four minutes on this workstation; this is not a runtime guarantee.

## 1. Validate sources and preparation without effects

A new workspace with the exact source bundle can reconstruct the source-only preparation at its frozen path. **Do not run this over an existing preparation directory.** Existing source/preparation snapshots must not be overwritten.

```text
.venv/bin/python scripts/analyze_macromap_program_context.py --output work/macromap-design/preparation-runner-v3
.venv/bin/python -m unittest discover -s tests -p "test*macromap_program_context.py" -v
```

Expected preparation summary: `source-only-readiness.json`, SHA-256 **`eae76811c5b46388965f4f2c27624e650748b9d7adc216f75bff05b389e1f507`**. It is source-only, not execution permission. All ten preparation files replayed identically. Additional byte-identical preparation copies do not change the frozen v3 input selection.

The separate runner validates only unless `--execute` is supplied:

```text
.venv/bin/python scripts/run_macromap_program_context.py --plan config/macromap-program-context-execution-freeze.json --expected-plan-sha256 7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df
```

Expected status is `execution_plan_validated_no_expression_read`. The root pre-effect pass had 71 targeted tests. The dated full-workspace snapshot had 637 reported tests and seven existing missing-history/R-runtime errors, documented in [pre-effect verification](../reports/macromap-prefreeze-verification.json). Do not repair or rehash historical results to make that snapshot an all-pass claim.

## 2. Execute in the fresh frozen destination

The original plan uses `work/macromap-design/primary-execution-v1`. That path must not exist. The runner rejects reused destinations, aliases, literal parent traversal and reused cache namespaces. It authenticates source/code/runtime/preparation pins before and after numerical work. Failures remain explicit; no failed fold or phenotype is removed to obtain an available primary estimate.

```text
.venv/bin/python scripts/run_macromap_program_context.py --plan config/macromap-program-context-execution-freeze.json --expected-plan-sha256 7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df --execute
```

Expected completion file: `execution-complete.json`. This inventories the private artifacts and keeps every prescribed endpoint/failure status. It does not by itself constitute independent numerical verification. Every real output remains under ignored `work/` until an aggregate-only reporter is separately checked.

## 3. Independent read-only verification

Hash the newly generated completion receipt. Use that actual hash, not the original run's hash merely because numerical outputs agree. Supply the frozen verifier hash **`31b69fca6bcbe39e9204d3816203b39c5cb96762ef2e56ec9562dede77c71a7a`**. A verifier receipt destination must be fresh and outside the immutable input bundles.

The original successful command was:

```text
.venv/bin/python scripts/verify_macromap_program_context.py --preparation work/macromap-design/preparation-runner-v3 --preparation-sha256 eae76811c5b46388965f4f2c27624e650748b9d7adc216f75bff05b389e1f507 --execution-directory work/macromap-design/primary-execution-v1 --plan config/macromap-program-context-execution-freeze.json --plan-sha256 7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df --execution-complete-sha256 e85445627230c7bdae025df48c80e01c34952da0c476237bfe099f351cee1408 --verifier-sha256 31b69fca6bcbe39e9204d3816203b39c5cb96762ef2e56ec9562dede77c71a7a --output work/macromap-independent-review/postfreeze-verification-v1.json
```

For a new run, replace only the completion-receipt hash and fresh receipt destination as appropriate. The original receipt hash was **`e85445627230c7bdae025df48c80e01c34952da0c476237bfe099f351cee1408`**. The [public verification receipt](../reports/macromap-independent-verification.json) is byte-identical to the inspected aggregate-only private report; it contains no individual identities or model arrays.

The frozen scope includes all 80 paired-model cells, all expected predictions/losses, all primary/sensitivity endpoint points and conditional prediction intervals. Raw normalization independently reconstructs eight hash-selected samples over all 58,243 genes. Source-to-score reconstruction covers fold0 for both protocols/times/scopes, all nine programs and four hash-selected training plus four held-out pairs per context. All response points/statuses are checked; response-interval values are independently replayed only for CIL6h, both protocols, all programs/quantities/scopes. Other interval values and other-fold raw-to-score calculations are not claimed independently reconstructed.

The frozen independent optimum checks use gradient/coefficient tolerances, not a required `optimizer_success=True` flag. A separate post-fit transparency diagnostic records the solver flags returned by the same frozen routine. It must not alter that routine, its thresholds or any primary fit. The [diagnostic](../reports/macromap-independent-optimizer-diagnostics.json) retains 120 true and 40 false solver-success flags; all 160 numeric criteria pass. None of the training features triggers the fixed scale fallback. It captures fresh runs of the identical helper, not discarded original SciPy result objects. The13 dedicated native tests pass. The original command was:

```text
.venv/bin/python scripts/audit_macromap_independent_optimizer.py --plan config/macromap-program-context-execution-freeze.json --plan-sha256 7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df --execution-complete-sha256 e85445627230c7bdae025df48c80e01c34952da0c476237bfe099f351cee1408 --verifier-sha256 31b69fca6bcbe39e9204d3816203b39c5cb96762ef2e56ec9562dede77c71a7a --output reports/macromap-independent-optimizer-diagnostics.json
```

Use a fresh report destination for a rerun; do not overwrite the published diagnostic. It rejects changed model/code/plan/runtime pins and authenticates all 80 fitted-array payloads before running the frozen helper. A rerun's creation timestamp may differ; no solver message/status code is available from this helper.

## 4. Clean replay without a scientific change

Never overwrite the original plan or output. Copy the frozen plan into a NEW ignored JSON file and change **only `output_directory`**, for example to `work/macromap-design/replay-execution-v1`. Preserve every other field, including the scientific freeze timestamp, source/code/software pins and method contract. Record original/new hashes and the mechanical creation time separately; this is not a backdated scientific amendment.

The original authorized mechanical clone is described by [the replay receipt](../reports/macromap-replay-verification.json). Its plan SHA-256 was **`4f9b74a8f99c5cd7167b3167f43aefb826b6f3ce78bf99b705472bc5e101f13d`**. The exact original replay command was:

```text
.venv/bin/python scripts/run_macromap_program_context.py --plan work/macromap-design/root-mechanical-replay-plan-v1.json --expected-plan-sha256 4f9b74a8f99c5cd7167b3167f43aefb826b6f3ce78bf99b705472bc5e101f13d --execute
```

A local replay plan under ignored `work/` is not supplied by Git; reconstruct it from the public original as described above. Before execution, confirm that dictionary comparison finds only the `output_directory` difference. Do not add fields to make a failed scientific check pass.

For both completed runs, enumerate the actual `artifact_sha256` inventories and independently rehash every listed file with streamed SHA-256. Reject extra/missing paths or unexpected differences. The original comparison rehashed534 listed files (267 per execution): **266 corresponding artifacts were byte-identical**, including the normalized cache, all score/model arrays and aggregate results. `execution-start.json` differed only in plan hash/start time. The completion receipt differed only in those run fields, completion time and the start-receipt manifest hash. All endpoint dictionaries were identical. Timestamped receipts are not expected to be byte-identical.

## 5. Publish aggregate results only

Post-fit reporting is separate from the frozen science. It authenticates the plan, completion, selected payloads and successful verification before exporting complete aggregate panels and figures. It must preserve both identity scopes, both predictive variants, all registered conditions including unavailable PIC, all quantity components and fit flags. No model probabilities, source/line joins, donor responses or raw/normalized expression belong in Git.

The [aggregate manifest](../data/derived/macromap-program-context/reporting-manifest.json) pins 11 deterministic gzip tables, a compact verification summary and four figure files. With the manifest itself this is13 data files plus 4 figures. The [presentation replay receipt](../reports/macromap-presentation-replay-verification.json) rehashes both copies: all 17 artifacts are byte-identical. Root inspected both rendered PNGs and complete output schemas. The [stage record](../reports/macromap-stage-verification.json) preserves targeted tests, privacy checks and the full-workspace test snapshot.

The final reporter has 25 passing native tests. Rendering uses Matplotlib3.11.2, Pillow12.3.0, NumPy2.5.3 and zlib1.3.2 with Python3.12.14. No package was installed for reporting. The reporter reads only approved authenticated aggregates/QC; it does not reopen model arrays, source expression or per-line predictions.

The original first-publication command was:

```text
.venv/bin/python scripts/report_macromap_program_context.py --plan config/macromap-program-context-execution-freeze.json --plan-sha256 7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df --execution-complete-sha256 e85445627230c7bdae025df48c80e01c34952da0c476237bfe099f351cee1408 --verification-receipt work/macromap-independent-review/postfreeze-verification-v1.json --verification-receipt-sha256 3cfaf66e53481883f04916935e7c1945fe664caf32b5cf29f962fb192d665e96 --public-verification-summary reports/macromap-independent-verification.json --public-verification-summary-sha256 3cfaf66e53481883f04916935e7c1945fe664caf32b5cf29f962fb192d665e96 --data-output data/derived/macromap-program-context --figure-prefix reports/figures/macromap-program-context --replay-verification reports/macromap-replay-verification.json --replay-verification-sha256 683c9068ae62019bd4bae8d5ddaa087e7da638691d023468bbefe66d0a477877 --optimizer-diagnostics reports/macromap-independent-optimizer-diagnostics.json --optimizer-diagnostics-sha256 586168f8e9c8cf6e246aba7e235445297e22cf0d57f7536d2207bfa9028b1232 --publish
```

These public output paths already exist in a checkout containing the results. Do not overwrite them. For a local presentation-only replay, create a NEW ignored parent directory while leaving its `data/` child and all figure files absent. The original successful replay used:

```text
.venv/bin/python scripts/report_macromap_program_context.py --plan config/macromap-program-context-execution-freeze.json --plan-sha256 7a26c93bde72a7f4ffeb4256da6f06d50e2df1b3a8f69b6da73a9f19bd6127df --execution-complete-sha256 e85445627230c7bdae025df48c80e01c34952da0c476237bfe099f351cee1408 --verification-receipt work/macromap-independent-review/postfreeze-verification-v1.json --verification-receipt-sha256 3cfaf66e53481883f04916935e7c1945fe664caf32b5cf29f962fb192d665e96 --public-verification-summary reports/macromap-independent-verification.json --public-verification-summary-sha256 3cfaf66e53481883f04916935e7c1945fe664caf32b5cf29f962fb192d665e96 --data-output work/macromap-reporting/presentation-replay-v1/data --figure-prefix work/macromap-reporting/presentation-replay-v1/macromap-program-context --replay-verification reports/macromap-replay-verification.json --replay-verification-sha256 683c9068ae62019bd4bae8d5ddaa087e7da638691d023468bbefe66d0a477877 --optimizer-diagnostics reports/macromap-independent-optimizer-diagnostics.json --optimizer-diagnostics-sha256 586168f8e9c8cf6e246aba7e235445297e22cf0d57f7536d2207bfa9028b1232 --publish
```

Use a new namespace if that replay path already exists. Compare the actual13 data files and four SVG/PNG files, not only manifest claims. Neither command reruns scientific scoring, models or bootstrap endpoints. The reporter rejects reused destinations; an unmanifested partial output is not a completed report. Its hash-specific bounded privacy review is not a promise of universal handling of every malformed future schema or hostile concurrent filesystem actor.

Pooled response medians/means are descriptive without CIs. Predictive intervals and protocol-specific response intervals remain conditional and pointwise. Reporting introduces no response-P screen and does not replace the primary with the stronger nonoverlap sensitivity.

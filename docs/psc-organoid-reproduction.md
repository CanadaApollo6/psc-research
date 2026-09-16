# Reproducing the PSC organoid IL-17A analysis

The [result](../reports/psc-organoid-il17.md) comes from the immutable [execution freeze](../config/psc-organoid-il17-execution-freeze.json), SHA256 `89650c5f72875cd26913ada15f6874194b4f4694fd05b182418b1bf00d44a58b`, frozen September 16, 2026 at 00:41:14.393412 UTC. Pre-effect commit: `52013ec`. Do not edit the frozen protocol, source qualification, primary/verifier scripts, tests, reviews or lock to make a replay pass. A substantive change needs an explicit amendment and a separate output namespace.

## Environment and data boundaries

Use the project's `.venv/bin/python`, not another Python kernel. The unchanged [blood environment lock](../requirements-psc-blood.lock) is reused. The recorded runtime is Python 3.12.14, PyDESeq2 0.5.4, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1, anndata 0.13.3.post0, formulaic 1.2.2 and joblib 1.6.0. The reporting manifest additionally records Matplotlib. No new dependency was installed for this stage.

The primary plan specifies four PyDESeq2 workers. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1` and `PYTHONHASHSEED=0` in each process environment. These settings are pinned in the plan. Commands below assume them and the repository root.

Original inputs and private donor/fit arrays are deliberately not in Git. A complete numerical replay needs their exact bytes or must rebuild them from the public releases. Merely opening the committed aggregate results is not a new independent verification. Allow several GB of free disk and 8–16 GB RAM for source preparation; source transfer in the recorded acquisition was 879,495,941 bytes. No FASTQ, account, researcher contact or private data are needed.

## 1. Rebuild qualified source inputs only in a fresh checkout/workspace

Keep the versioned qualification JSON/report unchanged. The acquisition script has its own explicit acquisition actions; numerical qualification itself is offline. New download receipts will have new retrieval times. If an upstream body no longer matches a frozen source checksum, stop rather than silently substituting it or rehashing the plan.

```sh
.venv/bin/python scripts/qualify_psc_organoid_inputs.py --acquire --acquire-original
.venv/bin/python scripts/qualify_psc_organoid_inputs.py --qualify --qualify-original --work-dir work/psc-organoid-inputs --output work/psc-organoid-inputs/reconstructed-source-qualification.json
```

Do not run these over the completed workspace's closed source artifacts. The reconstruction output is deliberately ignored and does not overwrite `data/derived/psc-organoid-input-qualification.json`. `--finalize` is not needed to rebuild the count inputs and is not used here to replace the historical qualification/validation receipt. That receipt remains a dated record, including its original separate numerical checks.

The frozen analysis consumes:

| Ignored input | SHA256 |
|---|---|
| `work/psc-organoid-inputs/raw-rna-pseudobulk.tsv.gz` | `2c4fe6fb5728ee18e6c0542285f62f2c49efbed518f5822d34ca245c657da4bb` |
| `work/psc-organoid-inputs/raw-rna-features.tsv` | `b6b1148b3f91de1f9078e49e876485e814b72a963ab45816e3a6673a05b8f010` |
| `work/psc-organoid-inputs/raw-rna-libraries.tsv` | `d3f0fddf8a269c89ee2af7a54c8a895539b34051cda28d2a0c7ac1b8e54630c0` |

The cell-source joins are separately pinned in the freeze. All 36,601 original feature rows and exactly the 46,343 published retained cells are used. The merged SCT-consistent matrix is never supplied as raw RNA counts. The complete-original-raw-universe versus historical filtered-RNA distinction is part of the estimand, not a repair to be guessed during replay.

## 2. Validate into new directories

Both output and work directories must **not exist**. Never reuse primary, validation or replay directories, including aliases. In the recorded run, preflight found 17,897 eligible genes and a 17,463-gene all-positive ratio basis without computing effects.

```sh
.venv/bin/python scripts/analyze_psc_organoid_il17.py --plan config/psc-organoid-il17-execution-freeze.json --plan-sha256 89650c5f72875cd26913ada15f6874194b4f4694fd05b182418b1bf00d44a58b --output-dir data/derived/psc-organoid-il17-validation --work-dir work/psc-organoid-il17-validation --validate-only
```

These commands use the exact published plan, not a copied plan with new method or input hashes. The source-generation work directory must match the plan's input paths. A fresh location for inference outputs does not alter the scientific plan.

## 3. Execute a clean replay

The committed `data/derived/psc-organoid-il17/` is the primary aggregate run. Do not overwrite it. Use new replay directories:

```sh
.venv/bin/python scripts/analyze_psc_organoid_il17.py --plan config/psc-organoid-il17-execution-freeze.json --plan-sha256 89650c5f72875cd26913ada15f6874194b4f4694fd05b182418b1bf00d44a58b --output-dir data/derived/psc-organoid-il17-replay --work-dir work/psc-organoid-il17-replay
```

The public directory contains six complete-feature CSV.gz tables, `organoid-qc.json` and **`organoid-output-manifest.json`**. The within-group file contains two full feature universes. The manifest inventories the seven payloads; it is the eighth public file. Low-count, zero, unsupported and flagged features remain present. The fit arrays, coefficients, donor changes, all allocation/leaveout values, logs and selected-library map stay in ignored work.

## 4. Independently verify without overwriting a receipt

```sh
.venv/bin/python scripts/verify_psc_organoid_il17.py --plan config/psc-organoid-il17-execution-freeze.json --plan-sha256 89650c5f72875cd26913ada15f6874194b4f4694fd05b182418b1bf00d44a58b --output-dir data/derived/psc-organoid-il17-replay --work-dir work/psc-organoid-il17-replay --verification-json reports/psc-organoid-independent-verification-reproduced.json
```

The receipt must be NEW and outside protected primary/input/work/code directories. Without `--verification-json`, only an aggregate summary is printed. The verifier must not alter primary outputs or private arrays.

The original receipt has 397 successful checks. It independently reconstructs the complete inferential arithmetic and the fixed-dispersion conditional coefficient panel. It does not independently fit dispersion estimation or certify biological independence/calibration. Nine conditional optimizer success flags are false even though all 16 selected genes meet the separate frozen numerical criterion. Those facts must not be replaced with “all optimizers converged.”

The recorded [replay receipt](../reports/psc-organoid-replay-verification.json) compares all eight public and all 13 private artifacts by SHA256 and bytes. Compare each replay file to its primary counterpart; do not regenerate expected hashes from a changed primary. Replay is computational reproducibility, not replication in another cohort.

## 5. Render aggregate results without private inputs

This stage can run using committed aggregate tables and the historical independent receipt alone. It performs joins, summaries and plotting; it does not rerun verification or fit new models. Use fresh output/figure destinations separate from each other and from the primary directory.

```sh
.venv/bin/python scripts/report_psc_organoid_il17.py --input-dir data/derived/psc-organoid-il17 --verification-json reports/psc-organoid-independent-verification.json --output-dir data/derived/psc-organoid-il17-summary-replay --figure-prefix data/derived/psc-organoid-il17-figure-replay/psc-organoid-il17
```

The reporting wrapper verifies the receipt's public-manifest hash and every payload, preserves exact source labels/values in the candidate and diagnostic CSVs, and produces a finite machine-readable summary. No source/core result file is rewritten. The figure uses only aggregate public estimates and pointwise intervals. Candidate panels are explicitly post-result displays, not new confirmatory tests.

## 6. Tests and preservation

```sh
.venv/bin/python -m unittest discover -s tests -p 'test_psc_organoid*.py' -v
.venv/bin/python -m unittest tests.test_verify_psc_organoid_il17 tests.test_report_psc_organoid_il17 -v
.venv/bin/python -m unittest discover -s tests -v
```

The first two commands cover 84 source/analysis +24 independent-verifier +17 reporting tests at this milestone. The complete repository still has seven known historical missing-cache/R-runtime errors. The [stage record](../reports/psc-organoid-stage-verification.json) contains the exact full-run count and error names. Do not “repair” frozen prior scientific outputs to make unrelated historical tests pass.

After reproduction, preserve the eight-donor limit, procedure controls, count-universe difference, model/Welch disagreement, all stage flags and unavailable boundary contrasts. Neither q values nor model-generated scores are clinical recommendations.

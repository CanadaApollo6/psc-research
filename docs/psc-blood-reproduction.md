# Reproduce the PSC blood reanalysis

## Scope and records

Read the [findings](../reports/psc-blood-replication.md) and [prospective protocol](psc-blood-replication-protocol.md) before interpreting the tables. This is public-data reuse, not a diagnostic validation or causal claim.

- [Source URLs, retrieval dates, sizes and hashes](../config/psc-blood-sources.json).
- [RNA execution plan](../config/psc-blood-rnaseq-plan.json), [array execution plan](../config/psc-blood-array-plan.json), and [pre-effect code/method freeze](../config/psc-blood-execution-freeze.json).
- [Mechanical integration receipt](../config/psc-blood-integration-plan.json): source-result hashes and column aliases became knowable after the assay runs; conjunction methods and code were already frozen.
- [Selection-time protocol snapshot](archive/psc-blood-replication-protocol-selection.md) preserves the hash referenced by the original [path selection](../config/psc-blood-path-selection.json). The current protocol incorporated the 30-sample count rule, missing-P family handling, HC3 reference distribution and dispersion fallback **before effects**. It is not the same byte snapshot as the initial selection narrative.
- [Runtime lock](../requirements-psc-blood.lock), [baseline limits](../reports/research-runtime-baseline.json), [independent checks](../reports/psc-blood-independent-verification.json), [replay receipt](../reports/psc-blood-replay-verification.json) and [method review](../reports/psc-blood-method-review.md).

The first execution freeze was **2026-09-15T22:56:48.344790+00:00**. Its SHA-256 is `62ce2b14fe425ea33531de0a938395542edc7de0e9b660ecaa9274dabb47a236`. This is a local prospective record, not an externally registered protocol.

## Environment

Run from the repository root. The new analysis uses Python **3.12.14**, PyDESeq2 **0.5.4**, numpy **2.5.3**, pandas **3.0.5**, scipy **1.18.1** and statsmodels **0.15.0**. The full lock includes the repository's existing scientific dependencies. It does not install an R runtime or historical caches.

```sh
uv venv --python 3.12.14 .venv
uv pip sync --python .venv/bin/python requirements-psc-blood.lock
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
mkdir -p work/research-path
```

The fit uses four local CPU workers. No GPU, paid compute, model API, account enrollment or researcher contact is required. On this machine the first RNA and array commands took about 14.4 and 12.2 seconds, respectively; times are not portable performance promises.

## Restore or check inputs

Original response bodies belong in ignored `data/raw/public-data-discovery-options/`. The manifest includes the three GSE177044 count exports to preserve the full source-partition audit, but the **only RNA input to the primary fit is `GSE177044_raw_countsPSC.csv.gz`**. GEO metadata and source articles are pinned separately.

The complete `data/derived/gse84161-gene-map.csv.gz` export is already versioned. Use **all 91,748 original relationship rows**, not its previous platform-specific eligibility flags. The historical upstream raw cache is absent; this analysis independently recalculates identity using the complete pinned export and does not claim fresh upstream annotation validation.

```sh
# Cached-input integrity check; no network requests:
.venv/bin/python scripts/acquire_psc_blood_inputs.py --check-only

# On a fresh checkout, omit --check-only to restore missing public bodies:
.venv/bin/python scripts/acquire_psc_blood_inputs.py
```

The restoration helper verifies exact size and SHA-256, refuses to overwrite a changed existing input, stops if upstream bytes differ, and records successful fresh retrieval dates separately under `work/research-path/restored-input-receipts.jsonl`. It does not update the historical source pins. The cached check was executed successfully here; an extra full network re-download was not needed. GEO response bodies can change: a mismatch requires investigation or a separately versioned source release, not a silently revised frozen plan.

Full matrices are modest public processed inputs. The pinned GSE177044 metadata declares GRCh37, STAR v2.6.1d and FeatureCounts/subread v1.6.4. The original count-annotation release was not independently resolved; no release number is invented. Stable-ID lookup against Ensembl 116 is not coordinate lift-over. See the [processing receipt](../reports/psc-blood-source-processing.json). No FASTQ reconstruction was used. The GSE119600 `RAW.tar` is an annotation archive, **not subject expression**; do not substitute it for the pinned normalized series matrix.

## Run the fixed analyses

```sh
.venv/bin/python scripts/analyze_psc_blood_rnaseq.py --plan config/psc-blood-rnaseq-plan.json --plan-sha256 da3a43d0bd73e1cee649358b83179372f2695d584bb8fcac98b9bdc43d401be1 --output-dir data/derived/psc-blood-rnaseq --work-dir work/psc-blood-rnaseq > work/research-path/rnaseq-run.log 2>&1

.venv/bin/python scripts/analyze_psc_blood_array.py --plan config/psc-blood-array-plan.json --plan-sha256 5c562f5045d55dbcf4744733a9de0bbc024e1b4514d4b43f2be95a68b48a1ef3 --output-dir data/derived/psc-blood-array --work-dir work/psc-blood-array > work/research-path/array-run.log 2>&1

.venv/bin/python scripts/summarize_psc_blood_replication.py --plan config/psc-blood-integration-plan.json --output-dir data/derived/psc-blood-replication > work/research-path/integration-run.log 2>&1

.venv/bin/python scripts/report_psc_blood_diagnostics.py > work/research-path/diagnostics-run.log 2>&1
.venv/bin/python scripts/plot_psc_blood_replication.py
.venv/bin/python scripts/package_psc_blood_results.py
```

The first two commands independently validate the frozen source plans, identities and complete sample joins. The integration receipt pins the raw primary Wald/Welch P fields, positive PSC-minus-reference effects, every consumed result/crosswalk file, and the executing integration script. The runner's general manifest-completeness guard is not universal; this exact receipt was separately checked. Similarly, its duplicate-stable-ID edge case is absent in this verified source, which has 63,677 valid unique IDs. A future dataset requires explicit source-feature audit handling; do not reuse these wrappers as unconstrained generic importers.

The RNA command preserves 63,677 rows, including 48,632 low-count exclusions and four unavailable Wald tests. The array command preserves its complete 20,749-gene family in all four contrasts. The shared family remains 11,074, including three unavailable constituent tests. No result file is limited to significant genes.

### Public artifacts and local patient-level work

| Location | Contents |
|---|---|
| `data/derived/psc-blood-rnaseq/` | All source features, full primary/HC3/subgroup results, leave-plate diagnostics, model QC and output hashes. Lossless `.csv.gz` copies are versioned; the larger exact uncompressed CSVs are generated locally and ignored by Git. |
| `data/derived/psc-blood-array/` | Complete probe annotation/join audit, Entrez gene universe, all four contrasts, aggregate sample/method checks. |
| `data/derived/psc-blood-replication/` | Complete shared-gene results, all-source RNA coverage, three conjunction summaries, all 93 BH candidates with flags, and all seven fixed benchmarks. |
| `reports/figures/` | Aggregate PNG/SVG plots. |
| `work/psc-blood-rnaseq/`, `work/psc-blood-array/` | Sample joins, design matrices, size factors and subject-level expression. Not versioned. |

Compression is a packaging transformation only. `compressed-artifacts.json` verifies that decompression exactly recovers each frozen uncompressed output SHA-256. If only the versioned gzip outputs are available, the complete regeneration commands above restore the original CSV paths needed by the integration receipt. The gzip files can also be read directly by pandas for inspection.

The outer RNA stderr log matters: joblib worker numerical warnings are not captured by the main-process warning list. Do not infer a warning-free fit from `warnings: []`. The warnings and convergence flags are summarized without changing the primary result. No post-effect convergence exclusion or normalization fallback was introduced.

## Independent verification

```sh
.venv/bin/python scripts/verify_psc_blood_replication.py --run --execution-freeze config/psc-blood-execution-freeze.json --execution-freeze-sha256 62ce2b14fe425ea33531de0a938395542edc7de0e9b660ecaa9274dabb47a236 --array-dir data/derived/psc-blood-array --rna-dir data/derived/psc-blood-rnaseq --rna-work-dir work/psc-blood-rnaseq --replication-dir data/derived/psc-blood-replication
```

The independent script does not import the three analysis scripts. It reconstructs the arrays from the original positive linear intensities, original GPL probe mappings and source adult metadata. It independently checks **82,996** Welch comparisons, full-family adjustments, all source joins/normalization and complete conjunction corrections. The fixed RNA check subset contains **20 SHA-256-ranked genes plus all seven original benchmarks**; full, subgroup and leave-plate HC3 fits are compared with statsmodels. All available NB P values are checked arithmetically against reported statistics, with Cook/missingness/convergence accounting. **There is no independent NB likelihood refit.**

All **236 checks** passed at the declared tolerances. These are computational checks, not 236 biological replications. Verifier implementation repairs (a bracket typo and draft-schema/coverage-key adaptations) are disclosed in its report; no frozen analysis result or numerical tolerance was changed to obtain agreement.

## Clean-output cached-input replay

The executed replay wrote both arms into separate ignored `work/psc-blood-replay/*-aggregate` directories, with distinct `*-patient` work directories. It then made a mechanical integration receipt pointing to those **replayed** assay outputs and verified their existing pins. It did not merely run the summary on the first-run files.

All **12 checked artifacts**, including every full gene-level table, were byte-identical. The integration JSON's only change was `plan_sha256`, as expected from different input paths. The scientific counts and table hashes agree exactly. This is a local same-environment result; numerical differences on another platform still require the declared tolerance and semantic checks. The replay did not use network access in its analysis path.

## Tests and historical environment failures

```sh
.venv/bin/python -m unittest discover -s tests -p 'test_psc_blood*.py' -v
.venv/bin/python -m unittest discover -s tests -v
```

The **67 new analysis tests pass**. The complete suite reports **322 tests**, with **319 individual test successes, three test errors and four class-setup errors**. The same seven errors occurred in the baseline before these changes. Missing items are historical ETS2 XML/ZIP caches, the old GSE84161 SOFT cache and `work/gse84161-normalization/runtime/bin/Rscript`. They are not repaired by installing the new Python lock. No new catalog parsing was introduced; the earlier catalog and frozen science remain unchanged.

Do not report the full historical suite as passing in this environment. Source downloads, numerical qualification, a successful replay and independent calculations also do not establish causality or clinical utility. The result remains the limited association statement in the findings report.

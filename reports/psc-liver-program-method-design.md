# Pediatric liver fixed-program method design

## Status and the narrower question

This is a **source-blinded method and guarded-runner preparation**, not a statistical freeze or a real analysis. No real expression payload bytes, library totals, normalized values, scores, effects or clinical covariates were inspected or computed here. Tests use generated synthetic values. Source-only preparation reads identifiers, source-membership metadata and the accepted sample join.

Root approved a **narrower ETS2/comparator plus constitutive ECM tissue-context question** before any liver effects. It has the nine accepted original/nonoverlap programs plus `reactome_ecm_organization`: **ten programs × three contrasts = 30 Core hypotheses**. This decision is recorded as `agentmsg_8bf94a77-6693-40fa-b017-395f3ce33d56` and in `work/psc-liver-program-plan/root-source-panel-decision.json`.

The epithelial IL17A-response source remains **HOLD**. It is not an empty program, zero score, curated-pathway substitute or member of this narrower executable panel. This design **does not fulfill the blocked epithelial IL17-response-transfer question**. The original ETS2/comparator membership preparation was queued before MacroMap effects. This work is **not independent mechanistic confirmation of MacroMap**.

Root also approved the prospective statistical choices in `agentmsg_d84e9242-be6a-4f66-9ebb-cb9d839f0ca6`, recorded in `work/psc-liver-program-plan/root-method-decisions.json`. These include ordinary one-zero-variance Welch, the fixed bootstrap seed and the source-patient-sample independence assumption. The proposal remains **unfrozen and unauthorized for real execution**. Both independent reviews passed. The proposal still carries no root acceptance record or root-frozen compiled-membership digest. Source readiness and synthetic success cannot create that permission.

## Source and identity boundaries

- Study: **GSE303271**, complete released **41,096-row × 64-column** RNA count interface.
- Groups remain **17 PSC / 17 ASC / 30 AIH**. No group recoding or merging occurs.
- The material is whole liver biopsy. It is not a sorted-cell assay.
- The source article asserts patients. There is no independent person/visit crosswalk. We call the units **source patient-samples**, not verified independent donors.
- Sample-linked age, sex, treatment, disease stage, IBD and batch covariates remain unavailable. No group summaries are assigned to individuals and no adjusted model is invented.
- The paper declares GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4. The exact historical counting GTF/provider/release/checksum and assignment settings remain unresolved.
- Core, `exact_symbol_ensembl_bijection`, is primary. `plus_reciprocal_entrez` is the separate Strict identity sensitivity. Full-reference identity checking does not reconstruct historical biological loci.
- There is no alias rescue, case folding, ambiguity collapse, Entrez-based rescue of an official symbol, or zero filling of absent genes. Published Entrez/source IDs are provenance, not alternate routes selected by coverage or effects.
- The accepted identifier, liver and organoid files remain byte-pinned and unchanged. No old producer is rerun. Adult GSE159676 expression HOLD and NoPSC atlas PARK are unchanged.

Source URLs, retrieval dates, releases and transformations remain in the pinned qualification manifests and reports. This method makes no network request or new retrieval. The fixed reference lineage is Ensembl 116 (June 2026), with the full human relation audit already closed.

### ECM source and source-only mapping

ECM is the official MSigDB `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION`, M610 / R-HSA-1474244, MSigDB 2026.1.Hs / Reactome 95. The original denominator is **321 official symbols**, not 340 original source-ID records and not a mapped subset. Each official symbol has one unit weight. All 340 source IDs remain in the source package as provenance. Literal **`ETS2` is absent** from the official 321-symbol source. No member was removed and no new ETS2 exclusion is applied. `UNSIGNED` marks absent perturbation direction; it does not mean induced or activated.

A source-only exact-symbol join to the closed full-reference feature audit finds:

| Identity tier | Mapped unique members | Original denominator | Annotation gate |
|---|---:|---:|---|
| Core | 301 | 321 | Pass |
| Strict | 299 | 321 | Pass |

Core excludes 20 members with multiple Ensembl relationships. Strict excludes the same 20, plus one shared-Entrez relationship and one absent Entrez relationship. All **642** member/interface statuses and indices agree with the independent mapper. These are **identifier coverage results**, not expression, detection, pathway activity or liver effects. All official source members retain an explicit mapping status. Missing feature indexes are `null`, not zero. The complete mappings and all twenty tier/program specifications remain in ignored work artifacts.

The external adapter binds a pinned, metadata-only source-worker handoff to the exact qualified manifest and CSV. It checks the manifest's qualified program/version, species, denominator, unsigned semantics, distinct member/provenance namespaces and unit weighting. The complete HTML/JSON source pins are also checked locally, without retrieval. The closed legacy upstream metadata/source pins are checked; the real RNA payload pin is deferred to authorized execution. Only this qualified external adapter is registered. An arbitrary source path, RNA payload or future signed source is refused. A new source requires a new qualified and reviewed adapter; a claimed `qualified` flag alone is insufficient.

Pinned source package:

- `config/psc-liver-external-program-sources.json`: `2428cf428187581d992dc2b507305b259465db6259aec19b8b658614c0773305`.
- `data/derived/psc-liver-external-program-membership.csv`: `8cf6573a97d01cc8e83558634bf62b6d48d9be015e67f50d15868071c18b7202`.
- `work/psc-liver-external-programs/method-schema-handoff.json`: `69d7653d5829837eed563cef2c5d22eca1394c03b62a2942909f18fef50e3652`.

## Prospectively approved computation; execution still blocked

### Released-universe normalization and scores

For feature `i` and source patient-sample `j`:

```
library_total[j] = sum(count[r,j] for every one of the 41,096 released rows)
CPM[i,j] = 1,000,000 * count[i,j] / library_total[j]
log_expression[i,j] = log2(CPM[i,j] + 1)
score[program,j] = mean(log_expression[i,j] for the fixed mapped members)
```

The code uses the numerically stable equivalent `log1p(CPM)/log(2)`. Every approved real endpoint is a **positive equal-member mean**. Source companion labels remain as published; gRNA1-up is not inverted for a narrative. A score difference is a difference in mean log-transformed transcript abundance. It is not a gene-level log fold change, pathway/TF activity estimate or experimental cytokine response.

Unannotated and nonprogram rows contribute to library totals. Identifier-only background pools are not normalization inputs. There is no abundance matching, learned control weight, gene-wise variance standardization, observed detection filter or variance filter. A mapped gene measured as zero remains a member. An absent gene does not become a measured zero.

Input counts must be finite, nonnegative integers with the frozen shape. Counts and totals must be at most `2**53`; every total must be positive. The limit is checked before integer-to-float conversion. Totals use bounded exact integer blocks, so `[2**53, 1]` cannot be accepted after float rounding. Corruption, nonfinite values, a zero library or unsafe integer precision fails the run. There is no deletion, repair, rounding or sample selection.

### Fixed coverage and nonoverlap

Original intended denominators from accepted source preparation remain unchanged. Coverage uses the exact inequality `5*mapped >= 4*original` and at least ten unique mapped members.

`ets2_g1_dn_without_comparators` remains the primary gRNA1-down set minus **only the original four comparator programs**: inflammation, interferon gamma, oxidative stress and apoptotic signaling. ECM is not added to this subtraction. Its gate remains the parent 927-member coverage gate plus at least ten retained members. It does not get a new 80% denominator.

Coverage is determined once from identifiers, not by group, expression, bootstrap draw or omission. Failure leaves an unavailable endpoint in the approved family. It does not trigger a source substitution, different identity route or smaller denominator.

Generic balanced signed-score functions have synthetic tests, but **no signed endpoint is approved or executable in this ten-program question**. They are not an IL17 substitute or a commitment for a later qualified source.

### Contrasts and Welch inference

The fixed ordered contrasts are:

1. **PSC minus AIH**, primary; the focal program is ETS2 gRNA1-down.
2. PSC minus ASC, secondary.
3. ASC minus AIH, secondary.

For source-sample score vectors, use unequal-variance Welch mean differences, sample variances (`ddof=1`), actual Welch–Satterthwaite degrees of freedom, a two-sided t survival-function P value and a **pointwise 95% t interval**. No pooled df or normal critical value is substituted. The three contrasts have rank two, but all three remain in the fixed family.

Exact constant vectors are detected before variance calculation. This avoids spurious variance and tiny P values caused by rounding `mean([2.1]*30)`. Scaled variance contributions avoid unnecessary df underflow. There is no epsilon variance floor or arbitrary small-variance threshold.

Root chose **`welch_if_one_positive`**. If exactly one group has zero variance and the other has positive variance, use ordinary Welch and flag the zero-variance group. Its df is the actual Satterthwaite df. If both groups have zero variance, keep the finite descriptive difference, but leave raw P, df and CI unavailable. The alternative conservative policy exists only in generic synthetic helpers; it is refused by the real-plan validator.

Bad inputs and unavailable inference are distinct. Invalid or nonfinite counts/scores fail data integrity. With valid finite scores and a finite descriptive difference, numerical variance or inference failure retains that difference and the full hypothesis row, but leaves **raw P/CI/df null**. Nonconstant variance underflow has an explicit numeric status; it is not called exact zero variance. Tail-probability underflow has `unavailable_numeric_tail_underflow` and a retained flag, never an observed P of zero or an epsilon substitute. There is no fabricated observed P of one for a degenerate test. A valid ordinary test can naturally return P=1 at a zero t statistic.

### Multiplicity and missingness

The sole discovery family is **all 30 Core endpoint/contrast rows**. Report both BH and BY, without selecting the more favorable adjustment. BH has dependence assumptions; BY handles arbitrary dependence of valid marginal P values. Neither repairs confounding, false independent-unit assumptions or invalid marginal tests.

An unavailable approved test has:

```
p_raw = null
p_for_correction = 1
p_was_substituted = true
q_bh = null
q_by = null
```

Only the complete correction vector uses the substitute 1. Its numeric adjusted placeholders are retained privately for audit, not presented as observed P or q values. The Welch layer explicitly classifies numeric inference failures as unavailable. An invalid numeric value supplied as if it were a valid raw P to the correction layer is an error, not silently recoded missingness. Strict retains the complete diagnostic grid and raw defined quantities, but has no discovery q values or discovery eligibility.

### Fixed bootstrap and omission diagnostics

Root approved **10,000** within-diagnosis source patient/sample resamples with replacement, **PCG64 seed `2026091602`**, group order **PSC, ASC, AIH**, and original source-column order within each group. Intervals use linear 2.5/97.5 percentiles. Each group draw is reused across every program, both identity tiers and all three contrasts. This preserves their joint dependence. Scores, memberships and normalization stay fixed. Constant draws are retained, not filtered by Welch variance rules. Exact constant decimal means avoid a spurious degenerate interval away from zero. Bootstrap outputs are pointwise sensitivity intervals, **not additional P values or independent confirmation**. Degenerate intervals are retained and flagged; they are not proof of zero population uncertainty. These intervals do not include annotation uncertainty, unknown clinical confounding or unknown person dependence.

Every one of the 64 source units is omitted for every endpoint, tier and contrast: **3,840 planned omission rows**. A unit outside a contrast gets an explicit unchanged row. For an included unit, recompute group summaries and Welch quantities from the remaining fixed scores. Do not renormalize other samples, remap genes, rerun coverage or drop an influential unit from the main result. Strict sign reversal means opposite nonzero signs; zero is not called a reversal. Compare signs directly, not by an underflow-prone product. Omission P values are diagnostic only and have no discovery q values.

**Public omission stability uses only in-contrast units:** 47 for PSC−AIH, 34 for PSC−ASC, and 47 for ASC−AIH. `relevant_unit_denominator`, defined-effect counts, effect ranges, maximum changes and reversal counts all use that same relevant subset. The unchanged third-group rows remain private; they cannot inflate a robustness denominator. Missing effects do not imply zero reversals.

Root accepted an explicit **source-patient-sample independence assumption**, not a verified unique-donor claim. A known unresolved repeat person blocks this implementation; it does not cause automatic aggregation or sample selection. The statistical plan itself still needs root freeze and separate real-execution permission.

## Guards, output contract and current commands

Current permitted commands:

```sh
.venv/bin/python -B -m unittest discover -s tests -p test_psc_liver_programs.py -v
.venv/bin/python -B scripts/analyze_psc_liver_programs.py --write-proposal
.venv/bin/python -B scripts/analyze_psc_liver_programs.py --preflight
```

`--write-proposal` records the existing prospective root panel and method decisions, but writes `state=draft_not_frozen` and no real-execution approval. It cannot write a root freeze. `--preflight` uses accepted metadata only. It never hashes or loads a real RNA payload. Its expected payload checks are explicitly deferred until separately authorized execution.

A future `--execute` requires **all** of:

- A designated root plan JSON and its explicitly supplied exact SHA-256.
- `state=frozen`, root role/reference/time, and a separate real-expression execution permission.
- The explicit acknowledgement `ROOT_AUTHORIZED_REAL_EXPRESSION_EXECUTION`.
- The exact root ten-program panel/order, three contrasts, Core/Strict roles and narrowed question.
- Both root decision references, the exact approved method/seed/order, accepted source-patient-sample assumption, and no known unresolved repeats. A different seed or the conservative one-zero policy is refused.
- Exact native Python/NumPy/SciPy versions, implementation pins and a matching accepted independent synthetic/guard review receipt.
- All closed metadata and qualified external source pins, and the exact preflight compiled-membership digest.
- An unused directory strictly under `work/psc-liver-program-plan/runs/`, with protected-path and symlink checks.
- Exact raw-count and accepted count-cache hashes, checked **after authorization and before numeric loading**.

The only plan paths accepted are `config/psc-liver-program-plan.json`, `config/psc-liver-program-freeze.json` and the owned unfrozen `work/psc-liver-program-plan/proposed-plan.json`. A count matrix, arbitrary source/QC JSON or symlink cannot masquerade as a plan. Source readers also have a fixed qualified-path whitelist. Pinned input leaves and ancestors cannot be symlinks. All ancillary readers share a pre-read alias guard, including designated plans and proposal implementation hashes. Their inputs cannot be hard links to the reserved RNA payloads. This check uses filesystem metadata, not payload bytes. These controls prevent accidental misuse; a root-owned JSON is a workflow authorization, not cryptographic identity proof or a security sandbox against a malicious Python caller.

No current plan meets the execution guards. **No real `--execute` was run.** Negative guard tests use generated temporary fixtures and sentinel loaders, never real values.

A future run writes only a fresh ignored private directory. It retains the exact root-plan bytes, source-unit IDs, scores, library totals, all omissions, bootstrap draws/differences, compiled memberships and correction vector privately. The aggregate-public candidate contains group/endpoint summaries and complete statuses, no source-unit IDs or per-unit values. It is not automatically copied into Git or a versioned report.

Completion is explicit. Pre-execution guard refusal creates no result directory. A later input race, data-integrity failure or incomplete computation has a failed-run receipt and a complete unavailable endpoint inventory, not a completion receipt. A directory without a valid completion receipt is incomplete. There are no partial-family discoveries or silent bad-draw deletions.

### Machine-readable statuses

- Modes distinguish source-only preflight, synthetic validation and separately authorized real execution.
- Twenty tier/program specifications produce **60** endpoint/contrast rows. Core correction has exactly **30** rows; Strict has none.
- Coverage, descriptive effect, Welch, bootstrap and omission availability are separate.
- Missing numeric values use standard JSON `null`, never NaN/Infinity or an invented zero.
- Source-only rows have `not_executed_source_only` and no correction vector. Passing preflight is not passing inference.
- Actual data-integrity failures are not converted into ordinary missing endpoints. Inference-only numeric failures retain finite descriptive effects, explicit statuses and correction-only substitution.
- The source-only receipt retains IL17 HOLD and the unfulfilled-question/nonconfirmation limits.

## Validation and handoff boundary

Native environment: Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1. No packages were installed. The new test suite covers mathematical oracles, integer/coverage boundaries, zero and nonfinite input cases, exact constants, fractional df, complete correction inventories, shared resampling, every-unit omissions, source binding, changed pins and refusal before payload loading.

The source-only ten-program preflight is complete. ECM identifier coverage is 301/321 Core and 299/321 Strict. The closed nine-program mapping remains unchanged. Current compiled metadata digest, including the ECM ETS2-absence record: `17a291f47fbb0ec0e815ac2f592ba4d0f657f3f5e3cbe4689dae5ecc44c90698`.

Final checks passed:

| Check | Result |
|---|---|
| Native author suite | **109 tests** |
| Separate statistical/guard review | **108 synthetic tests**, seven boundary probes, source-alias and ancillary-reader probes |
| Independent verifier suite | **61 tests** |
| Full artificial pipeline comparison | **584,546 scalar comparisons** |
| Joint resampling checks | **640,000 RNG indices**, **30,000 constant draws** |
| Independent boundary / entry-guard / private-writer checks | **62 / 42 / 51** |
| Identifier/source compilation checks | **11,372**, including all 642 ECM mapping rows |
| Source-only replay | Four planning/preflight artifacts byte-identical; no versioned file changed by replay |

Both final reviews bind the exact script/test snapshot. The statistical/guard receipt is `work/psc-liver-program-plan/independent-review/final-review-receipt.json` (SHA-256 `3c454520d3f8936465831913aa3f115185cbde06055aae616f01839a5cdc58e0`). The independent full-pipeline/source acceptance is `work/psc-liver-program-independent-review/accepted-review.json` (SHA-256 `a3b12cf8735a99cdf82440a6ea52fa96b62a1c7d9b046d96a2ea7e41f7927055`). Both state `real_execution_reviewed=false`. These are computational checks, not biological replication or clinical review.

The task's native synthetic/metadata tests were run. The repository-wide suite was not rerun under this stage's no-real-expression boundary. Existing catalog import logic was not changed. No unrelated historical test failure was repaired.

Final file hashes, checks and preservation records are in `work/psc-liver-program-plan/completion-receipt.json`. **No real RNA payload was opened, normalized or scored. No real effect was computed.** The held IL17-response-transfer question remains unanswered. Root must accept the reviews, freeze the prospective plan and separately authorize real execution before any real normalization, score or effect.

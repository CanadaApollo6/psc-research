# Reproducing PRKD2 public-data recovery

The [report](../reports/prkd2-public-data-recovery.md) explains the result. This stage adds source recovery, target availability and marginal evidence-coverage diagnostics. It does not rerun or overwrite the earlier PRKD2 statistical models, and it emits no new H4 posterior.

Run commands from the repository root in the existing project environment. Python uses the already installed pandas, NumPy, SciPy, pysam and PyArrow. The independent numerical check uses the existing pinned R runtime; no package installation or model API call is required.

## Frozen scope and source records

- [Recovery plan](../config/prkd2-recovery-plan.json): six targets, 42 eligible Catalogue contexts, the unchanged full region and evidence/model gates. SHA-256 `8f212c4d24e87c62648f6802edafc36288f833c35eafd9d42efc9136b3aa708f`.
- [Follow-up plan](../config/prkd2-recovery-followup-plan.json): independently verified rs313839 control coordinate and six-position reference-panel audit.
- [Harmonised-copy plan](../config/prkd2-recovery-harmonised-plan.json): the complete PSC source archive, publisher MD5 and the union-of-build-intervals extraction rule. Build controls precede interpretation of the converted coordinates.
- [Author-source plan](../config/prkd2-recovery-author-plan.json): only the nine contexts whose Catalogue files returned unavailable-file responses. The four IBDverse blood-monocyte archive members are selected from pinned author annotations before new association inspection. The transfer limit remains 300 MB per context.
- [Initial source receipts](../config/prkd2-recovery-sources.json) and [additional receipts](../config/prkd2-recovery-extra-sources.json): URLs, dates, response headers, cached byte ranges, source failures and hashes. They distinguish unavailable files from measured variant absence. Metadata Git commits are pinned separately from changing web pages.

All six full-edit targets retain both builds and source alleles. Existing regional references and reciprocal chains come from the preceding [full-allele analysis](prkd2-reproduction.md). They must be restored from their pinned receipts to rerun sequence normalization. None of the source rows, earlier fit inputs, component vectors or posteriors is replaced by this stage. Session-cookie headers are omitted from the HTTP receipts; `headers_omitted` records their field names without retaining their values.

## Portable checks on the versioned results

The repository includes complete retained PRKD2 association vectors, the Catalogue point/regional extracts, target and identity tables, and the independent R coverage table. The following verifier checks these versioned outputs and recorded acquisition outcomes without requiring the large whole-genome archives:

```bash
.venv/bin/python scripts/verify_prkd2_recovery.py
.venv/bin/python -m unittest discover -s tests -v
```

The portable check explicitly reports how many ignored raw-cache receipts it did not inspect. It does not claim to reacquire the whole original studies or independently verify every reference base in a fresh clone. To recalculate the evidence-coverage numbers using R, independently of the Python implementation:

```bash
work/coloc-runtime/r/bin/Rscript --vanilla scripts/verify_prkd2_recovery_coverage.R data/derived/prkd2-recovery-independent-coverage.csv
.venv/bin/python scripts/verify_prkd2_recovery.py
```

The R calculation uses all original PSC marginal evidence, source-specific complete RNA beta/SE vectors where supplied, exact source-row deduplication and the same fixed effect priors. Missing DICE SE remains unavailable. It also independently repeats the distribution-convention check below and writes `prkd2-recovery-independent-format.csv` beside the coverage table. No posterior is calculated. The separate evidence-coverage percentages do not establish shared causality.

## Analysis replay with the original local cache

The analyses additionally need the existing reference/chain cache and, for the conditional-export audit, the three complete author text files. With that cache restored, these commands regenerate the small analysis tables from the retained source rows:

```bash
.venv/bin/python scripts/analyze_prkd2_recovery.py
.venv/bin/python scripts/audit_prkd2_recovery_secondary.py
.venv/bin/python scripts/analyze_prkd2_recovery_ibdverse.py
.venv/bin/python scripts/audit_prkd2_recovery_ibdverse_format.py
work/coloc-runtime/r/bin/Rscript --vanilla scripts/verify_prkd2_recovery_coverage.R data/derived/prkd2-recovery-independent-coverage.csv
.venv/bin/python scripts/verify_prkd2_recovery.py --with-cache
```

Use a separate work copy or compare the recorded output hashes before and after replay. The analysis scripts write deterministic CSV/gzip/JSON outputs and retain original source fields beside normalized identities. The final stage record reports the actual replay and prior-file preservation checks performed in this run.

The executed isolated replay copies scripts, configuration and compact inputs into a new ignored work directory, reads the original raw cache, runs all five calculation commands sequentially, and checks all 27 output files byte for byte:

```bash
.venv/bin/python scripts/replay_prkd2_recovery.py
```

## Complete source verification

The two original whole-genome source archives occupy about 4.31 GB compressed. DICE's uncompressed stream is about 27.17 GB, so verification must stream it rather than materializing a giant text file. The complete IBDverse selected members require about 562 MB compressed; the 393 GB parent ZIP was not downloaded. The harmonised PSC archive adds 324 MB. Byte-range caches and complete assembled archives can coexist locally and therefore use more disk space than those source-size totals.

```bash
.venv/bin/python scripts/verify_prkd2_recovery_archives.py
.venv/bin/python scripts/index_prkd2_recovery_ibdverse.py
.venv/bin/python scripts/verify_prkd2_recovery_ibdverse.py
```

The original-source verifier uses external `gzip`, `rg` and an independent field parser to rescan every DICE and BLUEPRINT row through EOF, comparing every retained byte. The IBDverse verifier uses Python's standard ZIP reader over cached intervals and a separate byte-field parser; the primary extractor uses explicit raw-DEFLATE blocks. Both require complete member decompression/CRC agreement and compare all selected PRKD2 and target-position rows. Central-directory replay checks all 6,047 member metadata records. The frozen inventory CSV retains its original CRLF record endings; whitespace checks recognize CRLF rather than rewriting the pinned bytes.

The [DICE receipt](../config/prkd2-recovery-dice-stream.json), [BLUEPRINT receipt](../config/prkd2-recovery-blueprint-stream.json), [BLUEPRINT transfer ledger](../config/prkd2-recovery-blueprint-range-acquisition.json), [IBDverse query ledger](../config/prkd2-recovery-ibdverse-queries.json) and [PSC harmonised ledger](../config/prkd2-recovery-harmonised-query.json) identify the exact bytes. Restoring only a JSON receipt cannot recreate an ignored raw input. Reacquisition in a fresh work copy must use the recorded URL/range and reject a differing hash, genome build or remote file identity. Changing public HTML and retired endpoints can prevent exact recovery of historical metadata snapshots; versioned selected metadata and source extracts remain the recorded analysis inputs.

## Acquisition entry points and execution details

1. `prkd2_recovery_sources.py` records bounded small-file, header, directory and index requests. The extra ledger uses the same reader with a separately assigned manifest/directory; concurrent writers never share a ledger.
2. `fetch_prkd2_recovery_dice.py` retrieves and scans the complete original classical-monocyte gzip archive. Neither a TBI nor a CSI index was supplied at the tested URLs. It retains all PRKD2 rows and all rs112445263 rows across genes, with no significance filtering.
3. `fetch_prkd2_recovery_blueprint_ranges.py` retrieves the original BLUEPRINT archive in exact 8-MB ranges using two workers and two-second request-start spacing. The whole-file download was initially slow and interrupted; its partial bytes/hash are preserved. Range assembly agrees with that original prefix and an independent source prefix. `fetch_prkd2_recovery_blueprint.py` reads the complete assembled file through EOF.
4. `fetch_prkd2_recovery_catalogue.py r7` queries the fixed SNP/indel positions across genes; `r8_beta` reads Parquet identifier columns and all candidate PRKD2 groups under the fixed row/byte caps. The 26 r7 and five successful r8 outcomes are complete. The nine recorded r8 HTTP failures are preserved as access outcomes. For historical verification, use the dedicated verifier rather than re-executing the failed network attempts into the saved ledgers.
5. `audit_prkd2_recovery_controls.py controls --offline` and `audit_prkd2_recovery_controls.py reference --offline` verify the completed true-control/site extracts. The reference stage retains VCF site fields only, not individual genotypes.
6. `fetch_prkd2_recovery_harmonised.py` retrieves the complete PSC copy in bounded ranges and verifies its publisher MD5 and gzip EOF. All union-region rows retain original and harmonised columns. The secondary audit also examines the 64 retained records without a harmonised variant ID, testing both source allele orientations against the verified reference; it does not quietly drop them.
7. `prefetch_prkd2_recovery_ibdverse.py` changes only transfer scheduling, using two workers with one manifest writer. `fetch_prkd2_recovery_ibdverse.py --offline` subsequently reads all four complete frozen members, checks full CRC/size/EOF and writes deterministic source extracts. The original single-worker transfer was interrupted after verified ranges had been cached; the original extractor and its hash remain unchanged.

## Preserved query correction and source limitations

The initial r7 reader queried the correct five missing-evidence positions but used **46,702,450** as its auxiliary control. The actual verified rs313839 coordinate is **46,718,300**. The initial reader remains hash-pinned as executed. A separately recorded correction queries the real control in **all 26** new r7 contexts; final target tables use those corrected control files. The primary missing-SNP queries, cohort selection and original data were unaffected. The unrequested auxiliary point is not used as an rs313839 result.

One exact duplicate BLUEPRINT association is preserved in the source extract and collapsed only in the evidence calculation. Two original DICE symbolic `<CN0>` records remain explicit identity exclusions. Out-of-window or unnormalized gene rows remain in the full RNA denominator where a valid beta/SE is supplied. DICE's rounded statistic is never converted into an invented SE. Original and reprocessed estimates remain separate.

Nassiri's inspected exports contain conditional forward/backward results without allele pairs. Natri provides public LIMIX data, but the unindexed 11.03-GB immune archive exceeds this stage's fixed transfer budget. These are data-type/resource boundaries, not biological nulls or claims that no public author data exist. The available IBDverse nominal exports do not themselves provide an original-study LD/covariance matrix or complete per-signal vectors.

The additional [distribution-convention plan](../config/prkd2-recovery-ibdverse-format-plan.json) was frozen after an extreme S100A8/9 beta/SE-to-P discrepancy appeared and before fitting a matching residual-df convention. It checks all 14,649 original PRKD2 rows across all four contexts against the same integer grid 1–500. R and Python independently infer 69/66/73/2 residual degrees of freedom from the numerical relationship; these are not observed participant counts. The normal-approximation coverage in the last context remains preserved but does not qualify model interpretation. No covariate count or original sample N is invented.

The [verification record](../reports/prkd2-recovery-verification.json), [cache record](../reports/prkd2-recovery-cache-verification.json), [original-archive check](../reports/prkd2-recovery-archive-verification.json), [IBDverse check](../reports/prkd2-recovery-ibdverse-verification.json) and [final stage record](../reports/prkd2-recovery-stage-verification.json) distinguish source/algorithm reproducibility from biological evidence.

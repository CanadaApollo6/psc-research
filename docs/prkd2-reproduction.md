# Reproducing the PRKD2 signal analysis

The [report](../reports/prkd2-signal-analysis.md) is the interpretation; CSVs, source extracts and JSON manifests are the computational record. Run commands from this repository's root. The private repository contains only public-source derivatives and research documentation for this stage. Raw caches, reference-donor genotypes, large reference assemblies, fitted R objects and work logs remain ignored.

## Plans, sources and software

- Original [plan](../config/prkd2-plan.json): frozen September 15, 2026 at 03:49:32 UTC, SHA-256 `6e3d484699d7f77c76eaede640b131f89f7d388a9961c1f209ebdaf4afaae51b`.
- [Full-allele extension](../config/prkd2-expanded-plan.json): frozen 04:11:06 UTC, SHA-256 `c169474ec10b9a99a91121085031627f0d0075bfb86ddf17f40f0fcd18520092`, after original baseline coverage failures and before expanded posteriors. It is post hoc and retains original results.
- [General source manifest](../config/prkd2-sources.json), [initial reused-source inventory](../config/prkd2-source-preparation.json), [GWAS requests](../config/prkd2-gwas-queries.json), [RNA requests](../config/prkd2-qtl-queries.json), and [reference requests](../config/prkd2-reference-queries.json) preserve URLs, byte ranges, dates, headers, stream sizes and hashes. Nominal indexed queries hash their full all-gene regional streams before PRKD2 filtering; those hashes are not whole remote-archive hashes.
- [BLUEPRINT full-component receipt](../config/prkd2-bf-query-QTD000021.json) and [DICE receipt](../config/prkd2-bf-query-QTD000504.json) document complete compressed streams through EOF, including 2,690,735,125 and 782,795,307 bytes. Full PRKD2 source vectors are versioned under `data/derived/prkd2-query-rows/`. Complete whole-genome BF archives were streamed, not retained. A fresh retrieval must check these whole-stream hashes before claiming the same source version.
- [Extended references](../config/prkd2-expanded-sources.json) include Ensembl GRCh37 chromosome 19, publisher CHECKSUMS, padded region sequences, and a second indexed public-genotype extraction retaining complete biallelic edits. GRCh38 reuses the previously checksum/decompression-verified Ensembl release 112 primary assembly. The exact earlier [reference record](../config/prkd2-prior-GRCh38-reference.json) is retained for reproduction; it does not contain the assembly itself.
- Use the existing pinned [R runtime](../config/coloc-runtime.json) and [explicit package lock](../requirements-coloc-linux-64.lock): R 4.4.3, coloc 5.2.3, susieR 0.14.2. Python uses the existing project environment, including pandas, NumPy, SciPy, pysam 0.23.3, bcftools 1.21 and Matplotlib. No statistical package upgrade was made for this analysis.

The eQTL Catalogue's [current access documentation](https://www.ebi.ac.uk/eqtl/Data_access/) directs users to downloadable/indexed files; this analysis does not depend on the retired REST API. The initial failed network/DNS source attempt is retained in the general manifest. GCST004030's previously pinned study record is reused; the live study API was not called again for this stage.

## Offline verification with the original local cache

The execution wrappers verify already-completed stages by hash; they do not silently overwrite a completed attempt. The separate replay utility creates a fresh work directory and reruns the statistical stages against the same prepared inputs. An original raw-cache file missing from a fresh clone must be restored or reacquired with the pinned source receipt before these cache-verification commands can pass.

```bash
.venv/bin/python scripts/verify_prkd2_inputs.py
.venv/bin/python scripts/verify_prkd2_statistics.py
.venv/bin/python scripts/verify_prkd2_normalization.py
work/coloc-runtime/r/bin/Rscript --vanilla tests/test_prkd2.R
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/replay_prkd2.py prkd2
.venv/bin/python scripts/replay_prkd2.py prkd2-expanded
.venv/bin/python scripts/plot_prkd2.py
```

The independent statistical verifier reads all versioned prepared/source-vector tables and checks the baseline formulas, full-component weights, source PIPs and complete planned grid. The raw-input verifier additionally needs the public-genotype subset and LD matrix, and checks all 4,727 source PSC rows, all 13,766 distinct RNA variants, effect signs, platform counts and the 4,513² correlations. The replay reruns the six fits rather than borrowing previous RDS objects; its input LD matrices still require the original cache or source-based reconstruction.

The normalization verifier uses an exact copy of the [normalization queue](../config/prkd2-normalization-checks.json), regional reference sequences and the prior full-GRCh38 reference receipt. If rebuilding only the work metadata in the same source cache:

```bash
mkdir -p work/prkd2 work/ubash3a-raw-read-audit
cp config/prkd2-source-preparation.json work/prkd2/source-preparation.json
cp config/prkd2-normalization-checks.json work/prkd2/expanded-normalization-checks.json
cp config/prkd2-prior-GRCh38-reference.json work/ubash3a-raw-read-audit/reference.json
```

These are metadata copies, not replacements for reference bases, associations or genotypes. Restoring a JSON receipt alone cannot recreate its input bytes.

## Acquisition and preparation order

The executed acquisition/preparation entry points, with their historical hashes, are preserved in the source and execution manifests:

1. `fetch_prkd2_sources.py` verifies/retrieves small sources, indexes, full credible-set files and compressed BF prefixes. It expects the historical reused-source inventory in `work/prkd2/source-preparation.json`; the versioned copy above makes this explicit. Earlier shared raw sources have their own original fetch/provenance records.
2. `fetch_prkd2_regions.py` acquires the complete GWAS envelope, both nominal RNA regions, and initial SNP reference genotypes. The stored indexed query extracts are sufficient to inspect the association rows without downloading whole nominal archives.
3. `fetch_prkd2_bfs.py` streams each entire ordinary-gzip component archive through EOF, retaining all PRKD2 vectors. Stopping when the first PRKD2 rows appear would not establish completeness. The separate per-dataset receipts prevent parallel transfers from overwriting one another's status.
4. `prepare_prkd2_inputs.py` performs original SNP matching, alias handling, full-evidence denominators, prior rs313839 reproduction and signed reference LD. `execute_prkd2.py baseline` and `execute_prkd2.py fit` generate the original statistical inputs/results.
5. `fetch_prkd2_expanded_sources.py` and `prepare_prkd2_expanded_inputs.py` reconstruct the separately frozen full-allele inputs and LD. `verify_prkd2_normalization.py` validates all queued normalized edits. `map_prkd2_bfs.py QTD000021` and `map_prkd2_bfs.py QTD000504` map complete published vectors to the same verified edit identities without changing their BF fields.
6. `execute_prkd2_expanded.py baseline`, `execute_prkd2_expanded.py fit`, **`execute_prkd2_signals.py compare`**, and `execute_prkd2_expanded.py compare` are the executed comparison paths. Each records complete input hashes before its own stage and verifies they are unchanged afterward.
7. `audit_prkd2_missing_evidence.py` quantifies omitted weights and records the two exact indexed missing-variant requests. `audit_prkd2_missing_alleles.py` tests all forty reported RNA omissions against every regional PSC source row. These audits do not alter analysis eligibility or results.

A clean-room acquisition should be performed in a new checkout/work copy with new retrieval/execution records, preserving these historical manifests rather than deleting or rewriting them. Public source contents and HTTP behavior can change. Refuse a byte/hash mismatch instead of claiming an exact reproduction. Large archived inputs may require the documented older source snapshot; the current files have not been assumed identical.

## Preserved convention corrections

**Explicit null column:** installed coloc 5.2.3 treats the last BF column as a null model during its overlap calculation. Published `lbf_variable` vectors contain only real variants. The [convention record](../config/prkd2-component-conventions.json), selected before any component comparison, therefore appends a named zero-log-BF null column when calling coloc, while real-variant overlap is calculated separately. `run_prkd2_signals.R` is the corrected original-scope entry point. The unexecuted comparison function in `run_prkd2.R` is retained because that file was already pinned by baseline/fit receipts; **do not use `execute_prkd2.py compare`**. The expanded entry point already includes the explicit null. Synthetic tests verify real last-column variants, reordering, a missing dominant variant and 89%/91% coverage.

**Input-verifier Z convention:** the first independent verifier correctly reproduced raw association inputs and all LD entries but initially compared raw signed Z with the finite-N-adjusted Z returned by `kriging_rss`. Reading the installed implementation confirms `z * sqrt((n-1)/(n-2+z*z))`. The [verifier correction](../config/prkd2-input-verifier-convention.json) records that failed assertion and its source hashes; the corrected check agrees within 5.4×10⁻¹⁴. Fit inputs and statistical results are unchanged.

**Complete-source omission audit:** a preliminary uncommitted audit restricted its candidate GWAS universe to non-SNV rows despite checking some RNA SNPs. The final scripted audit considers all 4,727 GWAS rows in both possible REF orientations and repeats all forty targets. Its [record](../reports/prkd2-missing-allele-verification.json) preserves the scope correction. It was not used to modify statistical inputs or select another cohort.

The frozen expanded baseline CSV retains legacy column names `eligible_gwas_snvs` and `eligible_qtl_snvs` from the original writer. In that expanded file they count all qualified variants, including non-SNV edits; the variant-class audit and report give the explicit breakdown. The source outputs remain byte-for-byte preserved.

The [statistical verification](../reports/prkd2-statistical-verification.json), [input verification](../reports/prkd2-input-verification.json), [normalization verification](../reports/prkd2-normalization-verification.json), both [original](../reports/prkd2-replay-verification.json) and [expanded replay](../reports/prkd2-expanded-replay-verification.json) records, and [final stage verification](../reports/prkd2-stage-verification.json) distinguish numerical reproducibility from biological inference. The final result remains coverage-limited.

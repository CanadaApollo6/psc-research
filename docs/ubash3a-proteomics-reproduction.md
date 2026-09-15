# Reproducing the public CD4 proteomics analysis

The [analysis report](../reports/ubash3a-proteomics-analysis.md) gives the scientific result. This document describes the stored evidence and the order required to reproduce direction 3. It does not extend the experiment or authorize interpreting a global false-discovery estimate as validation of a novel protein. Use a separate reproduction workspace for acquisition, conversion and new searches; use the replay command below to check the existing completed cache.

## Preserved inputs and environment

Use the repository root as the working directory. The [specificity plan](../config/ubash3a-proteomics-specificity-plan.json), [reference lock](../config/ubash3a-proteomics-reference-lock.json), [spectrum selection](../config/ubash3a-proteomics-spectrum-selection.json), [search plan](../config/ubash3a-proteomics-search-plan.json) and [Comet parameters](../config/ubash3a-proteomics-comet.params) define the completed scope. The fourteen raw files have a combined source size of 2,387,832,397 bytes. Their archive metadata and computed checksums are preserved separately from the converted spectra.

Large downloads, full search outputs and vendor binaries are ignored local files. A clone contains small results, protocols, source manifests and verification code; it does not contain a complete executable data cache. The exact local snapshot is needed for byte-identical offline replay. The UniProt service is live: a later release can legitimately return different sequences, and must not be silently substituted for the recorded 2026_03 snapshot. Reacquisition is a new provenance event, not recreation of the original retrieval timestamps.

The [observed Python dependencies](../requirements-ubash3a-proteomics.txt) were already installed in the local research environment. The [environment record](../config/ubash3a-proteomics-reproduction-support.json) records their versions and the Python version. The [tool installation record](../config/ubash3a-proteomics-tool-install.json) pins the publisher asset URLs, sizes and hashes for Comet 2026.02.2 and ThermoRawFileParser 2.0.0-dev. Its installation-time status is historical; later qualification reports establish the actual runtime checks. Binaries are not committed or redistributed. Review the publisher/vendor terms when obtaining them for another environment.

## Source acquisition

The exact service wrappers used for this run are retained in [the source-helper directory](../scripts/ubash3a_proteomics_sources/fetch_skill_source.py). They call the unchanged PRIDE/UniProt skill clients, adding response-header receipts. Their installed-client paths reflect this host; a different installation must resolve and record its own client path and hash. [Request templates](../config/ubash3a-proteomics-source-requests/uniprot-human-page000.request.json) differ from the original requests only by making output paths repository-relative.

In a separate reproduction workspace, copy those templates into `work/ubash3a-proteomics-specificity/`. Fetch the first UniProt page with `scripts/ubash3a_proteomics_sources/fetch_skill_source.py`, then run `scripts/ubash3a_proteomics_sources/fetch_uniprot_pages.py`. The second helper follows the returned next links, checks every release/total header, rejects repeated pages, checks the terminal page and concatenates the saved FASTA pages. It intentionally requires the historical release and count. The recorded all-at-once export failed and is retained as a failure; it is not a successful source.

The complete GENCODE 50 translation file is checked against its publisher MD5 and local SHA-256. The original cRAP FTP file is preserved with source size, modification time and local SHA-256. All twelve deposited pepXML files were read only far enough to capture their complete search-setting elements. Their original peptide identifications are not new evidence from the reanalysis. The [source manifest](../config/ubash3a-proteomics-sources.json) retains the URLs, range responses, successful receipts and failed retrievals. Earlier conditional protein models and the complete PRIDE file inventory remain hash-pinned inputs from prior stages.

`scripts/fetch_ubash3a_proteomics_spectra.py` obtains only the fourteen frozen raw files, checks every expected byte total and computes SHA-256. No publisher raw-file checksum was supplied. The downloader saves individual receipts and a complete manifest. It must not be pointed at a different experiment while retaining the existing selection record.

## Reference checks and search database

The primary and independent specificity implementations are:

```text
.venv/bin/python scripts/analyze_ubash3a_proteomics_specificity.py
.venv/bin/python scripts/verify_ubash3a_proteomics_specificity.py
```

The first scans full sequences with literal and I/L-equivalent matching, then reports digestion context without discarding nontryptic competitors. The second uses Biopython, overlapping regular expressions, a separate cleavage-gap calculation and independent XML parsing. The four targets are fixed by the prior conditional translation; they are not learned from these spectra.

`scripts/prepare_ubash3a_proteomics_references.py database` builds exact-sequence groups with all source memberships and checks candidate decoy collisions. The original database FASTA and membership stream are pinned in the database record. Rebuilding the gzip membership stream can change timestamp metadata despite identical content; the independent search verifier compares all source records, headers, sequences and memberships. Do not regenerate a reference lock or replace the database in the original completed run merely to refresh a timestamp.

## Conversion and search qualification

The original Thermo spectrum converter is used with MS2-only, vendor-centroided mzML output. Its development-build metadata writer is disabled because of the documented repeated scan-count loop and mislabeled term. The [separate metadata reader](../scripts/ProteomicsRawMetadata/Program.cs) uses one pass through the vendor API to record all MS orders and every MS2 scan ID. Its source, dependencies and compiled artifact are pinned in the [repair record](../config/ubash3a-proteomics-metadata-repair.json).

The metadata reader was built with the existing .NET SDK 10.0.401 and the vendor assemblies from the pinned parser installation, without downloading new package dependencies. To rebuild in a separate workspace, keep intermediate and compiled output under `work/`:

```text
dotnet build scripts/ProteomicsRawMetadata/ProteomicsRawMetadata.csproj -c Release -o work/tools/thermo-metadata-validator -p:BaseIntermediateOutputPath=../../work/ubash3a-proteomics-spectra/dotnet-obj/
```

A rebuild is a new binary provenance event. Record its runtime and artifact hashes before conversion; do not overwrite the original repair record to make a different binary appear historical. Both conversion and metadata reading rely on the vendor library. Exact scan agreement is a software-consistency check, not an independent instrument measurement.

```text
.venv/bin/python scripts/qualify_ubash3a_comet.py
.venv/bin/python scripts/qualify_ubash3a_spectra.py
.venv/bin/python scripts/run_ubash3a_proteomics_search.py
```

The Comet qualification uses seven generated controls plus an excluded precursor; none is a T29 or public spectrum. Conversion checks the complete mzML scan list, finite binary arrays, centroiding, retention times, precursor values, charges and activation/analyzer annotations against the metadata export. Per-file raw, converted-output and log hashes precede scoring.

The completed run overlapped downloading, conversion and search under the [prospective execution amendment](../config/ubash3a-proteomics-execution-order.json). Each whole file had to pass before scoring. The same fourteen files, database and parameters were retained. No candidate-score inspection or pooled error calculation began until all fourteen searches finished. A sequential replay may use the simpler command order above.

## Analysis and independent verification

After all fourteen files complete:

```text
.venv/bin/python scripts/analyze_ubash3a_proteomics_spectra.py
.venv/bin/python scripts/review_ubash3a_candidate_fragments.py
.venv/bin/python scripts/verify_ubash3a_proteomics_search.py
.venv/bin/python scripts/audit_ubash3a_proteomics_scan_coverage.py
.venv/bin/python scripts/build_ubash3a_proteomics_source_manifest.py
.venv/bin/python -m unittest discover -s tests -v
```

The analysis retains one winner per original file/scan, resolves reported target–decoy ties conservatively and estimates spectrum and unique I/L-equivalent peptide q-values with the frozen `(D + 1) / T` rule. Converted scans, engine-loaded queries and scans with reported winners are counted separately. All four candidates are inspected among all five exported assignment positions, irrespective of global q-value. XCorr ranks are dense and can repeat; export position is retained separately.

The [counting-scope record](../config/ubash3a-proteomics-counting-scope.json) pins the source code showing that Comet's loaded count is the query-vector size. The configured 5,000 batch threshold can yield slightly more loaded queries when active preprocessing workers finish; the observed 5,005-query batch does not represent a changed parameter or a truncated source file. The scan-coverage audit accounts for every converted scan ID, all exported XML queries and scans without a reported assignment. Missing output is not automatically classified as an unsearched scan, because the export cannot distinguish every filtering/no-competitor case.

The fragment review uses independent monoisotopic masses and the original centroid peaks with a fixed 0.02-Da matching window for base b/y ions. It retains the best exported noncandidate target for the same query where available. It does not reproduce Comet's complete intensity/neutral-loss score and cannot itself calibrate novel-peptide false discovery.

The independent verifier checks all 615,403 source records and the 308,191 distinct protein sequences, parses complete pepXML outputs, compares every exported assignment with the text output, and recalculates both q-value levels using separate tabular operations. XML does not retain decoy prefixes on alternative protein names; the independently parsed PIN output supplies those labels, with full membership multisets, modifications and mass/score fields checked. The PIN row's single target/decoy label is insufficient for mixed membership, so the existing conservative rule uses every protein membership. This is computational verification of the same experiment, not independent donor replication.

The [first parser repair](../config/ubash3a-proteomics-output-parser-repair.json) documents repeated rank one and the distinction between exported Comet scan ordinals and native instrument IDs. The [verifier format repair](../config/ubash3a-proteomics-verifier-format-repair.json) preserves the failed parent-label assumption and two additional publisher source files. Original failures, code hashes and primary output hashes remain available. Reproduction must use the repaired parsers; neither repair changes the searched files, database, parameter file or error threshold.

The small specificity and primary spectrum tables are deterministic from their pinned inputs. The complete spectrum/peptide winner tables stay under ignored `work/`; their hashes are recorded in the summary. Replay the completed local cache with:

```text
.venv/bin/python scripts/replay_ubash3a_proteomics.py
```

This writes 19 specificity/metadata files and five spectrum-summary files to a fresh ignored directory and compares their bytes with the originals. It also rebuilds the two large winner tables and requires their hashes to remain identical. An existing replay directory is preserved; use `--output work/another-replay-directory` for another verification. Database construction, spectrum searches and fragment review are checked separately rather than rerun by this command.

Verification and execution records preserve the distinction between a qualified method, an observed assignment and a biological interpretation. A missing exported target is not an exhaustive targeted spectrum test and does not prove protein absence.

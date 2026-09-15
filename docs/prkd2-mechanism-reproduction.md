# Reproduce the PRKD2 regulatory follow-up

This guide reproduces the completed September 15, 2026 [mechanism report](../reports/prkd2-regulatory-mechanism.md). The scope is the newly recovered rs112445263 edit, previous rs313839 control, one qualifying comparison, and a fixed ENCODE panel. It is separate from earlier PRKD2 shared-signal calculations and the original-source recovery.

## Environment and immutable inputs

Use the existing project `.venv/bin/python`. The recorded client is AlphaGenome 0.9.0 at commit `aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d`. Runtime dependencies include numpy, pandas, scipy, requests, pyarrow, openpyxl and matplotlib. The model itself is the remote ALL_FOLDS service; the response does not disclose its build identifier.

The [plan](../config/prkd2-mechanism-plan.json), [baseline freeze](../config/prkd2-mechanism-freeze.json), [input lock](../config/prkd2-mechanism-inputs.json), [ENCODE selection](../config/prkd2-mechanism-encode-selection.json), [computed-audit qualification](../config/prkd2-mechanism-encode-audit-resolution.json) and [execution manifest](../config/prkd2-mechanism-execution.json) are immutable records of this run. Do not edit them to make a future result pass. The comparison shortfalls, 50-gene intersection and missing monocyte ATAC track remain part of the result.

The model output directory is `data/predictions/20260915T185853Z-prkd2-mechanism/`. It contains three tidy score files, 36 original component matrices with their axes, and REF/ALT arrays for both anchors. The versioned [model-run manifest](../config/prkd2-mechanism-model-run.json) records output identities, time stamps, shapes and hashes. Numeric arrays are not in Git.

## Replay existing local inputs

From the repository root:

```bash
.venv/bin/python scripts/summarize_prkd2_mechanism.py
.venv/bin/python scripts/audit_prkd2_mechanism_training.py --offline
.venv/bin/python scripts/analyze_prkd2_mechanism_peaks.py --offline
MPLCONFIGDIR=/tmp/psc-mechanism-mpl .venv/bin/python scripts/plot_prkd2_mechanism.py
.venv/bin/python scripts/verify_prkd2_mechanism.py
.venv/bin/python scripts/replay_prkd2_mechanism.py
.venv/bin/python scripts/build_catalog.py
.venv/bin/python -m unittest discover -s tests -v
```

The first three commands deterministically regenerate tables from pinned inputs. The replay command runs them again in an isolated directory under `work/`, compares 17 output files to the recorded results and fails on a difference. No replay command calls the model or network. The verifier independently reconstructs all component-to-table values, RNA ranks, original BED overlaps and public reference LD; it also checks source DNA and preservation of earlier scientific files.

Figure PDFs may include renderer-generated metadata and are not claimed to reproduce byte-for-byte. Their source tables, fixed crops, PNG/PDF hashes and plotting code are pinned in [figure provenance](../reports/figures/prkd2-mechanism/figure-provenance.json). Figures show both anchors' fixed ±2,500-bp windows regardless of effect size, and include all prespecified monocyte readouts.

After documentation is finalized, run:

```bash
.venv/bin/python scripts/verify_prkd2_mechanism_stage.py
```

That stage check validates all new cached source files and request receipts, the versioned model manifest, deterministic-result manifests, local document links and prior scientific-file preservation. It writes [the completion record](../reports/prkd2-mechanism-stage-verification.json).

## Cache and source boundaries

A fresh clone cannot replay ignored caches by itself. Restore exact files and check hashes using the [cache manifest](../config/prkd2-mechanism-cache-manifest.json), [reference-source ledger](../config/prkd2-mechanism-reference-sources.json), [peak receipts](../data/derived/prkd2-mechanism-encode/peak-analysis.json), [training audit](../data/derived/prkd2-mechanism-training/training-audit.json) and [public-LD receipt](../data/derived/prkd2-mechanism-results/reference-LD-provenance.json). Earlier cached inputs include:

- `data/raw/prkd2-GRCh38-region.txt`, the independently checked regional Ensembl DNA;
- `data/raw/GRCh38.p13.genome.fa.fai` and `data/raw/match-model-gencode-v46.feather`;
- `data/raw/qtl-alphagenome-supplement-tables.data`, the publisher XLSX training workbook;
- `data/raw/prkd2-expanded-EUR-dosages.tsv.gz`, exact public phase-3 EUR dosage rows, ignored because they contain sample-level data;
- the prior complete regional GWAS/alignment derived inputs used to select comparisons.

The preparation scripts are `prepare_prkd2_mechanism_sources.py`, `prepare_prkd2_mechanism_inputs.py`, `prepare_prkd2_mechanism_encode.py` and `qualify_prkd2_mechanism_encode.py`. They preserve requests, immutable locks and bounded download receipts. All ENCODE metadata requests use the installed `life-science-research:encode-skill` REST helper. Its location and checksum are recorded; a machine without that installation must restore the helper before new ENCODE requests. Source metadata JSON is parsed and reserialized, not raw HTTP bytes. Peak files are exact downloaded gzip bytes and match publisher MD5s.

Do not rerun `run_prkd2_mechanism.py --run` merely to replay a completed analysis. Successful predictions are already cached. A new model execution needs a separately dated scope and output directory and may differ even with the same client. The existing runner's default dry run only validates preparation and checks whether a key is configured. Credentials are read privately through the existing loader; never put them in arguments, Git or output.

## Interpretation boundaries

Training experiment/file membership is checked across all published human rows. The donor/biosample extension is intentionally bounded to monocyte RNA, DNase, ATAC, CAGE and histone training sources, using all replicates of those experiments. It does not certify independence of an unmatched donor or the current serving model. Source peak overlaps are regulatory context, not allele-specific experiments or verified PRKD2 links.

The summary retains one unmatched anchor, one singly matched anchor, missing splice outputs and missing assay contexts. Reference LD is descriptive and uses full GRCh37 allele identity where source rsIDs are absent. Neither reusing that reference nor aligning effect signs supplies study-specific covariance, a new shared-cause posterior or evidence for treatment effects.

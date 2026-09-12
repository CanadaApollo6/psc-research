# PSC research

Independent computational research into **primary sclerosing cholangitis (PSC)** using public genetics, gene-expression data, and regulatory-sequence models.

## First research question

**Can AlphaGenome help connect PSC-associated DNA variants to the genes and cell types they affect?**

Start with three published examples—**UBASH3A, ETS2, and PRKD2**—to assess the method, then investigate less-resolved regions. The intended result is a reproducible set of biological hypotheses with evidence and uncertainty, rather than a treatment recommendation. Public-data analysis can proceed independently; researcher feedback can improve it without being a prerequisite.

The starting reference is [Goode et al., Nature Communications (2024)](https://www.nature.com/articles/s41467-024-53602-w). Its Table 1 summarizes 19 association signals across 15 regions. Our import preserves the published values, including discrepancies documented in the [evidence log](docs/evidence-log.md).

## First results

The latest [direct-junction follow-up](reports/ubash3a-junction-followup.md) checks the UBASH3A strand discrepancy and measures public GEUVADIS junction reads. **The positive-strand boundaries are supported, and a published workflow-setting mismatch could explain the reversed labels. The actual historical run settings remain unconfirmed.** The public-data comparison was too sparse: only 18 of 360 genotyped European donors met the fixed coverage threshold, with no CC donors retained. An exact [data request](docs/ubash3a-data-request-draft.md) is drafted and unsent.

The earlier [measured-RNA follow-up](reports/qtl-followup.md) compares the model with public genotype–RNA associations. **PRKD2's expression direction agrees in BLUEPRINT and a separate DICE cohort. UBASH3A's total-RNA direction disagrees in both resting and activated DICE CD4 cells.** A 471-sample blood dataset contains the predicted +29-nucleotide splice boundary and direction, subject to the strand-provenance question above.

The [first AlphaGenome pilot](reports/first-alphagenome-pilot.md) and [mechanism audit](reports/mechanism-audit.md) remain available unchanged apart from forward links. Original GWAS rows support **C as the working PRKD2 rs313839 risk allele**. The subsequent measured coefficients clarify the molecular direction without silently correcting conflicting labels in the original literature. These are molecular research checks, not treatment findings.

## What works now

- A checksum-verified downloader for six small public source files.
- A reproducible import of Table 1's **19 signal summaries**, molecular colocalisation evidence, and **45 public GEO sequencing-library metadata records**.
- A dated [study and company watchlist](docs/watchlist.md) and [research contact directory](docs/contacts.md).
- A [research protocol and backlog](docs/research-plan.md), including the next AlphaGenome steps.
- Four verified build-38 variant inputs, a fixed scoring protocol, completed model requests and a report of all outcomes.
- A source-backed allele/transcript audit, one additional splice prediction, a documented coordinate correction, and [concrete researcher review questions](docs/mechanism-review-questions.md).
- Eight reproducible public QTL queries, four expression comparisons, two splice associations, a dataset-wide strand audit and explicit accounting of missing measurements.
- A reference-motif and workflow-strand check, public single-variant genotype matching, and a donor-level two-junction analysis that preserves its inconclusive coverage result.

The full signal-summary table contains the most probable variant per signal, not every credible-set member. Its original coordinates remain in build 37. The separate [four-variant pilot input](data/derived/benchmark-variants.csv) has verified build-38 coordinates and reference/alternate alleles. Its two original direction exclusions remain unchanged; the later [allele audit](data/derived/allele-audit.csv) is separate. No liver-atlas expression matrices have been analyzed.

## Reproduce the starting data

Use Python 3.10 or newer. The download is under 1 MB in the initial snapshot; the large sequencing matrices are not part of this setup.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/fetch_sources.py
python scripts/build_catalog.py
python -m unittest discover -s tests -v
```

The downloader caches files in ignored `data/raw/` and fails if a pinned checksum changes. To verify an existing cache without accessing the network:

```bash
python scripts/fetch_sources.py --offline
```

Small derived tables are versioned in `data/derived/`; the import summary is in [reports/data-audit.json](reports/data-audit.json). See [data provenance](docs/data-sources.md) for attribution and source terms.

## AlphaGenome approach

Prefer the hosted API and precomputed Atlas for the pilot. Local GPU hardware is not required to call the hosted service. The [official client](https://github.com/google-deepmind/alphagenome) documents API-key access, noncommercial use, output restrictions, and model capabilities. [Atlas](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/) provides precomputed single-nucleotide variant predictions.

The first run used the hosted scoring API with four requests; Atlas remains an option for larger future variant sets. Genome builds, alleles and cell-type availability were checked, and the [pilot protocol](config/alphagenome-pilot.json) was fixed before predictions. Model output remains subject to its own terms; it is not automatically covered by any future repository license.

Use Python 3.12 for the recorded model environment. Install the pinned dependencies in your virtual environment and reproduce the verified inputs:

```bash
python -m pip install -r requirements-alphagenome.lock
python scripts/fetch_sources.py --manifest config/reference-sources.json
python scripts/normalize_benchmark.py
```

Set `ALPHAGENOME_API_KEY` locally, or store only the key in the ignored file `data/private/alphagenome_api_key`. Do not put a key in a command argument, Git, issue, or chat. The API client reads it without printing it.

```bash
python scripts/inspect_alphagenome.py
python scripts/run_alphagenome_pilot.py        # validates inputs; no predictions
python scripts/run_alphagenome_pilot.py --run  # makes four prediction requests
python scripts/summarize_pilot.py data/predictions/20260911T223028Z
```

For a new run, replace the final path with the directory printed by the runner. Full predictions stay in ignored local storage; compact reports and provenance are versioned. The checked-in report describes the September 11 run. A changed server model can produce different future results.

## Reproduce the mechanism follow-up

After installing the model environment above, retrieve the additional pinned public sources and regenerate the allele audit:

```bash
python scripts/fetch_sources.py --manifest config/mechanism-sources.json
python scripts/fetch_sources.py --manifest config/association-sources.json
python scripts/fetch_sources.py --manifest config/splice-audit-sources.json
python scripts/fetch_gwas_audit_rows.py  # uses cached ranges when present
python scripts/build_allele_audit.py
python scripts/run_splice_followup.py   # dry-run with corrected v0.2 coordinate
python scripts/summarize_splice_followup.py data/predictions/20260911T225056Z-ubash3a-splice
python -m unittest discover -s tests -v
```

The summary command requires the saved local prediction folder. To obtain new predictions, use `python scripts/run_splice_followup.py --run` and summarize the printed folder. The original run used v0.1; its one-base endpoint mistake is documented alongside the correction, with the original protocol archived. Version 0.2 changes the interpretation coordinate, not the DNA prediction request. The report distinguishes the recorded results from anything produced by a future model update.

## Reproduce the measured-RNA follow-up

This analysis uses public summary statistics and needs no API key or GPU. Use Python 3.12 in an isolated environment:

```bash
python -m pip install -r requirements-qtl.txt
python scripts/fetch_sources.py --manifest config/qtl-followup-sources.json
python scripts/fetch_qtl_followup.py --offline
python scripts/summarize_qtl_followup.py
python -m unittest discover -s tests -v
```

The eight small, recorded association subsets are versioned. `--offline` verifies their hashes; the summary command also needs the pinned raw metadata downloaded in the preceding step. To repeat the remote indexed queries, run `python scripts/fetch_qtl_followup.py --refresh`. It requires exact agreement and preserves the recorded result if the archive changes. The raw multi-gigabyte association files are never downloaded in full. The Catalogue REST API has been retired; this reader uses its documented indexed archive.

## Reproduce the direct-junction follow-up

After installing `requirements-qtl.txt`, run:

```bash
python scripts/fetch_sources.py --manifest config/ubash3a-junction-sources.json
python scripts/fetch_ubash3a_junction_genotypes.py
python scripts/audit_ubash3a_strand.py
python scripts/summarize_ubash3a_junctions.py
python -m unittest discover -s tests -v
```

The [plan](config/ubash3a-public-junction-plan.json) was saved before count/genotype-group inspection; its checksum and timestamp are preserved in [run metadata](config/ubash3a-junction-run.json). The source manifest pins 35 files, including a 53 MB coverage-plot archive. Use `--offline` with the two fetch scripts to verify the existing cache; the genotype reader also supports `--refresh`. Sample-level public data stay ignored, and the aggregate report includes missing measurements and unestimable models.

## Interpretation

PSC susceptibility is different from progression of established PSC. A regulatory prediction can support a mechanism worth testing; it does not show that altering that gene will safely treat disease. Experimental and clinical validation would still be needed.

This repository contains public-source research material and unsent draft questions. Personal health records and actual private correspondence are kept outside it. The watchlist is a dated snapshot and does not refresh automatically.

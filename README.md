# PSC research

Independent computational research into **primary sclerosing cholangitis (PSC)** using public genetics, gene-expression data, and regulatory-sequence models.

## First research question

**Can AlphaGenome help connect PSC-associated DNA variants to the genes and cell types they affect?**

Start with three published examples—**UBASH3A, ETS2, and PRKD2**—to assess the method, then investigate less-resolved regions. The intended result is a reproducible set of biological hypotheses with evidence and uncertainty, rather than a treatment recommendation. Public-data analysis can proceed independently; researcher feedback can improve it without being a prerequisite.

The starting reference is [Goode et al., Nature Communications (2024)](https://www.nature.com/articles/s41467-024-53602-w). Its Table 1 summarizes 19 association signals across 15 regions. Our import preserves the published values, including discrepancies documented in the [evidence log](docs/evidence-log.md).

## First results

The [first AlphaGenome pilot is complete](reports/first-alphagenome-pilot.md): four variants across UBASH3A, ETS2 and PRKD2. The model ranked UBASH3A and PRKD2 first in the selected contexts, but UBASH3A's expression direction opposed the published comparison. ETS2 results differed between the lead and highest-posterior variants. These are mixed consistency results on known examples, not new treatment findings.

## What works now

- A checksum-verified downloader for six small public source files.
- A reproducible import of Table 1's **19 signal summaries**, molecular colocalisation evidence, and **45 public GEO sequencing-library metadata records**.
- A dated [study and company watchlist](docs/watchlist.md) and [research contact directory](docs/contacts.md).
- A [research protocol and backlog](docs/research-plan.md), including the next AlphaGenome steps.
- Four verified build-38 variant inputs, a fixed scoring protocol, completed model requests and a report of all outcomes.

The full signal-summary table contains the most probable variant per signal, not every credible-set member. Its original coordinates remain in build 37. The separate [four-variant pilot input](data/derived/benchmark-variants.csv) has verified build-38 coordinates and reference/alternate alleles. Two variants still lack a usable source-based risk-direction comparison. No liver-atlas expression matrices have been analyzed.

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

## Interpretation

PSC susceptibility is different from progression of established PSC. A regulatory prediction can support a mechanism worth testing; it does not show that altering that gene will safely treat disease. Experimental and clinical validation would still be needed.

This repository contains public-source research material. Personal health records and correspondence are kept outside it. The watchlist is a dated snapshot and does not refresh automatically.

# PSC research

Independent computational research into **primary sclerosing cholangitis (PSC)** using public genetics, gene-expression data, and regulatory-sequence models.

## First research question

**Can AlphaGenome help connect PSC-associated DNA variants to the genes and cell types they affect?**

Start with three published examples—**UBASH3A, ETS2, and PRKD2**—to assess the method, then investigate less-resolved regions. The intended result is a reproducible set of biological hypotheses with evidence and uncertainty, rather than a treatment recommendation. Public-data analysis can proceed independently; researcher feedback can improve it without being a prerequisite.

The starting reference is [Goode et al., Nature Communications (2024)](https://www.nature.com/articles/s41467-024-53602-w). Its Table 1 summarizes 19 association signals across 15 regions. Our import preserves the published values, including discrepancies documented in the [evidence log](docs/evidence-log.md).

## What works now

- A checksum-verified downloader for six small public source files.
- A reproducible import of Table 1's **19 signal summaries**, molecular colocalisation evidence, and **45 public GEO sequencing-library metadata records**.
- A dated [study and company watchlist](docs/watchlist.md) and [research contact directory](docs/contacts.md).
- A [research protocol and backlog](docs/research-plan.md), including the next AlphaGenome steps.

**No AlphaGenome predictions or expression-matrix analyses have been run.** The signal table contains the most probable variant per signal, not every member of each credible set. It is not yet a model-ready variant file: its coordinates use build 37, and reference/alternate alleles still need verification.

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

We will first verify genome builds, alleles, and available cell-type tracks, then freeze the benchmark before querying predictions. Keep API keys outside Git. Model output remains subject to its own terms; it is not automatically covered by any future repository license.

## Interpretation

PSC susceptibility is different from progression of established PSC. A regulatory prediction can support a mechanism worth testing; it does not show that altering that gene will safely treat disease. Experimental and clinical validation would still be needed.

This repository contains public-source research material. Personal health records and correspondence are kept outside it. The watchlist is a dated snapshot and does not refresh automatically.

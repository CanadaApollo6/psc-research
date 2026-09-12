# PSC research

Independent computational research into **primary sclerosing cholangitis (PSC)** using public genetics, gene-expression data, and regulatory-sequence models.

## First research question

**Can AlphaGenome help connect PSC-associated DNA variants to the genes and cell types they affect?**

Start with three published examples—**UBASH3A, ETS2, and PRKD2**—to assess the method, then investigate less-resolved regions. The intended result is a reproducible set of biological hypotheses with evidence and uncertainty, rather than a treatment recommendation. Public-data analysis can proceed independently; researcher feedback can improve it without being a prerequisite.

The starting reference is [Goode et al., Nature Communications (2024)](https://www.nature.com/articles/s41467-024-53602-w). Its Table 1 summarizes 19 association signals across 15 regions. Our import preserves the published values, including discrepancies documented in the [evidence log](docs/evidence-log.md).

## First results

The latest [liver-cell context analysis](reports/liver-atlas-cell-context.md) checks **195,417 public cells and nuclei**, retaining 181,125 under fixed annotation rules. Donor-level measurements support **UBASH3A in CD4 T cells and ETS2/PRKD2 in monocytes** as contexts for the existing hypotheses. All five candidates, including PFKFB3 and BCL2L11, have usable measurements in the selected PSC populations. Whole-cell controls use a different chemistry, and the relevant nuclear immune-cell groups have too few eligible controls for our disease comparison. This establishes measurement context, not a PSC-specific change or treatment effect. A [public ETS2 benchmark inventory](docs/ets2-benchmark-next-inputs.md) prepares the next analysis of a published downstream program; no program score has been calculated yet.

The preceding [variant-coverage follow-up](reports/variant-coverage-repair.md) recovers **two exact PFKFB3 deletions** and finds **all four missing BCL2L11 variants in three larger OneK1K CD4 groups**. PFKFB3's strongest RNA components still fail the fixed coverage gate (52.23% and 88.55%). The four BCL2L11 positions are absent across all genes in all nine selected DICE archives; they were already measured in BLUEPRINT. In OneK1K, all twelve target associations are weak (P=0.117–0.952), and regional RNA coverage remains inadequate. Neither gene has supported evidence of a shared PSC–RNA mechanism. Exact source gaps, filtering limits, sample-count discrepancies and [unsent source questions](docs/coverage-repair-source-questions.md) are recorded.

The preceding [colocalisation assessment](reports/colocalisation.md) checks **PFKFB3 and BCL2L11 across all ten original CD4 contexts**. PFKFB3's eligible-SNP and secondary RNA-signal results favor different signals. Its original coverage failures and the source-supported platform-count amendment remain preserved; the new coverage audit does not overwrite them or report a new colocalisation posterior.

The preceding [CD4 RNA evidence report](reports/cd4-rna-evidence.md) checks **all nine variants across ten datasets and 228 fixed gene pairs**. Of 2,280 planned dataset-level comparisons, 693 are measured. Five associations pass the fixed multiple-testing correction, representing **rs7923054–PFKFB3** and **rs72837826–BCL2L11**. Model RNA tracks disagree on direction for both pairs. DICE subsets share donors, and both DICE and BLUEPRINT appear in the original PSC study, so this is not wholly independent validation. These observed RNA associations are preserved alongside the subsequent shared-signal assessment.

The preceding [fixed matched comparison](reports/matched-comparison.md) completed **nine predictions across three matched groups**. One prioritized variant scored above its comparison pair and two below it; there is no consistent separation in this small sample. Both BACH2 candidates lacked an eligible pair under the unchanged matching rules and were excluded before scoring. All nine variants have complete model outputs across their fixed gene and CD4-track universes.

The preceding [fine-mapping archive audit](reports/psc-finemap-and-comparison.md) imports **71,083 variant rows across 18 GWAS regions and eight molecular-QTL datasets**. SNP probabilities agree with saved configurations, but accompanying logs have substantial provenance discrepancies. Exact published signal-specific credible sets remain unavailable. These limits carry forward into the model comparison.

The [direct-junction follow-up](reports/ubash3a-junction-followup.md) checks the UBASH3A strand discrepancy and measures public GEUVADIS junction reads. **The positive-strand boundaries are supported, and a published workflow-setting mismatch could explain the reversed labels. The actual historical run settings remain unconfirmed.** The public-data comparison was too sparse: only 18 of 360 genotyped European donors met the fixed coverage threshold, with no CC donors retained. An exact [data request](docs/ubash3a-data-request-draft.md) is drafted and unsent.

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
- A full fine-mapping archive import, log/configuration audit, explicitly labeled singleton-set reconstructions, and a fixed prospective comparison protocol with candidate references and exclusions.
- Complete matching audits for 6,117 low-PIP rsIDs, frequencies and LD from 503 public EUR donors, nine successful predictions, and all 228 variant–gene rankings with 456 primary track scores.
- Ninety bounded public CD4 eQTL queries, a complete 2,280-row coverage/effect table, corrected-threshold associations, and explicit model-direction, cell-label and donor-sharing limitations.
- Twenty complete regional PSC/RNA comparisons, audited genotyping-platform sample counts, original and amended ABF baselines, multiple-signal sensitivity analyses, full published target-gene BF vectors and explicit signal-weight coverage failures.
- Exact PFKFB3 deletion normalization, an all-gene DICE coverage audit and twelve reannotated OneK1K target measurements, with unresolved genetic coverage retained.
- Two checksum-verified liver-atlas expression objects, public donor metadata, five-gene donor-level summaries at two cell-count thresholds, complete comparison accounting and a separate ETS2 experimental-input inventory.

The full signal-summary table contains the most probable variant per signal, not every credible-set member. Its original coordinates remain in build 37. The separate [four-variant pilot input](data/derived/benchmark-variants.csv) has verified build-38 coordinates and reference/alternate alleles. Its two original direction exclusions remain unchanged; the later [allele audit](data/derived/allele-audit.csv) is separate. The liver-atlas analysis provides observational cell context and does not measure these donors' genotypes or splice choices.

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

## Reproduce the fine-mapping audit and comparison preparation

```bash
python scripts/fetch_sources.py --manifest config/finemap-sources.json
python scripts/audit_finemap_archive.py
python scripts/prepare_controlled_comparison.py
python -m unittest discover -s tests -v
```

The 43 pinned sources include the 16.7 MB study archive and small reference/documentation files. The new scripts use the standard library. `--offline` on the download command verifies the cache. The complete compact [variant table](data/derived/psc-finemap-variants.csv.gz) preserves original identifiers and probabilities; it is not a normalized VCF or a per-signal credible-set manifest. The [candidate table](data/derived/psc-comparison-candidates.csv) includes every high-PIP inclusion, exclusion and reserve in the selected regions.

The [protocol](config/psc-controlled-comparison.json) fixes matching, model and reporting rules before scores. The [original freeze record](config/psc-controlled-comparison-lock.json) preserves that preparation phase and its then-incomplete status. The later [final input manifest](config/psc-matched-comparison-inputs.json) locks the eligible matched groups and common gene/track universes before the nine new predictions. The original pilot predictions are unchanged.

## Reproduce the fixed matched comparison

Use Python 3.12 with both `requirements-alphagenome.lock` and `requirements-qtl.txt` installed. After the fine-mapping preparation above:

```bash
python scripts/fetch_comparison_sources.py
python scripts/prepare_matching_covariates.py
python scripts/fetch_comparison_genotypes.py
python scripts/match_comparison_variants.py
python scripts/fetch_comparison_metadata.py
python scripts/freeze_comparison_inputs.py
python scripts/run_matched_comparison.py
python scripts/summarize_matched_comparison.py data/predictions/20260912T043400Z-matched-comparison
python -m unittest discover -s tests -v
```

The 85 matching-source pins include a 333 MB model gene annotation; genotype queries are bounded indexed regional extractions, with donor-level subsets kept ignored. Each fetch command accepts `--offline` to verify cached inputs. The metadata helper requests only track metadata if its cache is missing; it needs the locally configured API key in that case and rejects any change from the pinned snapshot. `freeze_comparison_inputs.py` verifies the existing freeze without replacing its timestamp. The runner defaults to a dry run and requests no predictions without `--run`.

The summary command requires the saved local September 12 prediction folder. It is not included in Git. A new API run with `python scripts/run_matched_comparison.py --run` creates a new dated folder; pass that printed path to the summary command. Preserve the historical reports before summarizing a later run, since the summary filenames are fixed. A later server model or live metadata change may prevent exact historical replay; the recorded inputs, compact outputs and hashes remain the original result. Public APIs may also return changed or differently serialized responses; the downloader stops on a pinned-hash mismatch rather than silently updating the inputs.

## Reproduce the CD4 RNA evidence

Use the pinned Python 3.12 environment above. The repository includes the frozen plan and all 90 small aggregate association subsets; it contains no individual-level CD4 data. Retrieve the 36 pinned supporting sources and reproduce the analysis:

```bash
python scripts/fetch_cd4_rna_sources.py
python scripts/prepare_cd4_rna.py
python scripts/fetch_cd4_rna.py --offline
python scripts/summarize_cd4_rna.py
python -m unittest discover -s tests -v
```

`prepare_cd4_rna.py` verifies the existing plan, timestamp and complete 228-pair registry. `fetch_cd4_rna_sources.py --offline` checks the source cache. `fetch_cd4_rna.py --refresh` repeats all 90 rate-limited indexed queries and requires exact agreement with the saved responses; cached results are never silently replaced. No AlphaGenome key or prediction call is required for this phase.

The [full report](reports/cd4-rna-evidence.md) includes all ten datasets, every missing-measurement category, the two inconsistent Treg metadata labels, source allele/carrier counts, and all screened associations. RNA coefficients are normalized-expression effects per ALT allele; their signs do not determine PSC risk or treatment direction. Prior model results and their fixed matching protocol remain unchanged.

## Interpretation

The [liver-atlas reproduction guide](docs/liver-atlas-reproduction.md) covers the pinned Python environment, the 2.66 GB public count download, metadata audit, donor-level analysis, figure and isolated numerical replay. `python scripts/verify_liver_atlas.py` checks the compact versioned results without requiring the large source cache. The frozen analysis plan and full output tables retain sparse groups, zeros and unavailable disease comparisons.

The [colocalisation reproduction guide](docs/coloc-reproduction.md) provides portable hash checks, figure/table regeneration, isolated numerical replay of both baselines, and the additional requirements and current cache-restoration limit for full LD refitting. The [complete report](reports/colocalisation.md) and [next input specification](docs/coloc-next-inputs.md) explain the findings and the bounded follow-up.

The completed [variant-coverage follow-up](reports/variant-coverage-repair.md) adds exact sequence normalization and bounded reannotated OneK1K queries. Its [reproduction guide](docs/coverage-repair-reproduction.md) and [Python dependency pins](requirements-coverage-repair.txt) explain cache verification, source limits and the separate coverage-only analysis. Run `python scripts/verify_coverage_repair.py` to check versioned results without requiring ignored source caches.

PSC susceptibility is different from progression of established PSC. A regulatory prediction can support a mechanism worth testing; it does not show that altering that gene will safely treat disease. Experimental and clinical validation would still be needed.

This repository contains public-source research material and unsent draft questions. Personal health records and actual private correspondence are kept outside it. The watchlist is a dated snapshot and does not refresh automatically.

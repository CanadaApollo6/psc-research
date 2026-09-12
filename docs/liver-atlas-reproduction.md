# Reproduce the PSC liver-cell analysis

Recorded analysis: 2026-09-12. See the [findings](../reports/liver-atlas-cell-context.md), [frozen plan](../config/liver-atlas-plan.json), [source manifest](../config/liver-atlas-sources.json) and [expression audit](../reports/liver-atlas-expression-audit.json). This phase uses public RNA counts and metadata. It needs no API key or GPU.

## Environment and files

The recorded environment uses Python 3.12.14 and the six exact analysis-package versions in [requirements-liver-atlas.txt](../requirements-liver-atlas.txt). Use an isolated environment and run commands from the repository root:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-liver-atlas.txt
```

The complete repository test suite also covers earlier phases; its additional requirements are documented in the README and preceding reproduction guides. The atlas-specific tests use the environment above.

The two H5AD files total 2,662,773,805 bytes, approximately 2.66 GB. Allow additional space for the environment and intermediate files. The reader processes sparse counts in chunks rather than loading a dense expression matrix. Downloaded matrices stay in ignored `data/raw/liver-atlas/`; the pseudonymous donor-level expression cache stays in ignored `work/`. Only aggregate expression tables, public source metadata, scripts, manifests and reports are versioned.

## Verify the small versioned results

After installing the atlas dependencies, this command checks the frozen plan, source manifest, script and requirements hashes; the four aggregate-output hashes; cell accounting; gene coverage; donor gates; and contrast statuses:

```bash
python scripts/verify_liver_atlas.py
```

A missing donor cache is allowed for this check and is reported explicitly. If the cache exists, its hash and logical consistency are checked. This verifies the recorded tables; it does not reproduce their derivation from RNA counts.

## Retrieve and reproduce from counts

```bash
python scripts/fetch_liver_atlas_sources.py
python scripts/prepare_liver_atlas_metadata.py
python scripts/analyze_liver_atlas.py --chunk-size 8192
python scripts/verify_liver_atlas.py --require-donor-cache
python scripts/plot_liver_atlas.py
python -m unittest discover -s tests -p 'test_liver_atlas*.py' -v
```

The downloader streams to temporary files, checks exact byte counts and SHA-256, and preserves a valid cache. Add `--offline` to verify an existing cache without a download. The two dataset URLs are version-specific; the collection metadata URL is a live endpoint whose serialization or content may change. A changed response stops reproduction instead of silently replacing the pinned snapshot. Restore the original cached response or document a separate source revision if that happens.

The metadata builder verifies the pinned collection, both H5ADs and the versioned GEO CSV before parsing them. It reads H5AD observation metadata, not expression values. Paper, preprint-supplement and author-code facts are separately curated with provenance in [liver-atlas-metadata.json](../config/liver-atlas-metadata.json); the metadata builder does not re-extract those facts from the original documents. The 29-row donor/object table is not an exact crosswalk to the 45 GEO libraries.

The expression analyzer verifies the original plan and dataset hashes, including when `--dataset KEY=PATH` overrides a local filename. It then checks all stored raw counts and sparse-matrix structure, preserves author annotations, applies the frozen exclusions, and aggregates at donor × object × assay × cell population. All-feature counts supply each donor's CPM denominator; each donor has equal weight in the reported medians. The two modalities and PBC remain separate. Primary and sensitivity thresholds are 20 and 50 cells per donor/population; descriptive contrasts require three donors per cohort in the same assay. Whole-cell cross-chemistry comparisons are blocked. Nuclear contrasts retain the residual processing-confounding warning.

The original plan was saved before candidate expression inspection at `2026-09-12T14:15:24.857808+00:00`. Its SHA-256 is `961ad6eec38c25b352e32e5af51f9af9a663da138286ba8815f44ea2049b4a6b`. Do not edit the plan to accommodate a different input or preferred result. A changed biological question requires a separately recorded analysis.

## Isolated numerical replay

An alternate output directory preserves the recorded results while rerunning the same analysis with different chunk boundaries:

```bash
python scripts/analyze_liver_atlas.py \
  --output-directory work/liver-atlas-replay \
  --audit-path work/liver-atlas-replay-audit.json \
  --chunk-size 4096
```

Compare the five output SHA-256 values by basename with `outputs_sha256` in the recorded expression audit. All four public CSVs and the donor cache should agree exactly; output paths and count-processing chunk totals in the replay audit differ. `verify_liver_atlas.py --require-inputs` independently reruns full input count checks, but it is not a substitute for comparing the regenerated output hashes.

The metadata table can likewise be regenerated using `--output-csv` and `--output-report` with paths inside `work/`. The CSV should agree exactly. The metadata report records the actual generation time and requested output paths, so compare its substantive fields rather than expecting an identical report hash.

The recorded analysis also has an independent direct-matrix check for six selected gene/population pairs: ETS2 and PRKD2 in monocytes, and UBASH3A in CD4 T cells, separately for cells and nuclei. Donor eligibility, median CPM and median detection fractions agree with the pipeline. These checks assess calculation correctness, not biological independence.

## ETS2 benchmark input inventory

This separate preparation downloads about 35 MB of public article/release material and inventories available gene-set, rank-statistic and spatial files:

```bash
python scripts/fetch_liver_atlas_sources.py --manifest config/ets2-benchmark-sources.json
python scripts/inspect_ets2_benchmark.py
```

Use `--offline` on the fetch command to check cached sources. The inventory helper verifies every source hash, the publisher ZIP MD5 and selected archive-member hashes. It does not execute the archive's R code, extract nested spatial matrices, calculate a program score or retrieve controlled-access RNA data. See the [next-input specification](ets2-benchmark-next-inputs.md) for the next bounded analysis.

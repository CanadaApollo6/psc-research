# Reproduce the colocalisation assessment

The repository contains the compact public association subsets, complete target-gene published Bayes-factor vectors, frozen plans, code and result tables. It does not contain reference-donor genotypes, large LD matrices, R fit objects or the installed runtime. No AlphaGenome key is needed.

## Verify a fresh checkout

Use Python 3.12. The first command requires only the standard library. It verifies versioned artifacts and reports ignored cache files separately; an absent cache is not a reproduced fit.

```bash
python scripts/verify_coloc_artifacts.py
python -m pip install -r requirements-coloc-python.txt
python scripts/summarize_coloc.py
python -m unittest discover -s tests -v
```

The summary uses recorded input hashes, the saved association/posterior tables and complete published target-gene BF vectors. It rechecks published credible-set PIPs, recomputes missing-variant weights independently, restores coverage-failed signal rows, and regenerates the tables, JSON and PNG/SVG figures. Raw inputs that are present must match their pins; absent ignored raw inputs are not needed to summarize the versioned results. A copy with all 103 ignored cache files absent reproduced the ten generated report artifacts and both baselines exactly with the pinned installed runtime; see the [verification record](../reports/coloc-verification.json). Test discovery also covers earlier project phases, whose other dependencies are described in the README.

For the exact recorded environment, see [Python pins](../requirements-coloc-python.txt), the [Linux R package lock](../requirements-coloc-linux-64.lock), and [runtime provenance](../config/coloc-runtime.json). The recorded R runtime is R 4.4.3, coloc 5.2.3 and susieR 0.14.2. With micromamba available, the explicit lock installs the recorded package artifacts:

```bash
micromamba create --yes --prefix work/coloc-runtime/r --file requirements-coloc-linux-64.lock
OPENBLAS_NUM_THREADS=4 work/coloc-runtime/r/bin/Rscript --vanilla tests/test_coloc.R
OPENBLAS_NUM_THREADS=4 work/coloc-runtime/r/bin/Rscript --vanilla tests/test_coloc_platform.R
python scripts/replay_coloc_baseline.py
```

The replay computes **both complete ABF baselines in temporary workspace storage** and requires exact agreement with all six recorded output files. It does not overwrite historical reports and does not need reference genotypes or LD matrices. Pass `--rscript /absolute/path/to/Rscript` to use another installation with the same pinned R packages. A changed runtime or source serialization can prevent exact replay; do not update the historical hashes merely to accept a changed result.

## Verify the full original local cache

The original working directory retains all additional inputs. These commands verify them without new network access or fitting:

```bash
python scripts/verify_coloc_artifacts.py --require-cache
python scripts/fetch_coloc_sources.py --offline
python scripts/fetch_coloc_regions.py gwas --offline
python scripts/fetch_coloc_regions.py qtl --offline
python scripts/fetch_coloc_regions.py reference --offline
python scripts/fetch_coloc_signal_bfs.py --offline
python scripts/execute_coloc.py baseline --offline
python scripts/execute_coloc.py fit --dataset GWAS --gene PFKFB3 --offline
python scripts/execute_coloc.py fit --dataset GWAS --gene BCL2L11 --offline
python scripts/execute_coloc.py fit --dataset QTD000031 --gene PFKFB3 --offline
python scripts/execute_coloc.py fit --dataset QTD000031 --gene BCL2L11 --offline
python scripts/execute_coloc.py compare --offline
python scripts/execute_coloc_platform.py baseline --offline
python scripts/execute_coloc_platform.py fit --gene PFKFB3 --offline
python scripts/execute_coloc_platform.py fit --gene BCL2L11 --offline
python scripts/execute_coloc_platform.py compare --offline
```

Replacing `--offline` with `--rerun` on an executor repeats that stage and requires byte-identical numerical outputs. It overwrites the stage's local output files while checking them; use the isolated baseline helper above for a replay that preserves historical files. The published-signal comparison needs the fitted GWAS R objects in the ignored cache. Its repeated result, as well as both isolated baselines, was verified against the recorded hashes.

## What full LD refitting needs

Full refitting from a fresh checkout additionally requires the two ignored, sample-bearing EUR genotype subsets in [the reference query record](../config/coloc-reference-queries.json), plus the chain files and supporting sources in [the source manifest](../config/coloc-sources.json). The fixed reader's existing-manifest path expects these recorded subsets to be present; it does not automatically restore a missing cache. That is a current acquisition limitation, not a claim that a fresh checkout can refit everything with one command.

The source record specifies the public indexed VCF URLs, exact intervals, 503 EUR donor selection, parser version, original query limits, output filenames and hashes. A reacquisition must use those exact rules and match the compressed and uncompressed subset hashes before reuse. Preserve the historical manifest instead of deleting it or replacing its hashes. `prepare_coloc_inputs.py` and `prepare_coloc_platform_inputs.py` then regenerate allele alignment and signed LD matrices, which must match the recorded input locks. The original and amended R fitting scripts record all settings and output hashes. This phase was run in the original cached workspace; a complete fresh-download LD refit has not been independently replayed.

Complete PFKFB3 published BF subsets are already versioned, so reproducing their use does not require another 3.26 GB transfer. Their query record includes full compressed archive hashes and end-of-file verification. The [original signal-comparison CSV](../reports/coloc-published-signal-comparison.csv) keeps the software's output unchanged; use the derived [complete status table](../reports/coloc-signal-status.csv) for interpretation of coverage failures.

# Reproduce the variant-coverage follow-up

The [report](../reports/variant-coverage-repair.md) preserves all original colocalisation results and adds a separate exact-identity/coverage audit. The [main follow-up plan](../config/coverage-repair-plan.json) was frozen before new targeted source queries. The [OneK1K selection](../config/coverage-repair-alternative-plan.json) and two access amendments precede target-effect inspection. The initial TSV 404s and metadata-reader failure remain recorded; the archived first reader is not the active implementation.

The repository includes public aggregate query subsets, derived identities, all coverage/exclusion tables and source hashes. Reference genotype rows, downloaded DNA/documents, Parquet byte ranges and the runtime are ignored. No API key is needed.

Verify the versioned files in a fresh checkout with the standard library:

```bash
python scripts/verify_coverage_repair.py
```

This checks all available pinned files and reports absent ignored inputs. It does not claim that missing source caches have been reproduced. The recorded runtime is Python 3.12.14 with the [dependency pins](../requirements-coverage-repair.txt), including PyArrow 25.0.1. Test discovery covers the earlier project phases as well; retain their README dependencies.

In the original workspace, all source caches are available. These commands need no network and require existing bytes to match their hashes:

```bash
python scripts/verify_coverage_repair.py --require-cache
python scripts/fetch_coverage_repair_sources.py --offline
python scripts/fetch_coverage_repair_regions.py reference --offline
python scripts/fetch_coverage_repair_regions.py catalogue --offline
python scripts/fetch_coverage_repair_parquet.py --offline
python scripts/analyze_coverage_repair.py
python scripts/analyze_coverage_repair_alternative.py
python -m unittest discover -s tests -v
python scripts/verify_coverage_repair.py --require-cache
```

The exact-identity analysis validates full reference sequences, reciprocal mappings, alternate haplotypes and explicit phase-3 REF/ALT records. A second implementation, bundled `bcftools norm`, checks normalization for both builds. The complete published PFKFB3 BF vectors already in Git provide the signal denominator; no large molecular BF archive or AlphaGenome rerun is required.

The OneK1K reader identifies every Parquet row group overlapping the fixed position window, reads chromosome/gene identifiers to select exact BCL2L11 groups, and only then reads the full association columns. The r8 footer omits min/max for string columns, so it cannot be used to prune directly by gene. The [second access amendment](../config/coverage-repair-alternative-parquet-query-plan.json) records this change. Each transferred range has a SHA-256, exact byte interval, remote size and ETag/Last-Modified identity; full nominal data remain capped at 300 MB and one million candidate rows per dataset. Identifier-only row counts are recorded separately. All source gene/region rows are retained before SNP eligibility or alias handling.

For a fresh source reacquisition, `fetch_coverage_repair_sources.py` can retrieve the bounded reference DNA, source documents and small credible-set files using the pinned request specifications. The two reference-genotype subsets additionally need the original index/panel/chain sources in `config/coloc-sources.json`. The existing indexed-query and Parquet readers verify a recorded cache; their recorded-cache paths do not automatically restore deleted subsets/ranges. Reacquisition must preserve the query manifests and match their compressed/range hashes, rather than delete manifests or replace pins to accept new bytes. This is an explicit cache-restoration limitation. Portable hash verification does not depend on those ignored inputs.

The original DICE “unfiltered” TFH/TH17/TH2 files were inspected only through 64 KiB header probes. Their 10.17 GB combined contents were not extracted, and source-level target absence was not asserted. The versioned [access record](../config/coverage-repair-original-dice-access.json) preserves that boundary.

All new results use the frozen effect/allele rules and fixed coverage gates. PFKFB3's added support is not a refitted GWAS component. OneK1K's coverage uses fixed-SD quantitative Bayes factors without an exact cohort-N assumption; no new H0–H4 posterior or estimated-SD analysis was run. Published credible-set absence and missing RNA records remain distinct from measured null effects.

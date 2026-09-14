# Reproducing the UBASH3A public-data hypothesis test

The [result](../reports/ubash3a-hypothesis-test.md) is an annotation-level structural finding and an unestimable NMD contrast. It is not evidence for RNA survival or stable protein. The original [method](ubash3a-hypothesis-test-protocol.md), [annotation input specification](../config/ubash3a-test-annotation-inputs.json) and [original-read amendment](../config/ubash3a-upf1-read-test-plan.json) define the endpoints.

## Runtime and data

Run from the repository root. The existing `.venv` already contained all required packages; this follow-up installed none. Its direct dependencies are pinned in [requirements-ubash3a-hypothesis-test.txt](../requirements-ubash3a-hypothesis-test.txt): the previous requests/Biopython/NumPy/Matplotlib lock plus pysam 0.23.3. The complete repository suite also uses dependencies of earlier analyses. This is not a claim that a fresh-host environment or all earlier studies have been reinstalled and replayed here.

Raw files are ignored by Git. A fresh clone therefore needs the public inputs listed in [ubash3a-hypothesis-test-sources.json](../config/ubash3a-hypothesis-test-sources.json), the eight BAM subsets in [ubash3a-upf1-read-sources.json](../config/ubash3a-upf1-read-sources.json), and the already pinned earlier reference caches. The first manifest contains 70 usable objects and one recorded failed GEO family request. Small output tables and the plotted figure are versioned.

The general fetcher uses `requests` for explicit public files, the installed NCBI Entrez research-skill client for Entrez, and the BioStudies research-skill client for study metadata. The clients' hashes and saved representations are recorded. Existing source pins are never silently replaced; a re-download must exactly match its historical SHA-256 and size. An incomplete `.part` file requires inspection before another attempt. Do not run two general-fetcher processes simultaneously because they share one manifest.

For reacquisition, make a temporary request list from the manifest's successful `sources` array, excluding the known failed request, and run the fetcher once. The installed research-skill paths on the original host were:

```bash
.venv/bin/python - <<'PY'
import json
from pathlib import Path
source = json.loads(Path('config/ubash3a-hypothesis-test-sources.json').read_text())
target = Path('work/ubash3a-hypothesis-test/reacquire.json')
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text(json.dumps({'requests': source['sources']}, indent=2) + '\n')
PY
.venv/bin/python scripts/fetch_ubash3a_test_sources.py \
  work/ubash3a-hypothesis-test/reacquire.json \
  --entrez-client /home/riels/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/ncbi-entrez-skill/scripts/ncbi_entrez.py \
  --biostudies-client /home/riels/.codex/plugins/cache/openai-curated-remote/life-science-research/1.0.3/skills/biostudies-arrayexpress-skill/scripts/rest_request.py
```

On another host, supply the corresponding client paths. Metadata that have changed since the frozen capture will correctly fail a hash check; that is not permission to update the analytical snapshot. The failed GSE329233 family SOFT request remains documented alongside its zero-result Entrez lookup.

## Original-read extraction

```bash
.venv/bin/python scripts/fetch_ubash3a_upf1_reads.py
# With complete local caches, verify without any network access:
.venv/bin/python scripts/fetch_ubash3a_upf1_reads.py --offline
```

The general fetcher retrieves the original BAI files. With those indexes present, the BAM extractor saves each header/prefix and uses indexed access to the frozen UBASH3A interval. It checks a 64-KiB HTTP range response, the GRCh38 chromosome-21 length, remote size/ETag/modification time, subset readback, and output hashes. No whole source BAM was downloaded or hashed. All overlapping alignments, including the lone unspliced reverse alignment, are retained before analytical quality checks.

Missing historical subset caches can be re-extracted into separate ignored staging directories. Every resulting file must match its old hash and size, and the remote identity must match the original extraction before missing files are restored. Existing changed files are rejected. Restoration receipts stay in `work/`; historical retrieval metadata remain intact. This restoration helper was added after the initial extraction, whose original fetch-code hash is preserved in the source manifest. A successful fresh-host restoration has not been claimed.

## Deterministic analyses

```bash
.venv/bin/python scripts/analyze_ubash3a_observed_isoforms.py
.venv/bin/python scripts/analyze_ubash3a_upf1_reads.py
.venv/bin/python scripts/summarize_ubash3a_test_support.py
.venv/bin/python scripts/verify_ubash3a_observed_isoforms.py
.venv/bin/python scripts/verify_ubash3a_test_support.py
MPLCONFIGDIR=work/ubash3a-hypothesis-test/matplotlib \
  .venv/bin/python scripts/plot_ubash3a_hypothesis_test.py
.venv/bin/python -m unittest discover -s tests -v
```

The first three analyses accept `--output-dir`. The read and support analyses intentionally consume the pinned primary annotation/reference outputs even when their own output directory is changed. Separate replay directories were compared byte for byte against every primary file in their audit manifests, including the audit JSON itself. The figure is exported as PNG and SVG; the latter normalizes only trailing line whitespace for a clean Git diff.

The original annotation run stopped when a GTF attribute contained an unquoted number. The [repair record](../config/ubash3a-hypothesis-test-execution-repairs.json) explains the parser extension and preserves before/after code hashes. The plan, coordinate endpoint, source set and biological classification rules were unchanged. The original annotation-run and read-analysis-run records remain historical snapshots rather than being overwritten with completion status.

## Interpretation boundaries and checks

- GTF/GFF coordinates are one-based and inclusive. The independent source checker constructs zero-based intervals and independently parses all 13 files, including anonymous model IDs and distant exons belonging to selected transcripts.
- All matching source rows, original labels and support strings remain available. Candidate-model counts are distinct from read counts and donors. The TRAILS matrix lacks one negative-strand locus model; that absence is not imputed to zero. Its precise matrix export/normalization step is unresolved.
- The UPF1 read rule uses CIGAR `N` for splicing, with 20 aligned bases on each evaluated side, and retains complete-chain versus downstream-only assignments separately. Deletions are not introns; low-quality, secondary and supplementary alignments are not qualifying evidence. No full target RNA is covered in these extracts, so no effect estimate or P value exists.
- All 12 diagnostic RNA windows are reference constructions. The short versions omit rs1893592; only the longer templates distinguish hypothetical A/C alleles. Peptide identity does not uniquely assign its originating RNA.
- Source URLs, retrieval times, file sizes, hashes, independent checks, replay results, the full-suite result and preservation of earlier tracked files appear in the [completion verification](../reports/ubash3a-hypothesis-test-verification.json).

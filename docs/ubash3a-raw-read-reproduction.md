# Reproducing the initial TRAILS read audit

The fixed scope is the complete NCD4 and MemCD4 libraries in [the plan](../config/ubash3a-raw-read-plan.json). Their source totals are 8,401,078 reads and 8,050,358,628 bases. These are amplified cDNA reads from one donor, not original-molecule or independent-donor counts.

## Inputs and runtime

Restore the small pinned source caches named in the raw-read metadata manifest and both analysis plans. They include the five original DRA XML files, tool-release/method metadata, the prior UBASH3A gene-reference response and UniProt record, and the two ENA single-run totals. Earlier derived exon/protein tables are versioned. The fresh Git repository alone cannot reproduce every check offline because original source caches and large sequence/alignment files are ignored.

The three new large inputs are in [the download manifest](../config/ubash3a-raw-read-downloads.json): two DDBJ compressed FASTQ libraries and the Ensembl release-112 GRCh38 primary assembly. Compressed transfer is 7.12 GB. Allow roughly 30 GB local space for reference, index, retained downloads and target alignments; the declared preparation budget was 50 GB. Minimap2 indexing peaked at 16.475 GB RAM. Mapping uses eight workers, 50-million-base batches and one library at a time. Do not launch two genome mappings concurrently on this project's 32-GB host.

The [runtime record](../config/ubash3a-raw-read-runtime.json) pins existing Python dependencies and Minimap2 2.31-r1302. Only Minimap2 was newly installed, after explicit approval, in ignored project storage. The [installation record](../config/ubash3a-raw-read-install-plan.json) includes the official release URL and publisher SHA-256. No system Python/package changes were made.

## Commands

From the research repository, with its existing environment and restored small caches:

```bash
.venv/bin/python scripts/fetch_ubash3a_raw_reads.py --reference-ranges
.venv/bin/python scripts/run_ubash3a_raw_read_alignment.py
.venv/bin/python scripts/analyze_ubash3a_alternative_ending.py
.venv/bin/python scripts/analyze_ubash3a_raw_reads.py
.venv/bin/python scripts/verify_ubash3a_raw_read_audit.py
.venv/bin/python scripts/summarize_ubash3a_raw_read_sensitivity.py
.venv/bin/python -m unittest discover -s tests -v
```

The reference-range option speeds transfer from the same pinned Ensembl file. Every range must have the exact requested boundaries and ETag. The assembled genome also passed the publisher's BSD checksum and gzip integrity, and its UBASH3A sequence agrees exactly with the preceding gene reference. This primary assembly contains 194 sequences; alternate haplotypes and patches are not part of the alignment competition.

The aligner validates every FASTQ record and complete bzip2 decompression, then checks equal input-read/output-query-group totals and one primary or unmapped record per output group. It streams whole-genome alignments while retaining all alignments for each query with any locus overlap, including competing secondary/supplementary mappings. Target queries are also reconstructed in their original orientation for the separate splice-flank sensitivity run. This retains source read sequences and qualities via the primary SAM record; header comments are not preserved in that selected FASTQ.

Each completed library has an `alignment.json` record in ignored work storage. It is written only after successful primary mapping, complete input checks, target BAM/index creation and the second target-read mapping. Classifier output requires both libraries' completion records. These records and raw caches must be kept to reproduce the read audit from existing alignments.

## Interpretation and verification boundaries

Primary eligibility and the 10-base-anchor sensitivity are fixed in [the protocol](ubash3a-raw-read-protocol.md). A matching local junction, a downstream chain, and a full splice chain plus aligned shared start are different categories. None alone certifies a transcript's physical ends or an error-free full coding sequence. Minimap2's `ts` is query-relative; reverse SAM alignments require sign conversion to infer genomic transcript strand, as verified with synthetic forward/reverse reads.

The [additive coding descriptors](../config/ubash3a-raw-read-coding-descriptors.json) were defined after the independent T29 translation calculation, while alignments were running but before target-junction or architecture results were inspected. They report coverage of its newly calculated stop coordinates and the catalog end. They do not alter eligibility, primary counts or the original fixed plan, and do not prove translation or RNA termination.

The separate verifier reconstructs junctions/anchors from textual CIGAR, checks all read eligibility and full/downstream matches, compares complete-library totals to ENA, and independently recomputes affine-gap scores for every available primary +29 candidate plus ten normal candidates per library. Source metadata joins and file hashes are checked too.

The sensitivity-summary command additionally counts exact normal, +29 and T29-specific downstream junctions across all nominated-read output alignments under both settings, including secondary/supplementary records. It preserves neighboring junction coordinates separately. Its source hashes and two table hashes are recorded in `raw-sensitivity-summary-audit.json`; these descriptive tables do not change the frozen eligibility rules.

For an exact classifier replay, run the classifier with `--output-dir` pointing to an ignored temporary folder and compare its four tables and audit against the versioned outputs. Alignment completion timestamps and compressed-file container bytes may change in a fresh alignment run; compare scientific counts and read identities as well as the recorded source/method provenance. No alternative aligner setting replaces a failed primary result.

The [stage verification](../reports/ubash3a-computational-followup-verification.json) preserves the complete reference/index provenance and checksums, source-cache verification, all-setting event/neighbor checks through the separate textual-CIGAR parser, replay hashes and prior-file preservation against the pre-stage commit. Its software-test result refers to this execution, not to an experiment.

The separate [protein-data inventory](../reports/ubash3a-proteomics-readiness.md) is metadata preparation. Rebuild it with `.venv/bin/python scripts/summarize_ubash3a_proteomics_readiness.py`; this does not download or search mass spectra.

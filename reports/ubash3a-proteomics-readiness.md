# Public CD4 proteomics: file inventory prepared

The next protein-evidence step has a concrete public input: [PRIDE PXD000376](https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD000376), *A Multi-Omic Analysis of Human Naive CD4 T cells*. The complete inventory contains **273 files totaling 92.44 GB** (decimal units). This is metadata preparation; no spectrum file was downloaded or searched.

| Source category | Files | Listed bytes | Formats |
|---|---:|---:|---|
| RAW | 249 | 91,304,122,023 | `.raw` |
| SEARCH | 24 | 1,135,736,484 | Peptide/protein search XML |
| Total | 273 | 92,439,858,507 | |

The [inventory](../data/derived/ubash3a-proteomics-readiness/PXD000376-file-inventory.csv) preserves file names, exact source categories, sizes and public locations. Three metadata pages contain 100, 100 and 73 unique file records, agreeing with the separate project file-count endpoint. All records identify the selected project. None supplies a nonempty checksum in the returned file-checksum field; this does not imply that checksums are unavailable through every other archive route.

## Context and limitations

The project describes primary naïve and resting memory CD4 cells, fractionated tryptic proteomics, and iTRAQ comparisons. Its instruments are LTQ Orbitrap Elite and Velos. The metadata specify an older human RefSeq database searched using Mascot/Sequest, followed by a 1% decoy-based peptide filter. These describe the original analysis; they do not establish error control for newly proposed T29 peptides.

The submission is labeled `PARTIAL`. The [paper](https://link.springer.com/article/10.1186/s12918-015-0225-4) describes the core multi-omic comparison as a single-donor study; the file count is not a donor count. File names alone have not been accepted as cell-condition, replicate or reporter-channel assignments. Both the archive and the paper contain the same “CD45RO depletion” wording in the memory-cell isolation description, so reading the paper did not resolve that wording. The paper also reports additional proteogenomic analyses; the precise database behind each search XML remains to be checked.

## Next execution gates

1. Check the [four candidate peptides](ubash3a-alternative-ending.md) against a versioned complete human protein/isoform search space, including I/L equivalence, other candidate products and contaminants. Absence from the twenty preceding modeled products is only the initial comparison.
2. Resolve file-to-sample/fraction relationships, search parameters and database contents. Assess access to readable spectra and conversion tools before selecting a bounded file subset or transferring the full deposit.
3. Freeze the candidate and background sequence database, decoys, enzyme/modification/mass-tolerance rules and error controls for novel assignments before inspecting target hits. Existing global FDR does not by itself settle a new candidate's credibility.
4. Require interpretable, discriminating fragment evidence and assess competing assignments. Report peptide evidence, candidate-product attribution and source-RNA attribution separately. A missing match in an earlier database that omitted the candidate is not evidence of biological absence.

The [source ledger](../config/ubash3a-proteomics-readiness-sources.json) records the mandatory PRIDE skill client, requests, snapshot format and hashes. These JSON snapshots were parsed and reserialized by the client, not captured as original HTTP body bytes. The live API documentation supplied the current project-file endpoint after an older endpoint did not yield a parseable JSON response; that failure remains recorded and is not treated as a zero-file result.

The [readiness audit](../data/derived/ubash3a-proteomics-readiness/PXD000376-readiness-audit.json) explicitly records zero searched spectra and untested human-proteome specificity. Direction 3 has begun with source preparation; experimental protein support remains unresolved. PRKD2's shared-signal analysis remains the fourth direction in the [roadmap](../docs/computational-followup-roadmap.md).

Rebuild the two preparation files with `.venv/bin/python scripts/summarize_ubash3a_proteomics_readiness.py`. Both replay byte for byte against the pinned snapshots.

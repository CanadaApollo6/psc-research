# Reproducing the UBASH3A consequence audit

The completed [report](../reports/ubash3a-transcript-consequences.md) answers a fixed sequence question. This workflow does not estimate splice frequency or disease effects. The [plan](../config/ubash3a-consequence-plan.json), [protocol](ubash3a-consequence-protocol.md), [source manifest](../config/ubash3a-consequence-sources.json) and [initial execution record](../config/ubash3a-consequence-run.json) preserve the pre-calculation choices. Completion is recorded separately in the [verification record](../reports/ubash3a-consequences-verification.json); the initial execution record remains unchanged.

## Inputs and runtime

Run commands from the repository root. The successful environment used Python 3.12.14, Biopython 1.88 and Matplotlib 3.11.2; the exact package pins are in [requirements-ubash3a-consequences.txt](../requirements-ubash3a-consequences.txt). The existing research environment was reused, adding the previously absent Biopython and requests packages and their dependencies. No original installed package version was replaced. The primary reconstruction itself uses Python's standard library.

A fresh Git clone includes the derived FASTA files, coordinate tables, reports and manifests, **but not the raw third-party response caches**. Full replay requires all 26 cache files listed in the source manifest plus the nine earlier files pinned in the plan, including the earlier expanded gene annotation. Their paths and hashes are explicit. The earlier literature cache is needed for the recorded Ge 2018 reuse. These are reproducibility inputs, not optional files to skip when unavailable.

Ensembl and UniProt calls use the installed research skills' `scripts/rest_request.py` clients. The fetch script locates them automatically, or accepts explicit `--ensembl-client` and `--uniprot-client` paths. New downloads require network access; offline verification does not. Public providers can change releases or JSON representations, so downloading a different response does not reproduce the frozen snapshot. Hash mismatches fail without silently changing historical pins; an updated-source analysis would need a separate amendment and outputs.

## Offline verification and replay

```bash
.venv/bin/python scripts/fetch_ubash3a_consequence_sources.py --offline
.venv/bin/python scripts/analyze_ubash3a_consequences.py --output-dir work/ubash3a-consequences/replay
.venv/bin/python scripts/verify_ubash3a_consequences.py
.venv/bin/python -m unittest discover -s tests -v
MPLCONFIGDIR=work/ubash3a-consequences/matplotlib .venv/bin/python scripts/plot_ubash3a_consequences.py
```

The primary reconstruction can write to a separate directory, preserving the published outputs. Compare all ten replay files to `data/derived/ubash3a-consequences/` byte for byte. The independent checker reads that published directory, reconstructs RNA by extending genomic exon intervals, translates with Biopython, and verifies the scenario tables, FASTAs, exon coordinates and protein-feature projections. Its timestamped verification record is rewritten when rerun. Plot regeneration rewrites the two figure files; figure appearance was visually checked in the original run.

The 21 new synthetic tests exercise independent known sequences and important failure cases: exon/strand boundaries, CDS stop conventions, first-stop translation, coding retention, 50/55-nt cutoffs, start proximity, exon length and ambiguous feature mapping. The full repository suite passed **165 tests**. The independent checker passed **69,547 scalar comparisons** over the real pinned sources and outputs. These are computational checks sharing source annotation; they do not prove a synthetic full-length RNA or protein exists in cells.

## Coordinate and feature conventions

- Genomic positions and reported cDNA/CDS/amino-acid positions are **one-based** unless a column explicitly says `0based`. Exon cDNA ranges use zero-based, half-open intervals. A cumulative exon-end boundary is therefore numerically equal to that exon's final one-based nucleotide.
- FASTA mature RNAs are represented in the **DNA alphabet**, using T. ORFs include the first in-frame stop; protein sequences omit the terminal stop symbol. Translation starts at the annotated start, not at a newly selected ATG.
- The primary NMD flag is `final_junction_boundary - stop_end_cdna_1based > 55`. The sensitivity flag substitutes 50. Both were fixed before calculation. The independent record also reports distances from the stop's first base, exposing the 49/51-nt boundary issue in ENST00000398367.2.
- `reference_residues_present` in the feature table counts **identical amino acids encoded by the same genomic codon triplet**. It is not a general sequence-similarity score. A changed codon can yield the same amino acid; canonical residue 465 does so here, explaining 70 exact-codon feature matches versus 71 unchanged amino acids within the phosphatase-like region.
- `novel_peptide_tail` means the translated sequence after the first difference from that transcript's reference peptide. It is a computed replacement sequence, not a claim that the sequence or splice event was previously unknown.

All first-run input/method hashes still match their pins. No analysis parameter, source sequence or original output was changed after seeing the results. The independent checker, visualizations and reporting were added after calculation as planned.

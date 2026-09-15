# UBASH3A candidate-ending peptides in public CD4 proteomics

**Completed within the fixed fourteen-file experiment: none of the four candidate-ending peptides appears among the exported assignments. Ordinary UBASH3A also lacks a qualifying assignment, limiting the biological meaning of the candidate nondetection.**

This is direction 3 of the [computational roadmap](../docs/computational-followup-roadmap.md). It asks whether public CD4 mass spectra support the distinctive ending of the conditional TRAILS T29 protein. The preceding original-read audit did not confirm T29's exact junctions in the two selected RNA libraries. A conditional translation and a stop in the final exon do not establish RNA survival or protein existence.

## Completed spectrum result

All fourteen searches finished successfully with the fixed database and parameters. **Zero assignments** match any of the four candidates, literally or with I/L equivalence, among **269,986 exported assignments**. All exported positions were examined, irrespective of rank or q-value; the engine exports up to five assignments per query. [Candidate results](../data/derived/ubash3a-proteomics-spectra/candidate-search-summary.csv), [complete summary](../data/derived/ubash3a-proteomics-spectra/search-summary.json).

| Fixed peptide | Exported assignments at any rank | Qualifying spectrum wins |
|---|---:|---:|
| RPNSSWKLGIFK | 0 | 0 |
| LGIFKIHLESGIR | 0 | 0 |
| IHLESGIR | 0 | 0 |
| IHLESGIRKP | 0 | 0 |

There were no candidate spectra to compare at fragment level. The [fragment-review record](ubash3a-proteomics-fragment-review.json) therefore contains zero candidate/competitor reviews, rather than an invented score or a peptide-detection claim. A missing exported assignment is not an exhaustive targeted comparison against every spectrum.

| Measurement or analysis level | Count |
|---|---:|
| Complete raw files | 14 |
| Converted original MS2 scans | 66,766 |
| Engine-loaded queries, each mapped to a distinct original scan here | 55,741 |
| Original scans with a reported peptide winner | 54,278 |
| Exported assignments across all positions | 269,986 |
| Target spectrum wins passing global spectrum q ≤ 0.01 | 11,611 |
| Distinct I/L-equivalent target peptides passing peptide q ≤ 0.01 | 2,978 |
| Ordinary-reference-compatible UBASH3A spectrum wins passing q ≤ 0.01 | 0 |

The spectrum and peptide thresholds are separate summaries, not a protein-level error estimate. All **55,741 loaded queries** appear in the XML output and map to the original instrument IDs. Of the 66,766 converted scans, 11,025 have no loaded query and 1,463 loaded queries have no reported assignment. All 12,488 scans without an assignment are accounted for by ID in the local audit. The engine does not supply a complete per-scan explanation for every preprocessing or no-competitor outcome. [Scan accounting](ubash3a-proteomics-scan-coverage.json), [per-file results](../data/derived/ubash3a-proteomics-spectra/file-search-summary.csv).

The database includes seventeen human UBASH3A reference records representing eight distinct sequences, including canonical P57075. The single [UBASH3A-compatible spectrum winner](../data/derived/ubash3a-proteomics-spectra/ubash3a-spectrum-winners.csv) has spectrum q ≈ **0.567** and peptide q ≈ **0.790**, far above the fixed threshold. It is not a positive control for UBASH3A detection. The absence of qualifying ordinary UBASH3A evidence makes this experiment a weak constraint on the proposed ending's biological existence, despite recovery of thousands of background peptide assignments.

**Decision:** the bounded public-protein step is complete and supplies no qualifying T29 protein evidence. Retain T29 as an unconfirmed conditional sequence; do not infer protein absence, RNA decay or a PSC mechanism. A stronger protein test would need informative UBASH3A coverage, suitable samples and calibrated novel-peptide validation. PRKD2 shared-signal analysis remains the next queued computational direction.

## Sequence-specificity result

All four previously fixed peptides are absent from the retrieved human-reference sequences, the common-contaminant collection and the twenty earlier modeled products. This holds for both literal matching and I/L-equivalent matching, which accounts for leucine/isoleucine mass indistinguishability. Complete protein sequences were scanned before applying cleavage rules; an inconvenient nontryptic competitor would not have been discarded.

| Fixed peptide | Length | Literal background matches | I/L-equivalent background matches |
|---|---:|---:|---:|
| RPNSSWKLGIFK | 12 | 0 | 0 |
| LGIFKIHLESGIR | 13 | 0 | 0 |
| IHLESGIR | 8 | 0 | 0 |
| IHLESGIRKP | 10 | 0 | 0 |

The [full match table](../data/derived/ubash3a-proteomics-specificity/all-sequence-matches.csv) contains sixteen rows: each peptide occurs in both conditional T29 allele products under both comparisons. All parent-sequence and frozen digestion checks pass. The short IHLESGIR is contained in two longer candidates, and all four come from the same short terminal region. These sequence definitions do not establish four independent protein observations or distinguish the hypothetical A/C products.

The [reference lock](../config/ubash3a-proteomics-reference-lock.json) was written at **23:10:03 UTC on September 14, 2026**, before new candidate matching. It covers:

| Source | Sequence records | Completeness check |
|---|---:|---|
| UniProtKB 2026_03, all retrieved human entries and available isoforms | 232,837 | All 466 pages; stable release/total headers; unique IDs; terminal pagination; 210,706 canonical entries plus 22,131 isoforms |
| GENCODE 50 comprehensive translations | 382,428 | Complete gzip stream and publisher MD5; GRCh38.p14, including annotated nonproductive and alternative translations |
| Original GPM cRAP collection | 116 | Original FTP file, published byte count and modification timestamp, local SHA-256 |
| Earlier conditional/reference UBASH3A models | 20 | Unchanged prior input hash |
| Conditional T29 allele products | 2 | Unchanged prior input hash |

These 615,403 source records yield **308,191 distinct complete protein sequences** after exact-sequence deduplication; all source memberships are retained. The candidate peptides have no matches to the inverse sequences that would generate them as Comet decoys. [Database record](../config/ubash3a-proteomics-search-database.json).

An independent implementation using Biopython, overlapping regular expressions and a separate cleavage-gap calculation agrees on all sixteen match rows and 272 fields. Independent XML parsing also agrees on all 250 deposited raw-file references. [Specificity verification](ubash3a-proteomics-specificity-verification.json).

This qualifies sequence targets within these references. It does not prove uniqueness across unannotated proteins, every population variant, other species or regions obscured by ambiguous residues. GENCODE translations and unreviewed records are possible sequence competitors, not necessarily observed proteins. Ambiguous residues remain in the source files; the search engine's handling is [recorded explicitly](../config/ubash3a-proteomics-sequence-handling.json).

## Experiment selection and source context

The deposited search headers map **all 249 archived raw files**, without cross-experiment file reuse. They contain 250 references because one separate iTRAQ labeling-check file is absent from the archive. That missing file is outside the selected experiment. All twelve search headers, original database names, tolerances and file mappings are retained in the [metadata table](../data/derived/ubash3a-proteomics-specificity/experiment-search-metadata.csv) and [file mapping](../data/derived/ubash3a-proteomics-specificity/experiment-raw-mapping.csv).

The prospectively selected experiment is **05BR2_nonNuclear_IOB**, with fourteen files and **2,387,832,397 listed bytes**. Selection used its explicit naïve-CD4/non-nuclear annotation, coherent file mapping and bounded transfer size, before new candidate matches or spectrum scores. It includes the source's pooled fractions and repeat injection. These do not establish independent donor replication. The publication's core study uses one donor; per-file donor keys are not supplied. No donor key links these spectra to the TRAILS RNA donor, so this is not a matched RNA–protein comparison. This search does not estimate a naïve–memory, genotype or disease contrast. [Frozen selection](../config/ubash3a-proteomics-spectrum-selection.json).

The [study's supplemental methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12918-015-0225-4/MediaObjects/12918_2015_225_MOESM25_ESM.doc) and deposited search settings support tryptic digestion and high-resolution HCD measurements. The unresolved CD45RO isolation wording is retained; the chosen experiment is explicitly labeled naïve CD4. Old RefSeq-based identifications are not treated as a new search of T29 or as evidence of protein absence.

## Frozen search and error controls

The [search protocol](../docs/ubash3a-proteomics-spectrum-protocol.md), [parameter file](../config/ubash3a-proteomics-comet.params) and [method freeze](../config/ubash3a-proteomics-search-plan.json) precede public-spectrum scoring. The database includes both T29 products, every earlier model, all retrieved human references and contaminants. Comet 2026.02.2 performs ordinary FASTA searches with internal concatenated target–decoy competition; it does not use a fragment index.

Settings include full trypsin, two missed cleavages, 7–35 residues, 600–5,000-Da peptides, 10 ppm precursor tolerance, isotope offsets 0/+1/+2, and high-resolution 0.02-Da fragment bins. Fixed carbamidomethyl C and variable oxidation M, phosphorylation STY, N-terminal Q cyclization and protein-N-terminal acetylation are included, with at most four variable modifications per peptide. The deposited experiment specifies 10 ppm whereas the general paper states 20 ppm. Comet binning and intensity processing differ from the old pipeline; this is a declared reanalysis, not an exact historical reproduction. [Comet's binning explanation](https://pmc.ncbi.nlm.nih.gov/articles/PMC4607604/), [decoy documentation](https://uwpr.github.io/Comet/parameters/parameters_202301/decoy_search.html).

Error estimates retain one winner per original file/scan, conservatively resolve reported target–decoy ties and collapse multiple charge queries. Best expectation values are ranked with tied thresholds, using `(decoy wins + 1) / target wins` and monotone q-values. Both spectrum-level and unique unmodified I/L-equivalent peptide summaries are retained. The primary background threshold is q ≤ 0.01.

Every fixed candidate is examined among all five exported assignment positions, regardless of its global q-value. Dense XCorr ranks can repeat across positions. A global 1% estimate does not establish 1% error for four overlapping novel targets. Candidate assignments require explicit fragment and competing-sequence review; novel-class error calibration remains limited. Normal/reference-compatible UBASH3A evidence is distinguished from assignments supported only by conditional protein models.

## Execution and technical qualification

The pinned Comet binary recovered seven generated control spectra with the expected peptides, modifications and target/decoy classes; all calculated masses agree within 0.001 Da and the expected fourteen or sixteen fragments were recovered. An out-of-range precursor was excluded. No public or T29 spectrum entered this qualification. These controls test software behavior, not experimental detection sensitivity. [Control record](ubash3a-proteomics-comet-qualification.json).

ThermoRawFileParser's selected build is a development release. All fourteen complete files were converted, yielding **66,766 MS2 spectra**, with **zero converter warnings**. The complete scan IDs/counts agree with the separate vendor-metadata export; independent XML checks pass for centroid arrays, precursor values, charge fields and instrument/activation annotations. All MS2 filters are high-resolution HCD on LTQ Orbitrap Velos. [Complete conversion record](ubash3a-proteomics-conversion.json).

Two upstream metadata-writer issues were identified in the pinned source: a complete-file count loop nested inside each scan iteration, and an MS1 display name copied onto the MS2 term. A separate one-pass reader exports the vendor scan counts and all MS2 identifiers; the independently parsed mzML must agree exactly. The source code, runtime and [technical repair](../config/ubash3a-proteomics-metadata-repair.json) are pinned. Both readers use the vendor library, so agreement is a computational consistency check rather than an independent instrument measurement.

The all-at-once UniProt export failed; complete paginated retrieval replaced it, with every failure and successful source receipt preserved. Slow raw-file transfers were overlapped with conversion and search only after each entire file passed its own integrity and conversion checks. The [execution-order amendment](../config/ubash3a-proteomics-execution-order.json) was recorded before any public spectrum inspection or score; it changed no selected file, search parameter or evidence rule. All fourteen searches completed before the [analysis start](../config/ubash3a-proteomics-analysis-start.json) and the first target-assignment inspection.

The [source manifest](../config/ubash3a-proteomics-sources.json) verifies all fourteen raw files and all 615,403 reference/model source records. It preserves 44 qualification-source receipts, including unsuccessful retrievals, plus the complete UniProt pagination record. All raw files match their listed total of 2,387,832,397 bytes and recorded local hashes. The [reproduction guide](../docs/ubash3a-proteomics-reproduction.md) distinguishes offline replay from acquiring a later live-service snapshot.

Comet's reported loaded count refers to its query-vector size. Its configured 5,000 batch threshold produced a 5,005-query batch when outstanding preprocessing completed; this is explicitly reconciled with the original protocol's batch wording in the [counting-scope record](../config/ubash3a-proteomics-counting-scope.json). Converted MS2 scans, loaded queries, exported queries and scans with winners must remain separate counts.

The initial summary parser stopped on repeated rank-one assignments before producing peptide summaries. The [documented repair](../config/ubash3a-proteomics-output-parser-repair.json) retains one deterministic, conservatively classified winner while preserving every tied assignment and its export position. It also resolves Comet's exported scan ordinal to the explicit original instrument ID. All 55,741 mappings agree with the complete mzML inventory, precursor values and retention times within output precision.

The first whole-output verifier then exposed missing decoy prefixes in XML alternative-protein names, including alternatives attached to a target parent. The [format repair](../config/ubash3a-proteomics-verifier-format-repair.json) uses the independently exported PIN memberships to recover those labels. Publisher source confirms the omission. Primary result bytes and all search parameters remain unchanged by that verifier correction; both failed checks and earlier code hashes remain preserved.

The [independent verifier](ubash3a-proteomics-search-verification.json) confirms all **615,403 source records**, all **308,191 distinct database sequences**, all **269,986 assignments**, 2,159,888 XML/text field comparisons and the PIN membership/modification/mass checks. It accounts for 3,479 queries with repeated rank one, 508 mixed target/decoy membership assignments and 831,997 omitted XML decoy prefixes. Both complete q-value calculations agree within 10⁻¹², including zero candidate assignments. These are computational consistency checks on one experiment, not biological replication.

All **238 repository tests pass**, including eight new output-format regressions. All **24 small primary outputs** replay byte for byte, and both large spectrum/peptide winner tables rebuild with identical hashes. [Replay verification](ubash3a-proteomics-replay-verification.json). The [stage verification](ubash3a-proteomics-stage-verification.json) records final hashes, test evidence and preservation of the earlier research.

## Interpretation boundary

Sequence specificity supports a measurable target definition. It does not establish that T29 exists. Even a matching ending peptide would require careful error control and would not by itself prove the complete RNA splice chain, allele identity, NMD escape or a role in PSC. A nondetection in this experiment would be limited by sample context, abundance, digestion, modifications, acquisition and these search settings; it would not prove biological absence.

The proposed laboratory RNA-identification/decay experiment remains unperformed. PRKD2 shared-signal analysis remains direction 4 of the roadmap.

# GSE255234 ETS2-study qualification

Audit date: **September 15, 2026**. This is source qualification, not a program analysis. Earlier plans and results remain unchanged.

## Decision

**No deposited input yet qualifies for the unchanged full-universe CPM procedure.** The workbook's historical annotation and complete count-denominator provenance remain unresolved. Public read metadata establish a possible reprocessing route, not a ready count matrix.

The **design identities of all three planned arms are now supported**: LPS−CON, H10−CON and H5−CON, each with four donor pairs. H5 dose/time is resolved by the original supplement. It must not be dropped merely because the earlier preparation lacked that evidence. Exact control harvest/vehicle details remain unreported, so an executable matched-control stratum is not yet certified.

This remains **Tier C primary-human macrophage context calibration**. It is neither direct ETS2 perturbation nor a PSC disease, treatment-response or clinical validation. Four donors and three related conditions are not twelve independent replicates.

## Design and control evidence

[Pradhan et al.](https://doi.org/10.1016/j.redox.2024.103191), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255234), and the [publisher supplement](https://ars.els-cdn.com/content/image/1-s2.0-S2213231724001691-mmc1.pdf) support:

| Arm | Dose and RNA-seq exposure | Donor pairs |
|---|---|---:|
| LPS−CON | LPS 1 µg/ml, alone, 6 h | 4 |
| H10−CON | Heme 10 µM, alone, 6 h | 4 |
| H5−CON | **Heme 5 µM, alone, 6 h** | 4 |

Supplementary Fig. S1, PDF page 2, and S2, page 4, explicitly state:

> hMDMs were treated with heme (5 µM, 10 µM) and LPS (1 µg/ml) alone for 6 h and then processed for bulk RNA-seq.

S2 also labels Con, LPS, Heme5 and Heme10. This resolves H5 using more than its name. The eight-page supplement was recovered from the publisher and independently from the Europe PMC supplementary archive; the two copies have identical bytes. They are copies of one source, not independent experiments.

GEO explicitly describes four individual donors. Each D1–D4 group has one CON/LPS/H10/H5 record. All 16 titles, sample keys and twelve planned pairs match the preserved preparation. No extra culture or technical split is documented. Methods 4.3 describes healthy-donor monocyte-derived macrophages, differentiated in RPMI with 5% AB serum and 25 ng/ml M-CSF. **Treatment** medium is 1% serum with 12.5 ng/ml M-CSF. These are different culture stages.

CON appears in the same experiment, supporting a shared source control. However, neither its exact harvest time nor solvent/vehicle equivalence is explicit. Heme stock preparation, culture/plate dates, randomization and library-batch assignments are also unavailable. Do not label these details verified or silently infer a human differentiation duration from the adjacent mouse protocol.

## Annotation and count-universe gate

The freshly retrieved [workbook](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255234/suppl/GSE255234_raw_counts_for_rna_seq_submission.xlsx) exactly matches the earlier 3,938,694-byte input and its SHA-256 `61d23d009b9b9a70a6127e0ae024df18f83fa3ed1a8079b7cbe86854ac5aac18`.

Structural inspection, without target-gene effects, finds:

- **64,253 unique unversioned ENSG identifiers × 16 columns**; no summary rows or version-removal collisions.
- All **1,028,048** entries are finite, nonnegative integers, with no missing values.
- **31,285 all-zero genes** and 674,699 zero entries are retained.
- Exact sample-key agreement; no duplicated sample keys or identical count columns.
- Positive deposited-row column totals, ranging from 24,203,383 to 60,454,799. These are not certified full-universe library denominators.

Methods 4.6 and every GEO sample instead describe **HT-seq v0.11.1, uniquely mapped reads, Hg38 RefSeq**. They provide no annotation file/release, Ensembl conversion, strandedness/counting arguments or export code. All deposited identifiers are Ensembl gene IDs; changing the prose label to Ensembl would be an unsupported correction.

S1's workflow separates raw counts from DESeq2 normalization and low-expression filtering. Its DE branch removes genes with normalized sum below ten across the two compared conditions; the main Methods says “all samples of interest.” The many retained zero genes show the workbook is not simply a table containing only genes that pass this filter. **They do not prove that every originally counted gene is present**, identify its annotation, or exclude another unreported export step. A row-count or identifier-set match to an unrelated reference release would not recover execution provenance.

Thus the [effective feature-universe gate](../config/ets2-validation-amended-rules.json) remains **failed/pending**, for every arm. Do not calculate the old full-universe CPM endpoint with this workbook, rename it “unfiltered,” or use only mapped signature genes in totals.

## Source QC and public reads

The paper names RNA-quality checks by gel or Bioanalyzer, TruSeq Stranded mRNA libraries and NextSeq 500 sequencing. GEO names FastQC 0.11.3, Trimmomatic 0.36, STAR 2.4.0.1 and RSeQC 2.6.4. No per-sample RNA/mapping QC results, acceptance thresholds, exclusion ledger or exact execution commands are released in the inspected linked sources. All sixteen submitted samples remain accounted for; no outcome-dependent exclusion was made.

[PRJNA1074058 / SRP488566](https://www.ebi.ac.uk/ena/browser/view/PRJNA1074058) has:

- **16 public runs**, SRR27896156–SRR27896171, one per GEO specimen; 16 distinct experiments/samples and no deposited run splits.
- **Paired-end** NextSeq 500 records and **32 generated FASTQ files**, with public URLs, source MD5s and byte sizes.
- **73,683,351,141 compressed FASTQ bytes** (73.68 GB) for all sixteen samples. The optional CON/LPS/H10 subset is 58,126,594,563 bytes; dropping H5 is no longer justified by missing dose/time.
- All 32 anonymous HTTP HEAD requests return 200 and matching Content-Length. **No FASTQ body was downloaded; source MD5s were not verified against read files.**

There is a new provenance conflict: paper/GEO says 75-base paired reads, but SRA reports mean mate lengths of **94.47–99.43 bases**. ENA's `read_count` equals SRA **spots/pairs** for these records, not the number of individual mates. Do not interpret twice that mean as a mate length or infer the cause of the disagreement. Public submitted reads are not proven to be untouched original sequencer output.

This establishes a bounded acquisition route only. A future reprocessing plan would need a pinned reference/annotation, explicit counting and strandedness procedure, read-level QC, complete-gene export and separate source-control decisions before program scoring. It would generate a **new** count input, not repair or reproduce an unknown historical pipeline. No processing was launched, and downstream scratch/reference storage is not included in the download total.

## A separate possible estimand—not a repaired CPM result

A **prospectively frozen available-universe rank calibration** could be a valid, narrow descriptive pilot: relative program position **within the 64,253 deposited features**, rather than expression per complete source count universe. Within-sample ranks do not need the unknown CPM denominator. This does not resolve annotation/export uncertainty; unknown feature selection can affect ranks and coverage.

Before such a pilot, freeze its rank/tie/zero rules, exact feature and program mapping, coverage/duplicate gates, control-only reference matching, all nine programs and all twelve source-paired comparisons. Address the control-time/vehicle limitation explicitly; do not claim a verified matched-vehicle effect. Preserve all arms and unavailable statuses. The result would describe an available-matrix experiment, not demonstrate ETS2 specificity or satisfy the old CPM protocol. **No rank protocol or execution freeze has been adopted here.**

**Most useful next action:** keep this small study as a documented reserve. If a GSE255234 pilot is selected later, decide prospectively between that limited rank estimand and a new, fully specified read-to-count pipeline. Do not start a 73.68-GB download merely to bypass unresolved provenance. Historical full-universe qualification still needs the original annotation and pre-filter export record; strict matching also needs control harvest/vehicle details.

## Provenance and verification

The [source ledger](../config/ets2-study-qualification-sources.json) pins URLs, UTC retrieval dates, bytes, SHA-256 values, failed requests, run/sample metadata and local audit outputs under `data/raw/ets2-study-qualification/`. Fresh GEO SOFT and article XML also match their earlier preparation hashes. Failed HTML/PDF routes are excluded from source evidence.

The finite source boundary was GEO metadata plus its sole supplementary workbook, the primary paper and sole linked supplement, SRA/ENA metadata, and accession searches in GitHub repositories, Zenodo and Europe PMC. GitHub/Zenodo returned no records; Europe PMC returned this paper. A Zenodo page-size error was retried within the public limit. Web-search service access was unavailable. This is not an exhaustive code search or a claim that unretrieved records do not exist.

An independent ZIP/XML workbook reader reproduces identifiers, all numerical structural checks, zero counts and every deposited-row column sum. Its initial header check needed a parser repair for an omitted blank A1 XML cell; no source or criterion changed. Reproduce with the project environment:

```bash
.venv/bin/python scripts/audit_ets2_study_qualification.py
.venv/bin/python scripts/verify_ets2_study_qualification.py
```

Python 3.12.14 and openpyxl 3.1.5 were used. Prior rules, preparation and benchmark artifacts retain their hashes. No program scores, new condition effects, model requests, read processing, contact, spending or clinical conclusions occurred. Published narrative and figure captions were incidentally visible during source review; this was not a blinded audit and those results did not determine selection.

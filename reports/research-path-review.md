# Research-path review: where more compute can change the PSC evidence

Reviewed 15 September 2026. This is a critical portfolio assessment, not a new biological result. It starts from `docs/evidence-log.md`, `ROADMAP.md`, `README.md` and the reports cited below. Later in this review, the parent authorized a bounded public MacroMap metadata qualification. No target expression effects, new predictions or clinical outcomes were inspected. No outreach or spending occurred.

## Main decision

**Prefer a genuinely separate PSC cohort and explicit disease comparators over another regulatory-model audit.** The newly proposed [GSE177044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177044) → [GSE119600](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119600) whole-blood route would be the first choice if its sample map, assay and independence checks qualify. Its aim should be external reproducibility and comparisons with other diseases, not intrinsic PSC specificity or a new clinical diagnostic claim. The blood-route source assessment is being performed separately; the sample numbers and access claims should come from that assessment rather than this review.

**MacroMap is the strongest documented next investment for larger-scale ETS2 context calibration.** It is a useful secondary route, not PSC disease validation. GSE255234 remains a small primary-human-cell fallback. Do not automatically escalate its annotation mismatch into full read reprocessing before trying the much larger, publicly released MacroMap expression matrix.

No currently qualified computational route guarantees a new PSC mechanism. The defensible near-term result is a measured, independently tested statement about the specificity or limits of a PSC-linked expression program. A reproducible nonspecific response is a substantive result. A new readiness report is not that result.

## What is strong already

- **Excellent source accountability.** Alleles, genome builds, source fields, uncertainty and missing evidence stay separate. The corrected UBASH3A coordinate and conflicting published alleles were not hidden.
- **Real negative tests.** The matched model comparison reports one anchor above and two below their comparison medians, with two evaluable regions and no accuracy claim. It did not replace unmatched BACH2 anchors after seeing effects. [Matched comparison](matched-comparison.md).
- **Biologically useful distinctions.** PRKD2 has allele-aligned RNA associations; its shared PSC cause is still unresolved. UBASH3A gene expression, a particular splice junction, a complete RNA, decay and protein production are kept distinct. [PRKD2 recovery](prkd2-public-data-recovery.md), [RNA-read audit](ubash3a-raw-read-audit.md).
- **Donors, not cells, are recognized as replicates.** The liver atlas and ETS2 benchmark preserve chemistry confounding, limited controls and repeated modalities. [Atlas](liver-atlas-cell-context.md), [program benchmark](ets2-program-benchmark.md).
- **Reusable analysis infrastructure.** Exact variant matching, donor aggregation, fixed gene programs and source-specific processing are now available. These reduce the cost of a well-chosen new hypothesis test.

These are strengths of the research, but test counts, archive bytes and repeated numerical replays are not additional biological replications.

## Recurring dead ends

1. **More precision about the same insufficient measurement.** GEUVADIS had only 18 eligible European donors and no CC donors. The later long-read experiment had no exact +29 reads; the protein experiment did not detect ordinary UBASH3A either. More conditional structures or unqualified spectra do not close that measurement gap. The +29 event itself was already published. [Junction follow-up](ubash3a-junction-followup.md), [raw reads](ubash3a-raw-read-audit.md), [protein search](ubash3a-proteomics-analysis.md).
2. **Repeatedly changing public sources without changing biological independence.** DICE subsets share donors; original and reprocessed data are not replications. Seven of eight selected PRKD2 ENCODE experiments share donors with training sources. ALL_FOLDS is not held-out-locus validation. [CD4 evidence](cd4-rna-evidence.md), [PRKD2 mechanism](prkd2-regulatory-mechanism.md).
3. **A larger calculation cannot identify missing covariance or recover missing alleles.** PRKD2 marginal coverage improved, but exact study model/covariance or full components remain missing. PFKFB3/BCL2L11 still have material coverage limitations. A converged external-LD fit is not a repaired study model. [IBDverse model audit](ibdverse-model-audit.md), [coverage repair](variant-coverage-repair.md).
4. **Scale can disguise a confounded design.** Whole-cell PSC/control chemistry is completely separated. Nuclear myeloid/CD4 controls number only two eligible donors. Batch integration cannot identify a disease coefficient when the design lacks overlap. [Atlas](liver-atlas-cell-context.md).
5. **Audit gates can become broader than the question.** Unknown historical software details can prevent exact reproduction without preventing every new analysis. Conversely, unknown donor pairing prevents donor-paired inference, and unknown assay sensitivity prevents strong absence claims. Separate these cases instead of treating every discrepancy as equally fatal.

Prospective freezes should prevent undisclosed selection, not prevent a clearly labeled new estimand. Preserve the old rules and unavailable outcomes. A new analysis can be narrower without pretending those old gates passed.

## Ranked paths

These ratings are qualitative judgments of the proposed work, not probabilities of success. Independence means a separate biological data source, not a second script.

| Rank / role | Concrete path | Scientific value | Feasibility | Independence | Potential new PSC insight |
|---|---|---|---|---|---|
| **1, conditional main project** | GSE177044 discovery/estimation followed by untouched adult GSE119600 replication and PBC/UC/CD specificity tests | High: directly tests disease reproducibility and specificity | Promising if sample/assay checks pass; modest expression-matrix computation | Potentially strong; cohort/donor overlap must be checked | Best near-term prospect; observational blood biology, not a liver mechanism or validated diagnostic |
| **2, context project** | MacroMap donor-held-out, multi-stimulus/time response of the fixed ETS2 and broad comparator programs | High for program transfer, precision and heterogeneity | Small enough public matrix for CPU work; metadata qualification below | Separate study; independent donors and source reuse need explicit support | Moderate mechanistic context; low direct disease specificity because samples are iPSC-derived macrophages, not PSC patients |
| **3, fastest small calibration** | GSE255234 LPS/control and H10/control donor program contrasts, with H5 explicitly unavailable if its background/time remains unresolved | Moderate; can show generic inflammation/stress rather than specificity | Small deposited workbook; annotation/denominator scope must be settled | Four explicit donor pairs in a separate study | Limited: four donors and no PSC comparison |
| **4, narrow descriptive alternative** | GSE84161 as-deposited arm-mean program shifts in G-573/G-432/LPS conditions | Moderate pharmacological context | The normalization method is complete; missing local caches may need restoration | Separate study, but anonymous cross-arm donor pairing unresolved | Limited; neither direct ETS2 perturbation nor PSC response |
| **5, high-value but input-gated** | PRKD2 qualified conditional/shared-signal test in original IBDverse monocytes | High for a susceptibility mechanism if adequate inputs arrive | Low now: compute cannot supply missing study covariance/model records | New cohort relative to BLUEPRINT/DICE; monocyte contexts still share participants | Potentially high, but no present shared-cause conclusion |
| **6, bounded measurement search** | UBASH3A exact +29-to-downstream RNA linkage in informative independent CD4 donors | High only if informative original reads exist | Low expected yield after prior sparse/negative searches; prequalify coverage | Further libraries from the same TRAILS donor do not replicate donors | Conditional RNA identity could be new; the +29 splice event is not new, and identity does not establish decay or disease causality |
| **7, exploratory disease-localization** | Within-donor PSC epithelial–myeloid/spatial niche contrasts | Interesting if genuinely new to the source study | Small donor counts; spatial panel/coordinate joins need qualification | Existing atlas modalities and sections overlap donors | Case-study or small-donor hypothesis, not robust PSC-specific population inference |
| **Do not prioritize** | Larger AlphaGenome sweeps, extra folding, pooled-cell PSC/control differential expression, or more colocalisation fits on unchanged inputs | Low incremental value under current gaps | Technically possible | Often repeats the same loci/sources | Output volume without the missing independent evidence |

### Important spatial limit

The existing liver collection metadata enumerate C73_A1/B1/C1/D1 and PSC011_A1/B1/C1/D1: **one control donor and one PSC donor**, not eight biological replicates. PSC011 is also in the selected nuclear atlas. See [liver-atlas-metadata-audit.json](liver-atlas-metadata-audit.json), its eight Visium object entries, and [library metadata](../data/derived/liver-atlas-libraries.csv). Four sections can support localization within a liver, not population PSC inference.

The four nested CosMx count archives in the [ETS2 inventory](../data/derived/ets2-benchmark-inventory.csv) remain uninspected. Their names do not establish donor counts, usable coordinates or a measured gene panel. No spatial cell-communication or PSC-specific effect should be promised from those names. The primary-cell atlas has only two PBC donors per modality, three unique overall, so that tempting disease-control comparison also remains very small.

## A completion-oriented main analysis

If the two blood cohorts qualify, freeze a compact question: **Does a fixed expression program reproduce in PSC, and is its shift distinguishable from PBC and IBD rather than generic systemic inflammation?** Existing published IFIT1/G0S2 results are benchmarks, not discoveries.

1. Establish unique participants, repeated visits, within-plate PSC/control overlap, assay type and a common measurable feature universe. Metadata and disease labels may be inspected for design; external target expression effects remain held out.
2. Use only discovery data to select any additional features, tune models or set thresholds. Start with fixed programs and a simple regularized linear baseline, not a broad machine-learning sweep. Split by donor; test processing-batch transfer where the design permits it.
3. Fix cohort-specific normalization and cross-platform score construction before external evaluation. Do not use external outcome labels for batch correction, feature selection or recalibration.
4. Evaluate the adult external PSC/control comparison, then every prespecified PSC/PBC/UC/CD specificity contrast. Shared controls and related tests are not independent validations. A direct between-disease contrast is needed; significance in PSC and nonsignificance in another disease is not a specificity test.
5. Report effect sizes and uncertainty, donor accounting, age/sex adjustment and medication/IBD limitations. Total-blood association and blood-cell-composition-adjusted association are different estimands; use the latter as a declared sensitivity where justified. Neither is a liver-cell causal effect.
6. Stop with the complete external result, including failed transfer or nonspecificity. Do not reopen the external cohort to optimize a better-looking signature. Clinical diagnostic claims would require a different validation and intended-use design.

The discovery cohort's UC samples must not support a clean PSC/UC inference if UC is restricted to separate plates. The independently measured disease groups are important precisely because more statistical adjustment cannot create missing design overlap.

## Useful narrower analyses, without an audit spiral

- **GSE255234:** a separately frozen deposited-feature-universe, within-sample-rank analysis could avoid claiming whole-library CPM or reconstruction of the undocumented RefSeq pipeline. This is a different, weaker estimand and still requires verified feature identity, controls and treatment time. Keep all arms/statuses and fixed programs. Four donors limit precision: an all-concordant two-sided sign test has minimum P=0.125, so present donor estimates rather than manufactured gene-level replication. [Readiness report](ets2-validation-readiness.md), [context search](../docs/ets2-validation-context-search.md).
- **GSE84161:** unknown pairing does not block an array-weighted mean(treatment) minus mean(control). With equal arms, this point contrast is invariant to a bijection between arms. It does block the original donor medians, paired uncertainty and donor leave-outs. A new descriptive plan could retain every as-deposited arm and eight mapping-qualified programs, omit donor-based P values, keep gRNA1-up unavailable, and assign neither disputed dose as correct. It would not pass the old paired protocol. [Normalization](gse84161-normalization-and-mapping.md), [supplement](gse84161-supplement-followup.md).
- **PRKD2:** separate essential statistical inputs from optional exact historical reproduction details, but do not discard the real missing-covariance, sample-scale or allele requirements. A reference-LD sensitivity can study assumptions under a new plan; it cannot certify a shared study signal. No additional model call repairs this.

## Resource and stopping policy

Use CPU-based expression-matrix work first. No GPU is justified by these proposed analyses. Raw read alignment should follow a demonstrated scientific need, not merely an inconvenient historical annotation description. Published pipelines and source checks exist, but several ignored raw caches and historical runtimes are absent in this checkout; a report of past normalization is not proof that its matrix is locally available now.

For each project, allow one bounded input-qualification pass. Before work, name the effect estimate the project will deliver, its independent replication unit and the data that would make the estimate unavailable. Keep all prespecified outcomes. Once the effect panel is complete, stop even if it is weak or opposing. Park input-gated branches until a genuinely new input appears rather than restarting the same search with a larger transfer budget.

## Completed bounded MacroMap qualification

This section records new **input qualification**, completed during the review after authorization. It does not report gene-program effects. The public expression matrix and sample metadata now support preparing a concrete context analysis; no scoring plan has been executed.

### What was verified

- [Zenodo 11563707](https://zenodo.org/records/11563707) lists the raw, TPM and filtered-TPM files separately, with open CC BY 4.0 access. The raw file is **174,122,135 bytes**; its publisher MD5 `3e585dfeb097368db7eaab9889f76f58` matches the complete download. SHA-256 is `eb9dd32b6f25ae9912e47dcc5e1f5170a98af919a42d48a4bca3a680fc5353d9`.
- A bounded 262,144-byte range first established the full header. The source is concatenated gzip/BGZF; a first-member-only decoder truncates that long header. The standard gzip reader recovered the complete line before any whole-file download. Both the initial truncated prefix and corrected header are preserved.
- The full raw file passes gzip EOF/CRC checks and contains **58,243 unique gene rows × 4,698 unique sample columns**. All **273,625,614** entries are finite, nonnegative and numerically integral. There are no duplicate gene IDs or malformed rows, 2,397 all-zero gene rows remain, and every library total is positive. These were whole-file structural checks, not target expression inspection.
- The expression paper is [Panousis et al., DOI 10.1038/s41467-025-61670-9, PMC12391345](https://pmc.ncbi.nlm.nih.gov/articles/PMC12391345/). The earlier repository context link, **PMC12391537**, is the companion splicing-QTL paper. The historical reference was not overwritten. The expression methods specify GRCh38, GENCODE v27, STAR 2.5.3a and featureCounts 1.5.3, with 75-bp paired-end RNA sequencing. Raw sequencing has a controlled-access route; the released count matrix does not require that access.
- Every released versioned gene ID matches the original repository's pinned GENCODE v27 annotation. That annotation has 58,288 records; the **45 reference-only IDs all carry `_PAR_Y`**. There are no matrix-only IDs. Preserve those absent reference records rather than imputing zeros or removing suffixes to collapse genes. The exact released gene universe is therefore known. This does not assert that every historically quantified feature was exported. The release is demonstrably broader than the paper's 14,060-gene filtered analysis universe.
- Publisher Supplementary Data 11 contains a 4,698-row, 21-column RNA-seq metadata sheet. **All 4,698** matrix names equal `SampleID + "_" + RunID`; there are no unmatched or extra records. There are 209 source-line keys. Of these, 205 have distinct mapped HipSci IDs; four separately named lines, covering 91 samples, carry literal `NOMATCH`. They must never be collapsed into one person.
- The paper states that its 217 starting iPSC lines came from **unrelated healthy HipSci donors**; 209 source lines remain in the release. This is a published independence statement, not a new kinship check. There are **124 modified-SmartSeq2 lines** and **85 NEB lines**, with 2,788 and 1,910 samples, respectively, across 57 sequencing runs. Each line uses one run and one protocol. Its recorded culture/differentiation fields remain constant across conditions.
- All 20 stimulation/time labels have a same-time basal control and **165–197 structurally matched line pairs**. Excluding the four `NOMATCH` lines leaves **162–193**. Every available pair matches recorded run, protocol and static culture fields. This reduces the protocol problem for within-line contrasts, but does not independently identify a protocol effect because different donor pools used the two protocols.
- **PIC is an important exception to biological control matching.** Its methods use lipofection and mention a mock-transfection control, but the release has no separately identified mock arm. PIC minus basal Ctrl cannot isolate poly(I:C) from delivery. A strict matched-background primary panel should retain its two comparisons as unavailable; a separately labeled delivery-plus-PIC comparison would be a different estimand. The other 18 labels have **169–193 mapped-line pairs**. CIL and LIL10 are mixtures, not isolated-component experiments. MBP is myelin basic protein, not a mock control.

No MacroMap/HipSci/Panousis or release-accession reuse was identified in the inspected Stankey main article; its known GSE47189 reuse remains visible. This is a bounded source-selection check, not proof that every donor or training source is disjoint. MacroMap's own published QC used expression patterns when excluding failed/mislabeled stimulation samples; the released cohort is not a blinded test of stimulation success.

### Practical analysis freeze

A new plan could now fix the following without another open-ended source search:

1. **Question:** does the fixed ETS2 program add reproducible stimulus-context information beyond broad inflammation, interferon, oxidative-stress and apoptosis programs? This is not an ETS2-dependence or PSC-diagnosis test.
2. **Units and scope:** use the 205 HipSci-mapped source lines as a conservative primary identity set and all 209 as an identity sensitivity. Keep all conditions of each line together. Retain the complete 20-row stimulation/time registry, with the two PIC rows unavailable for the strict background-matched endpoint. Precursors are separate maturation comparisons, not stimulation controls.
3. **Holdout:** split whole sequencing runs within protocol strata before reading target effects. This also keeps every line out of both training and test sets. Use training controls only for reference-gene bins/draws and any score scaling or model tuning. The exact split, minimum eligible groups and seed must be recorded in the execution plan.
4. **Primary measured result:** compute the unchanged program and comparator within-line stimulation-minus-time-matched-control scores, retain all line estimates, and summarize effects separately by library protocol. Preserve the nonoverlap ETS2 sensitivity. Use donors/lines or run blocks—not genes, samples or stimuli—as the resampling unit, according to the declared uncertainty model. Do not tune gene-set membership against favorable conditions.
5. **Optional new predictive endpoint:** compare a fixed broad-program response baseline with that baseline plus the ETS2 response score for predicting prespecified stimulus labels in held-out lines/runs. Predeclare multiclass log loss and equal line weighting. This tests extra *condition information*, not whether the residual is caused by ETS2. It is secondary to the complete response-effect panel and must not become an unlimited model search.
6. **Final assay gates:** map the unchanged programs to the verified released universe under fixed identifier rules; preserve coverage failures. Pin the released-universe count normalization and reference matching rather than claiming exact reconstruction of the authors' historical filtered analysis. Freeze the complete executable plan before scores. The existing historical raw-count gate is not silently amended by this recommendation.

The result should be one complete response panel, protocol-stratified uncertainty and the prespecified held-out comparison if included. A weak or redundant program is a valid outcome. This is substantially more informative than re-aligning four-donor GSE255234 reads solely to recover historical annotation details, while remaining less directly relevant to PSC disease association than the two-cohort blood project.

### Source and replay record

Detailed inputs and failed/successful retrieval receipts are in ignored `data/raw/macromap-qualification/`. Downloaded response bodies total **185,995,664 bytes**, below the 250-MB ceiling. The full matrix was downloaded only after its bounded header and source sample design qualified. No dependency was installed; the unavailable RDS decoder was bypassed using the original XLSX metadata.

Key local records:

- `metadata-ledger.jsonl` and `worker-methods/source-ledger.json`: source URLs, retrieval timestamps, response receipts, sizes and hashes.
- `matrix-download-receipt.json`, `matrix-structural-audit.json`, `gene-universe-audit.json`: complete-file and released-feature checks.
- `worker-methods/expression-header-sample-map.csv`, `worker-methods/condition-control-pair-counts.csv`, `worker-methods/methods-qualification-summary.json`: exact source sample joins, paired-design accounting, protocol and control caveats.
- `paired-design-summary.csv`: a separate metadata aggregation agrees on all 20 total paired counts and recorded culture/run/protocol matching.
- `qualify_matrix.py`: native project-environment whole-file acquisition and structural check; it emits no gene-specific expression results or program scores.

The public source pins most useful for restoring this qualification are:

| Source | Bytes | SHA-256 |
|---|---:|---|
| [Expression article XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12391345/fullTextXML) | 137,743 | `d3f69f885db807ff964f9158c7b7f581bea649225d6c76a37e4d3ef34f38cb8f` |
| [Supplementary methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-61670-9/MediaObjects/41467_2025_61670_MOESM1_ESM.pdf) | 7,140,183 | `00d11361e7a874452c84ac062409eb95667482c4a58e0681758a674283452951` |
| [Supplementary Data 11](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-61670-9/MediaObjects/41467_2025_61670_MOESM13_ESM.xlsx) | 526,439 | `97805429aaf777d8d13f625375a0dee124f130a2499fd6237bec279295bd89d8` |
| [Original-study GENCODE v27 annotation](https://raw.githubusercontent.com/andersonlab/macromap_eqtl/191fc67f84d631e5b4fc9a97993c142b348cf2c9/Data/gencode.v27.annotation_chr_pos_strand_geneid_gene_name.gtf) | 3,079,419 | `db11e855d959656bc3a25de65460ae12368ca13b1ae1a5493973f6ea433bdf41` |

GitHub was pinned to tree `191fc67f84d631e5b4fc9a97993c142b348cf2c9`. All retrievals occurred on 15 September 2026 UTC. Source numeric expression values were used only for file-level validity checks. No target-gene values, program scores, model effects, shared-cause posterior or clinical result was produced. No prior research file or frozen protocol was changed by this review.

# GSE84161 publication and inhibitor provenance

**Search date:** 2026-09-12. **Stopping point:** focused provenance pass completed after the direct GEO family records, the associated primary paper, its cited inhibitor papers, and one author-hosted full-text copy were checked. This is not an exhaustive literature or compound-identity census. No expression values, gene-level effects, signatures, or scores were inspected.

## Finding

GSE84161 is the primary human monocyte microarray subseries of GSE84162, the Senger et al. study reported as **“The kinase TPL2 activates ERK and p38 signaling to promote neutrophilic inflammation”**, *Science Signaling* 10(475):eaah4273 (18 April 2017), PMID [28420753](https://pubmed.ncbi.nlm.nih.gov/28420753/), DOI [10.1126/scisignal.aah4273](https://doi.org/10.1126/scisignal.aah4273). PubMed records the title, date, DOI, the first authors Kate Senger, Victoria C. Pham, Eugene Varfolomeev, and Jason A. Hackney, and Genentech Research affiliations for the author list. GSE84162 explicitly lists GSE84161 as a SubSeries and GSE84153 as the related human-neutrophil RNA-seq subseries.

This is a Tier B biological candidate for the PSC project: primary human monocytes, five stated buffy-coat donors, MEK and TPL2 inhibitor arms, LPS stimulation, and matched vehicle groups. It is an Affymetrix GPL570 array experiment rather than a public raw integer gene-count experiment. It therefore remains a separate array-calibration lead and is not ready for the frozen count-based validation procedure. The accession is the original source for Senger et al., so analysis of it does not independently replicate that paper. It is separate from the Stankey ETS2 source by accession and contributor group; public records do not establish individual donor overlap. No ETS2-program response has been measured here.

## Evidence ledger

| Question | Evidence status | Record-backed conclusion |
|---|---|---|
| Associated publication | **Explicit** | Senger et al., *Science Signaling* 2017, PMID 28420753, DOI 10.1126/scisignal.aah4273. |
| Series relationship | **Explicit** | GSE84161 is a SubSeries of GSE84162; GSE84162 also contains GSE84153. |
| Cell system | **Explicit** | Human monocytes purified from buffy coats by negative selection and cultured with M-CSF. GEO states five total donors. |
| Perturbations | **Explicit labels; chemical names absent** | GEO calls the arms “Tpl-2 SMI”/`TPL2i` and “MEK SMI”/`MEKi`, with vehicle controls and LPS arms. |
| Donor pairing | **Paired design supported; donor key unavailable** | Replicate numbers 1–5 each contain all six treatment/stimulation labels. The public sample records do not provide donor IDs or a separate replicate-to-donor key. |
| Exact TPL2i identity | **Unresolved** | No chemical name, catalog number, CAS, or structure was found in the accessible GEO metadata or article text. Do not substitute a Hall et al. compound or a commercial TPL2 inhibitor. |
| Exact MEKi identity | **Unresolved** | The cited Hatzivassiliou et al. (Nature 2013) study discusses several MEK inhibitors. A literature citation alone does not identify the GSE84161 compound; no individual candidate is assigned. |
| Pathway engagement | **Published qualitative evidence** | The paper reports that, in primary human monocytes, both inhibitors reduced LPS-associated pERK1/2 while p38 and MEK3/6 phosphorylation did not change in the monocyte assay. GEO says supernatants were saved for ELISA, but no ELISA result file is present in the inspected GEO metadata. |
| Original analysis code | **Unavailable in inspected public records** | Methods describe custom R/Bioconductor scripts, but no public repository or script accession was located in the paper/GEO records. |

## Sample structure and pairing

The GSE84161 series design states that monocytes were purified from five buffy-coat donors, treated for one hour with TPL2 or MEK inhibitor, stimulated or not with 10 ng/ml LPS, and harvested after six hours. The sample-level SOFT records provide six complete labels for each replicate-number group:

| Replicate label | Untreated | Tpl-2 SMI | MEK SMI | LPS | Tpl-2 SMI + LPS | MEK SMI + LPS |
|---:|---|---|---|---|---|---|
| 1 | GSM2228386 | GSM2228387 | GSM2228388 | GSM2228389 | GSM2228390 | GSM2228391 |
| 2 | GSM2228392 | GSM2228393 | GSM2228394 | GSM2228395 | GSM2228396 | GSM2228397 |
| 3 | GSM2228398 | GSM2228399 | GSM2228400 | GSM2228401 | GSM2228402 | GSM2228403 |
| 4 | GSM2228404 | GSM2228405 | GSM2228406 | GSM2228407 | GSM2228408 | GSM2228409 |
| 5 | GSM2228410 | GSM2228411 | GSM2228412 | GSM2228413 | GSM2228414 | GSM2228415 |

The table is a metadata map derived from the published sample titles and `replicate number` fields. The five groups plausibly correspond to donors because the series design says five donors and the source analysis adjusts for donor, but the correspondence is not explicit in the public sample annotations. Treat them as five candidate blocks. Confirmation only requires the sample-to-donor grouping; personal identities are neither needed nor requested.

The measurement is expression profiling by array on GPL570 (Affymetrix Human Genome U133 Plus 2.0). GEO lists a 30-sample `GSE84161_RAW.tar` containing CEL files and a processed series matrix of RMA intensity. The raw archive is an array file set; it does not provide an integer gene-count matrix. A separate assay-specific calibration and a donor-key decision would be required before any numerical program comparison.

## Compound provenance and concentration discrepancy

The accessible primary records support the following hierarchy:

1. **GEO explicit record.** The overall design gives **2.8 µM** for the Tpl-2 SMI and **0.31 µM** for the MEK SMI. Sample characteristics name the compounds only as `TPL2i` and `MEKi`.
2. **Senger article explicit record.** Figure 2D describes the human monocyte microarray as preincubated with **926 nM TPL2i** versus **103 nM MEKi**, followed by six hours of LPS stimulation in five healthy donors. The caption says the concentrations were selected to fully block ERK phosphorylation. These values do not equal the GEO series-design values (2.8 and 0.31 µM). The inspected sources do not resolve whether they represent a different assay subset, an earlier protocol, or a metadata discrepancy.
3. **Cited inhibitor literature.** The article cites Hall et al. (JBC 2007, DOI [10.1074/jbc.M703694200](https://doi.org/10.1074/jbc.M703694200)) for pharmacological TPL2 work and Hatzivassiliou et al. (Nature 2013, DOI [10.1038/nature12441](https://doi.org/10.1038/nature12441)) for MEK-inhibitor mechanism. Hall et al. describes named historical TPL2 compounds, and Hatzivassiliou et al. discusses several allosteric MEK inhibitors including PD0325901. Neither citation establishes that the same molecule was used in GSE84161.

Accordingly, the exact TPL2i and MEKi identities are **unknown**. The unavailable publisher supplement, which the article identifies as containing inhibitor selectivity and supporting figures, is the most useful unresolved source for chemical provenance. No compound identity should be backfilled from a structurally related or commercially available inhibitor. This limits chemical attribution; it does not prevent fresh array preparation using the source's arm labels.

## Engagement and computational provenance

The paper reports target-engagement experiments in primary human monocytes and neutrophils. For monocytes, the accessible article text states that TPL2i and MEKi reduced pERK1/2 after LPS stimulation, whereas p38 and MEK3/6 phosphorylation were not changed in the corresponding monocyte assay (Fig. S9). The article also points to Fig. S2E for the ERK-phosphorylation potency basis used to choose the inhibitor concentrations. This supports pathway engagement as published assay evidence; it is not a replacement for a donor-level transcriptomic readout.

GEO's processing statement says the array data were RMA normalized, one probe set per Entrez gene was retained by highest variance, unmapped probes were removed, and limma models included LPS, treatment cohorts, and donor as a covariate. The paper's methods describe custom R/Bioconductor analysis, RMA and limma for monocyte arrays, and donor-only residual normalization for heat-map visualization. The RNA-seq part of the umbrella paper names HTSeqGenie 3.16.1, GSNAP, and DESeq2, but that pipeline applies to GSE84153 neutrophil RNA-seq rather than the GSE84161 monocyte array.

No original R scripts, package lockfile, or code repository was found in the inspected paper, supplement endpoint, or GEO family metadata. The public record therefore supports method-level description but not exact code replay.

## Bounded search log and unresolved items

The focused pass used one direct record each for GSE84161, GSE84162, PubMed PMID 28420753, and the Crossref DOI record; one page of the author-hosted ResearchGate copy; and one citation cross-check for Hatzivassiliou et al. Search/find checks on the author-hosted article used `TPL2i`, `MEKi`, `PD0325901`, `GDC-0973`, `GDC-0623`, `compound`, `Table S1`, and `Reagents`. No multi-page compound census or broad new dataset search was run after the associated paper was identified.

The Science publisher article and supplement endpoints returned Cloudflare challenge HTML in the captured files, so those files are retained as failed retrievals rather than treated as article or supplement content. The ResearchGate page exposes indexed primary-paper text, including Methods and figure-caption material, but the attempted local download was HTML rather than a valid PDF. PubMed XML, Crossref JSON, GEO SOFT metadata, and the GEO CEL archive are the usable local provenance caches. The exact inhibitor identities, public donor key, and original analysis code remain unresolved.

### Sources

- [GEO GSE84161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84161) and [GEO GSE84162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84162)
- [PubMed 28420753](https://pubmed.ncbi.nlm.nih.gov/28420753/)
- [Senger et al. author-hosted article copy](https://www.researchgate.net/publication/316231923_The_kinase_TPL2_activates_ERK_and_p38_signaling_to_promote_neutrophilic_inflammation) (secondary host of the primary article; indexed full text used only for method and caption provenance)
- [Hall et al. JBC 2007](https://pubmed.ncbi.nlm.nih.gov/17848581/) and [Hatzivassiliou et al. Nature 2013](https://pubmed.ncbi.nlm.nih.gov/23934108/)

Local source sizes, hashes, failed retrievals, and query bounds are recorded in [`config/gse84161-literature-sources.json`](../config/gse84161-literature-sources.json).

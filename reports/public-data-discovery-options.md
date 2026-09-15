# Public-data alternatives for PSC research

**Source check: September 15, 2026 (UTC).** This is a discovery and input-qualification report, not a biological result. Read alongside [the evidence log](../docs/evidence-log.md) and [ROADMAP](../ROADMAP.md).

## Recommendation

**Run a fixed, disease-specific blood replication analysis first. Follow it with a donor-paired cholangiocyte perturbation analysis.** Both answer questions that the existing variant/model audits and generic macrophage calibration cannot answer.

| Rank | Question and source | Why it adds value | Readiness |
|---|---|---|---|
| **1** | Which PSC-associated blood programs replicate in a separately recruited cohort, and distinguish PSC from PBC, UC and CD? GSE177044 + GSE119600. | Direct disease comparisons, hundreds of participants, independent recruitment and explicit alternative-disease controls. Can separate reproducible PSC-associated signal from a generic inflammatory program. | **Complete expression and annotation inputs acquired and structurally checked.** Ready for a separately frozen analysis; substantial covariate limits remain. |
| **2** | Is the IL-17A response different in PSC versus non-PSC cholangiocyte organoids when the unit is the donor, not the cell? GSE239283. | Direct intervention in the relevant epithelial cell type; a formal difference-in-responses is stronger than comparing separate DEG counts. | Public integer-count matrix and complete eight-donor/16-library map verified. Full matrix not downloaded. |
| **3** | Does a predeclared epithelial–stromal injury program distinguish PSC from other liver diseases across pediatric biopsies and adult tissue? GSE303271 + GSE159676; limited spatial follow-up in GSE325861. | Direct liver tissue, separate patient collections, and a test of PSC specificity versus generic fibrosis. | Small public processed inputs verified by byte ranges. Weaker disease-control sizes, stage differences and sparse clinical covariates lower its priority. |

The first route has the best immediate compute-to-information ratio. The second has the strongest experimental mechanism leverage. Neither establishes a treatment recommendation. The third is preferable to another broad model-score sweep, but its confounding is more severe than the sample totals alone suggest.

**Compared with GSE255234:** the four-donor heme/LPS macrophage experiment can calibrate program behavior under generic stimulation, but it has no PSC case–control contrast. The blood route directly tests reproducibility and disease specificity in many more participants; the organoid route tests an intervention in PSC-derived cholangiocytes. These are stronger next questions for direct disease relevance. They do not invalidate the prepared ETS2 calibration, which remains useful for a narrower assay-context question.

## What was done before inspecting effects

The initial question record was saved at **2026-09-15 22:25:39 UTC**, before new expression values or target coefficients. It asked about PSC-specific cellular states, defined cholangiocyte perturbations, and separation from IBD or general liver injury. A subsequent source-selection record fixed the dataset-specific questions before matrix inspection. Search responses necessarily included published titles and study summaries; published findings are prior knowledge, not new discoveries.

This scout did not calculate differential expression, disease-group scores, genetic effects or model predictions. Complete numeric checks below concern file integrity, namespaces and global value distributions only. Candidate genes were not selected from the new matrices. A Methods-section text filter also returned a published organoid Results subsection after the organoid question had been recorded; those published results were not used to select a new gene list.

Serper was unavailable. Discovery used NCBI GEO/Entrez, Europe PMC, publisher-linked records and public source inventories. The search was finite, not an exhaustive statement about all available PSC data. No researcher was contacted, no account was created, and no paid or controlled-access route was used.

## 1. Disease-specific blood replication — recommended first

### Verified studies and access

**Norwegian PSC discovery arm:** [GSE177044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177044), Wacker et al., [paper and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10832281/), DOI [10.1016/j.jhepr.2023.100988](https://doi.org/10.1016/j.jhepr.2023.100988).

The full study has 1,035 deposited libraries. It contains **two case–control cohorts**, not a geographically balanced PSC-versus-UC cohort:

- Norway: 220 PSC participants, including 177 with UC and 43 without a concurrent UC diagnosis; 77 controls.
- Germany: 495 UC participants and 243 controls.
- The paper reports participant counts that match the deposited library count. This supports the reported participant-level design, but the anonymous labels do not independently verify person-level disjointness from every prior study. The exact matrix-to-GEO mapping is checked.
- The primary Norwegian file contains all **297** Norwegian records. Every record has age, `Sex` and `sequencing plate`; 214 PSC records have years since diagnosis.
- Across the full series, plate fields use two literal keys and sex/age completeness differs. Preserve the raw keys and record any normalization. A naïve pooled PSC-versus-UC classifier would confound diagnosis with recruitment cohort and plate structure.

Exact downloadable primary input:

- [GSE177044_raw_countsPSC.csv.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE177nnn/GSE177044/suppl/GSE177044_raw_countsPSC.csv.gz): **8,136,584 bytes**. Header: `Geneid,gene_name`, followed by 297 sample-title columns.
- Public [source analysis repository](https://github.com/ikmb/ucpsc-rnaseq) is available. Its tree was inspected, not its effect outputs.
- GEO supplies SRA links, but raw reads were neither required nor downloaded. These files are measured count matrices, not GWAS summaries or individual genotypes.

**Polish validation and disease controls:** [GSE119600](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119600), Ostrowski et al., [paper and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6510750/), DOI [10.1038/s41598-019-43699-1](https://doi.org/10.1038/s41598-019-43699-1).

All enrolled participants are described as Polish Caucasians. The deposited 370 arrays contain:

| GEO age-class label | PSC | Control | PBC | UC | CD | Total |
|---|---:|---:|---:|---:|---:|---:|
| Adult | 45 | 47 | 90 | 45 | 48 | **275** |
| Child | 0 | 0 | 0 | 48 | 47 | **95** |

Use the **275 source-labeled adult arrays** for the primary comparison. Do not pool pediatric IBD controls into an adult PSC analysis. “Adult” is the source classification, not a independently verified minimum age: the paper's PSC discovery age range starts at 16. GEO provides diagnosis and age class, but not sample-level exact age, sex, medication, leukocyte counts or PSC–IBD status. The paper describes all PBC participants as female and the PSC group as mixed sex. Sex and treatment confounding cannot be repaired by inventing sample covariates or inferring sex from a disease-selected expression signature.

Verified downloadable inputs:

- [GSE119600_series_matrix.txt.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119600/matrix/GSE119600_series_matrix.txt.gz): **69,056,692 bytes**, complete normalized expression matrix, **47,230 probes × 370 samples**.
- [GPL10558.annot.gz](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz): **7,290,886 bytes**, 48,107 unique platform rows. Every matrix probe has an annotation row; 31,265 matched rows have nonempty `Gene symbol` and `Gene ID`. The FTP annotation snapshot is dated August 22, 2016. It supplies source relationships, not a current reconciled gene identifier system.
- The [supplementary inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119600/suppl/) also lists a roughly 105-MiB non-normalized matrix. It was **not downloaded**.

**An important access correction:** the tempting 27,443,200-byte `GSE119600_RAW.tar` is **not subject expression data**. Its [file list](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119600/suppl/filelist.txt) contains two GPL10558 annotation files. A bounded original TAR-header request confirms the first member. The earlier provisional suggestion that this smaller TAR might supply array expression is superseded. No full TAR was downloaded.

### What “independent validation” means here

The Norway/Germany and Polish cohorts were separately recruited. However, **the GSE177044 paper already used the Polish study as a replication cohort**; its reference 20 is the GSE119600 paper. This is **independent-cohort reuse**, not untouched held-out validation. Repeating the original classifier or the source-selected IFIT1/G0S2 narrative is not a new mechanism discovery.

The proposed added value is a prospectively fixed, complete disease-specificity test with PBC, adult UC and adult CD controls, transparent missingness and all planned results. No claim is made that any resulting gene is newly discovered. The larger qPCR replication cohorts described in the Polish paper are not additional public whole-transcriptome matrices.

### Proposed donor-level design

1. Fix all tested genes/programs, mapping rules, contrasts and multiplicity before effects. Either test the full eligible gene universe with a fixed discovery-to-replication rule, or choose independently defined pathways before scoring. Do not build a signature on the Polish validation labels.
2. Use Norwegian PSC versus its 77 same-cohort controls as discovery. Model age, sex and sequencing plate where estimable. Keep PSC with and without concurrent UC separate as prespecified secondary contrasts. One participant is one replicate.
3. Use the Polish adult PSC-versus-control contrast for replication. Then report the same frozen measurements for PSC versus PBC, adult UC and adult CD. A replicated PSC-versus-healthy signal that is equally changed in every disease is not PSC-specific.
4. Keep the studies and platforms separate. Use gene-level effect direction or standardized pathway effects rather than pooling normalized matrices. Do not equate a sequencing log-fold change with an array standardized effect. Retain ambiguous and absent cross-platform mappings explicitly.
5. Report composition-sensitive and broader inflammatory programs alongside any candidate program. Deconvolution is an expression-based estimate, not measured cell counts or proof of a cell-intrinsic mechanism. Compare primary and composition-sensitive analyses without claiming residual confounding is eliminated.
6. Treat the Polish source's nine historical outlier exclusions and its PSC count discrepancies as historical context. The main paper's Table 1 has 46 PSC discovery participants whereas GEO has 45 PSC records. Do not guess which samples to remove; a new analysis can define its own blinded QC on the complete deposited set.

**Mechanistic interpretation:** this can identify a reproducible PSC-associated systemic program and challenge claims that it is merely a generic inflammatory response. It cannot by itself assign a cell of origin, establish genetic causality, or show whether the program causes PSC rather than follows disease, medication or cell-composition changes.

### Completed input qualification and compute estimate

All three full GSE177044 files were acquired under the subsequent explicit qualification authorization. They have **63,677 unique unversioned Ensembl gene IDs each**, finite nonnegative integer counts, and identical gene sets. The combined headers give 1,035 unique samples with an exact GEO-title match. The `gene_name` field has **7,039 duplicate rows beyond first occurrences** per file; symbols must not be blindly collapsed. The source genome is **GRCh37**, with STAR 2.6.1d, FeatureCounts 1.6.4 and nf-core/rnaseq 1.3 declared. Exact identifier-annotation version remains a source attribute to retain, not guess.

The Polish matrix passes complete gzip EOF, row-width, numeric and probe-ID checks. Its deposited values are positive, with global minimum **37.50227**, median **73.0333**, 99th percentile **6721.527**, and maximum **47771.99**. These are **linear-range normalized intensities**, not an already-log2 range. GEO states quantile normalization using `lumi`; the paper states quantile-normalized average bead signals with no background correction. Any new log2 transformation must be explicitly recorded as our analysis choice.

Expected compute: **4–8 CPU cores, 8–16 GB RAM, no GPU**; minutes to roughly two hours for ordinary gene/program models and modest permutation/sensitivity analyses, depending on the implementation. These are planning estimates, not measured performance. The primary three required compressed files total **84,484,162 bytes**; raw-read alignment is unnecessary.

## 2. Donor-paired PSC cholangiocyte IL-17A response

**Sources:** [GSE239283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239283); [Garcia Moreno et al. paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11150034/), DOI [10.1097/HC9.0000000000000454](https://doi.org/10.1097/HC9.0000000000000454).

### Access and biological units

The GEO design explicitly specifies **four PSC patients and four non-PSC patients**, each with vehicle and IL-17A organoid cultures: 16 libraries, not 16 independent donors. All are passage 3; treatment is **100 ng/mL recombinant human IL-17A for 24 hours**.

Both the GSM titles and the released cell metadata support these complete pairs:

- PSC: `PSC_4`, `PSC_6`, `PSC_8`, `PSC_9`.
- Non-PSC: `nonPSC_2`, `nonPSC_3`, `nonPSC_4`, `nonPSC_5`.
- Every key has `CNT` and `IL17`; the metadata contains **46,343 cells**, all assigned to those 16 source libraries and eight `Sample` labels.

The non-PSC donors are procedure controls with benign biliary/pancreatic indications, not healthy volunteers. The paper reports donor age, sex, sampling site and cirrhosis information; its larger overall clinical table must not be mistaken for the eight sequenced donors. Disease and control groups differ in age, clinical context and sampling method. Some source material is bile-derived rather than brush-derived.

Public inputs:

- [Integer matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239283/suppl/GSE239283_extrahepatic_organoids_matrix.mtx.gz): **142,381,971 compressed bytes**. Its MatrixMarket header specifies **21,562 features × 46,343 cells, 46,282,529 stored entries**. Header and dimensions only were read; whole-file integer validity has not yet been checked.
- [Cell metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239283/suppl/GSE239283_extrahepatic_organoids_metadata.tsv.gz): **2,239,262 bytes**, fully acquired. Despite its `.tsv.gz` name, the body is comma-separated. Original names are preserved.
- [Gene list](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239283/suppl/GSE239283_extrahepatic_organoids_genes.tsv.gz) and [barcodes](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239283/suppl/GSE239283_extrahepatic_organoids_barcodes.tsv.gz) are listed in the verified [directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239283/suppl/), about 132 and 214 KiB. Their full contents were not inspected.
- Original per-sample count-matrix, feature and barcode links are in every GSM record. The approximately 697-MiB series TAR was not downloaded. The merged matrix and metadata are processed measurements, not raw reads or an unprocessed historical pipeline reconstruction.
- GEO declares hg38, Cell Ranger 5.0.0, Seurat 4 and Harmony. Data availability guarantees the scRNA-seq deposit; code and additional data require an author request. Public donor-level NanoString/Olink inputs were not independently verified.

### New question, not a replay of DEG counts

Use a **donor × treatment pseudobulk model**. The main estimand is:

`mean(IL17 − vehicle in PSC donors) − mean(IL17 − vehicle in non-PSC donors)`.

Use a paired count model with donor blocking and a disease-by-treatment interaction, or fit donor-level changes and compare the four-versus-four donor groups. Do not include a redundant disease main effect in a donor-fixed-effects design. Use a fixed all-epithelial primary population; author-defined cell clusters are secondary, because treatment-dependent clustering can select an outcome-defined population. Report every donor's change, fixed broad inflammatory/stress comparators, and leave-one-donor-out sensitivity. More cells do not fix the eight-donor sample size.

The published paper already reports IL-17 responses. The added analysis should test the **interaction directly**, rather than infer a group difference because one group's DEG list is longer or one within-group P value is significant. Transfer a frozen response program to an independently collected PSC liver dataset only after defining the transfer endpoint; such tissue concordance remains observational. The previously audited liver atlas is not a new cohort merely because a new program is scored in it.

**Causal boundary:** IL-17A exposure is an experimental intervention in cultured organoids. Under the culture design assumptions, it can support an ex vivo response mechanism. The PSC-versus-non-PSC response difference remains partly observational because disease, donor history and culture selection were not randomized. It does not establish that IL-17 causes PSC, that a genetic variant mediates the response, or that IL-17 blockade benefits patients.

Expected compute: **4–8 CPU cores, 8–16 GB RAM**, roughly 15–90 minutes for streaming pseudobulk/QC and donor models; no GPU and no new sequence model. About **145 MB** of merged processed inputs are needed. A dense cell-by-gene matrix is unnecessary. Full matrix acquisition is a separate next action, not completed here.

## 3. PSC liver-state specificity across independent tissue collections

**Pediatric source:** [GSE303271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303271), linked to [the 2026 pediatric PSC paper](https://doi.org/10.1172/jci.insight.199226). The [raw-count matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303271/suppl/GSE303271_raw_counts.txt.gz) is **3,021,812 bytes**; its complete header names 64 sample columns. GEO describes diagnostic liver biopsies and declares GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4.

Crucially, the individual records are **17 PSC, 17 ASC and 30 AIH**, not 34 interchangeable PSC donors. Preserve `ASC` separately until the paper's classification is applied explicitly. Public GEO covariates are tissue and group; 64 unique sample titles do not alone exclude repeat sampling. A prospective PSC-versus-AIH primary analysis can retain ASC as a separate contrast rather than silently merge it into PSC.

**Adult comparison:** [GSE159676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159676), [publication](https://pubmed.ncbi.nlm.nih.gov/35278227/), with a [2,252,130-byte deposited matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159676/matrix/GSE159676_series_matrix.txt.gz). The source uses Affymetrix Human Gene 1.0 ST arrays on explants or diagnostic liver biopsies from the Norwegian biobank. GEO records have **12 PSC, 7 NASH, 3 source-labeled “Primary biliary cirrhosis,” 3 AIH, 1 haemochromatosis, 1 alcohol-related and 6 control-style samples**. The series prose says five controls; preserve the six deposited records and do not invent a donor explanation. PBC nomenclature should be normalized only in a separate reviewed field.

**Design:** fix one epithelial–stromal/fibrosis program and generic immune/fibrosis comparator programs from an external source. Estimate PSC-versus-AIH in pediatric tissue and PSC-versus-disease-control contrasts in adult tissue separately. Evaluate reproducibility at the donor level, retain age/stage/platform differences, and avoid pooling the two matrices. Compare disease tissue against disease tissue rather than only healthy controls. Cell-composition sensitivity can be informative but cannot prove a cell-intrinsic mechanism. With only three adult PBC and three AIH records, results against these individual etiologies are exploratory and imprecise.

**Spatial opportunity, not replication by spots:** [GSE325861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325861) supplies two pediatric PSC diagnostic biopsies profiled by a **custom 400-gene MERSCOPE panel**. A supplementary TAR is publicly listed, but it was not fetched and its complete contents/size are unverified. This is the same study family as GSE303271, has no disease-control group, and is not genome-wide. It could test whether a fixed measurable program localizes near ducts/fibrosis within those two donors. Cells and spatial regions are not independent patient replicates, and co-localization is not cell-to-cell causal signaling.

Expected bulk compute: **2–4 CPU cores, 4–8 GB RAM**, less than an hour for standard small-matrix analysis, excluding annotation qualification. The two compressed expression matrices total **5,273,942 bytes**. Stage, donor mapping and diagnosis covariates—not hardware—are the principal limits. This is a useful bounded pilot, but weaker than the first two paths for a strong PSC-specific mechanism claim.

The suggested [GSE61260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61260) lead was checked. It is a 134-sample liver age/BMI study of morbidly obese patients and controls; the inspected series record does not establish a PSC cohort. It should not be relabeled as independent PSC replication.

## Genetics alternative screened but not promoted

A parallel metadata-only scout found a public [FinnGen R13 manifest](https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv) listing `K11_CHOLANGI_STRICT` with 685 cases and 469,984 controls, alongside IBD-related phenotypes. This is a possible future source, not a qualified PSC-specific design. The exact [R12 endpoint definition](https://r12.risteys.finngen.fi/endpoints/K11_CHOLANGI_STRICT) combines cholangitis codes with rheumatic-disease reimbursement and excludes overlapping disease/reimbursement categories from controls. Release-matched R13 rules, effect convention, indexing, overlap covariance and a suitable PSC-versus-IBD-only comparison remain unverified. Its listed summary file was not fetched or separately HTTP-tested. Marginal overlapping FinnGen endpoints cannot be subtracted as independent case–case estimates.

The [GWAS Catalog exact-trait query](https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait=primary%20sclerosing%20cholangitis&size=100) returned five studies; only the already-used GCST004030 flagged full summary statistics. This bounded search does not prove no other data exist. The old meta-analysis includes Finnish replication recruitment, so full historical-meta-analysis independence from FinnGen also cannot be assumed. These gaps make genetics a lower-priority option than a ready disease-specific transcriptomic analysis today.

## Transfer, reproducibility and status boundaries

The original discovery limit was 100 MB. A subsequent explicit parent instruction authorized **up to 250 MB** for complete blood input qualification. Total recorded response bodies across this scout and the genetics scout were **126,129,047 bytes**, including bounded prefixes and an aborted 10,027,008-byte GEO-family transfer. The aborted adult-liver family file was neither decompressed nor used; the metadata-only `view=quick` route replaced it. No file exceeding 100 MB was fully downloaded. No raw sequencing archive was downloaded.

Local source caches and receipts are under `data/raw/public-data-discovery-options/`:

- `discovery-plan.json`, `source-question-selection.json`, scope/authorization amendments;
- `source-receipts.json`, `blood-input-acquisition.json`, complete input hashes and HTTP ranges;
- `selected-sample-metadata.json`, `matrix-header-join.json`, `blood-input-structure.json`;
- `organoid-donor-and-validation-age-map.json` and source-format summaries;
- `genetics-scout/` with the eight-endpoint independent genetics search;
- native-project acquisition and structure-check scripts. Session cookie headers were omitted from retained request receipts.

A premature agent-message transcription misstated gene counts, annotation fields and the array value range. An immediate explicit correction was sent before effect analysis; **the source files and `blood-input-structure.json` were correct and unchanged**. `structural-message-erratum.json` preserves this correction. The numerical facts in this report follow the checked JSON, not the superseded message.

Selected SHA-256 pins:

| Input | SHA-256 |
|---|---|
| Norwegian PSC count matrix | `d7ba09dab9da6ce78bace86c1e868f2312c1906d4d32d679f8cd7127caf6b6ff` |
| Polish series matrix | `bef0f89eab904937ba63e1cfe89ee33e4478e57d91b18f3d897950e45b62176c` |
| GPL10558 annotation | `c914fdbe1130906ce3b9c97f5a75280591c89c617c168fb687c641f683e76b45` |

**Completed here:** source discovery, public access checks, limited full-input acquisition and structural qualification. **Proposed separately:** expression models, program specificity tests, organoid interactions and spatial analyses. No new biological association, genetic posterior, perturbation result, clinical outcome or therapeutic claim was produced by this scout.

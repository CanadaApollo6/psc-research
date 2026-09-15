# PSC research roadmap

Updated September 15, 2026 after the [two-cohort blood reanalysis](reports/psc-blood-replication.md).

**Next compute priority: a donor-level PSC-by-IL-17A response interaction in cholangiocyte organoids, not a blood classifier or another model sweep.** The completed blood analysis found 93 dependence-sensitive BH conjunction candidates, zero BY candidates and no additional disease-control conjunction support. It does not establish PSC-specific mechanisms. [Independent calculations](reports/psc-blood-independent-verification.json) and a [clean replay](reports/psc-blood-replay-verification.json) agree with the frozen results.

This is the current execution roadmap. The [original four-direction roadmap](docs/computational-followup-roadmap.md), [original protocol and dated backlog](docs/research-plan.md), frozen analysis plans and completed reports remain the historical record. The proposed organoid and reserve calibration analyses have **not** been run. No request to a researcher or laboratory experiment is implied.

## Current priorities

| Priority | Question or branch | Next substantive action | Limits and prerequisites |
|---|---|---|---|
| **1** | PSC cholangiocyte response to IL-17A | Complete GSE239283 matrix/feature/barcode qualification, freeze a donor-paired pseudobulk interaction, then estimate every donor's change. | Four PSC and four procedure-control donors. Age, sampling and clinical differences remain; cells are not independent replicates. |
| Conditional follow-up | Disease-tissue relevance | Transfer a separately frozen response estimand to a qualified independent PSC liver cohort if the organoid result justifies it. | Do not use the same atlas as new replication or select a transfer endpoint after its effects. Small disease-control groups/stage differences remain. |
| **Reserve** | ETS2 program context | Use the now-qualified MacroMap release for fixed program/comparator responses and held-out line/run tests under a new plan. | Healthy donor-derived macrophages test stimulus context, not PSC diagnosis or ETS2 dependence. |
| External-data priority | PRKD2 shared genetic signal | Use qualified original model/covariance or complete-component inputs when available. | Existing expression associations and metadata do not fill the [remaining inputs](data/derived/ibdverse-model-audit/remaining-required-inputs.csv). |
| Bounded/laboratory question | UBASH3A complete RNA/protein product | Retain the specific missing molecular observations and any narrowly justified new library follow-up. | Coverage-limited absence and model scores cannot establish a complete product or its absence. |

## Next compute milestone: donor-paired organoid interaction

The [source comparison](reports/public-data-discovery-options.md#2-donor-paired-psc-cholangiocyte-il-17a-response) verifies GSE239283: four PSC and four non-PSC donor keys, each with vehicle and 100 ng/mL IL-17A for 24 hours. There are 16 libraries and 46,343 cells. The complete cell metadata is acquired; the approximately 142-MB matrix was header-checked, not fully acquired or analyzed.

1. Acquire and hash the complete matrix, feature and barcode files. Check integer counts, exact source library/donor joins and the released feature universe. Do not claim reconstruction of an unverified historical full library.
2. Freeze the count filter, identifier rules, normalization, all-epithelial primary population, donor-level estimand, contrasts and missingness rules **before new effects**. Author-defined treatment-sensitive clusters are secondary, not outcome-selected primary populations.
3. Estimate `mean(IL17 − vehicle in PSC donors) − mean(IL17 − vehicle in non-PSC donors)`. Use donor blocking and an interaction without a redundant disease main effect, or compare donor-level paired changes. Do not infer an interaction from unequal within-group DEG-list lengths.
4. Show every donor-level change, uncertainty, fixed broad inflammatory/stress comparisons and leave-one-donor-out sensitivity. Four-versus-four donors limits power regardless of cell count.
5. If the result warrants tissue follow-up, freeze that new endpoint and source choice separately. A concordant tissue association is not automatically causal validation.

The non-PSC donors are procedure controls, not healthy volunteers. Ex vivo IL-17A exposure is an intervention; PSC status and donor history are not randomized. A positive response interaction would not establish that IL-17 causes PSC or that IL-17 blockade benefits patients. Expected resources are about 145 MB of processed inputs, 4–8 CPU cores and 8–16 GB RAM; no new sequence model or GPU is needed.

## Reserve compute milestone: fixed ETS2 context calibration

The [MacroMap qualification](reports/research-path-review.md) now establishes 58,243 released gene rows, 4,698 samples and 209 source lines, with 205 mapped HipSci identities. Eighteen strict non-PIC stimulation/time comparisons have 169–193 mapped-line control pairs. Preserve the two PIC rows as unavailable for a strict matched-background endpoint because the required mock-control identity is absent. Hold whole runs/lines out, stratify protocols, and learn any scaling/reference matching only in training controls. No program scores or held-out effects have been calculated.

The [GSE255234 qualification](reports/ets2-study-qualification.md) resolves H5 as 5 µM heme for 6 h, alongside H10 and LPS, with four paired donors. It still does not satisfy the unchanged historical full-universe CPM provenance gate. Do not download roughly 73.7 GB of FASTQs merely to bypass that gate when a stronger reserve dataset is available. A separately declared released-universe rank estimand would be a new analysis, not a retrospective pass of the old one.

For the separate GSE84161 inhibitor arrays, compound identities are resolved but anonymous donor pairing and actual array-dose provenance remain distinct questions. An unpaired/arm-mean analysis would need its own question and limits rather than guessed pairs.

## Earlier branch-specific dependencies and completed foundation

The retained sections below describe specific missing study inputs and experimental observations. They do not override the new first compute priority or imply that contacts have occurred.

## Existing information to request from study teams

These requests concern existing files or metadata and could unlock computation without commissioning a new experiment. No outreach has been sent. Anonymous mappings, aggregate counts or appropriately qualified summary inputs should be used where sufficient; personal participant identities are not needed for the initial clarification requests.

| Study team | Specific assistance | What it unlocks | Prepared evidence |
|---|---|---|---|
| **IBDverse** | For classical/intermediate blood monocytes: actual covariate names/counts and model sample records, genotype/phenotype inclusion and QC, selected expression-PC outputs, executed software/container/configuration records, and signed study dosage LD/covariance with allele/order provenance or complete per-signal posterior/BF vectors and conditioning details. | Qualifying the original RNA model and separating its signals for the PSC comparison. The observed RNA-eligible counts 94/93 and conditional PC arithmetic 18/20 are not confirmed executed-run settings. | [Model audit](reports/ibdverse-model-audit.md) and [exact required-input table](data/derived/ibdverse-model-audit/remaining-required-inputs.csv). |
| **PSC genetic-study team** | Per-variant platform sample counts, cross-platform overlap/covariance or adequate study-matched LD with justified sample assumptions, and coefficient/SE scale and allele/QC provenance. | A defensible PSC side of the shared-signal model. Published platform totals alone do not reconstruct covariance. | [Original signal analysis](reports/prkd2-signal-analysis.md) and [remaining-input table](data/derived/ibdverse-model-audit/remaining-required-inputs.csv). |
| **GSE84161 / Senger study team** | Confirm the anonymous sample-to-donor mapping across all six arms and the actual inhibitor doses used for the deposited arrays. G-432/G-573 identities are already resolved. | The planned donor-paired ETS2-program pharmacology analysis of the normalized arrays. | The two current questions in the [supplement follow-up](reports/gse84161-supplement-followup.md). The [older unsent draft](docs/gse84161-source-questions.md) retains historical questions and should be narrowed before use. |
| **Original ETS2 study team** | Complete inhibitor differential-expression tables before symbol collapse and the exact filtering, annotation, duplicate/tie handling and rank-export code/environment. | More complete reproduction of the original inhibitor analysis. This is separate from, and need not block, a newly qualified independent context study. | [Unsent focused draft](docs/ets2-lee-request-draft.md) and [reproducibility brief](docs/ets2-lee-reproducibility-brief.md). |

PRKD2's original IBDverse associations now pass the fixed marginal coverage check, but covariance/model qualification still fails. The intermediate lead belongs to the author's rs313839 clump; clump membership is not a signed LD matrix or proof of shared PSC causality. No new H4 posterior has been calculated. A future reference-LD sensitivity analysis would need its own prospective scope and could examine assumptions without resolving the missing study facts.

## Work that needs suitable experimental measurements

For **UBASH3A**, the immediate laboratory question is whether the complete +29 RNA exists in primary CD4 cells and how its downstream structure affects decay. The [laboratory brief](docs/ubash3a-cd4-experiment-brief.md), [proposed analysis plan](docs/ubash3a-cd4-analysis-plan.md) and [portable packet](reports/ubash3a-cd4-lab-handoff.zip) are prepared. A collaborator would need to qualify the assay and controls, provide appropriately sourced samples, establish complete RNA linkage and measure turnover. Protein production would then require its own informative measurement. The proposed four-donor pilot is for feasibility, not a powered confirmation study; the experiment has not begun.

The additional public-data option is deliberately narrower: seek diagnostic reads in further relevant CD4 libraries and independent datasets under a new fixed selection and artifact-check plan. More libraries from the same TRAILS donor would expand cell-state coverage, not provide independent donor replication. Keep the original boundary hypothesis and the alternative-ending T29 candidate distinct. A new protein search should first qualify the experiment's ability to measure ordinary UBASH3A; shared peptides cannot uniquely identify an originating RNA.

For **PRKD2**, suitable allele/regulatory perturbation and gene-function experiments would eventually be needed to move from genetic association and model predictions to a demonstrated cellular mechanism. For **ETS2**, causal attribution requires a suitable direct perturbation experiment with documented controls and target engagement. An adequate existing public experiment could supply evidence for these questions; our completed searches have not established all of it. These needs are distinct from obtaining missing summary files.

Additional UBASH3A protein-structure modeling is deferred until RNA/protein existence has stronger support. Sequence reconstruction remains a conditional prediction and does not measure RNA survival or protein production.

## Completed work retained as the foundation

| Completed stage | Main result and its boundary |
|---|---|
| [UBASH3A original RNA reads](reports/ubash3a-raw-read-audit.md) | Two libraries, 8,401,078 reads: no exact +29 alignment and 80 qualifying normal-junction reads. Other relevant subsets and donors remain untested. |
| [UBASH3A alternative-ending reconstruction](reports/ubash3a-alternative-ending.md) | A 4,392-nt reference RNA, conditional 482-aa products and four candidate peptides. Neither surviving RNA nor protein is established by this reconstruction. |
| [UBASH3A public protein search](reports/ubash3a-proteomics-analysis.md) | No candidate among 269,986 exported assignments in the fixed fourteen-file experiment; ordinary UBASH3A also lacks a qualifying assignment, limiting the negative result. |
| [PRKD2 initial shared-signal analysis](reports/prkd2-signal-analysis.md) | All 72 component/prior/sample-size settings fail the fixed coverage criterion in the original two-cohort comparison. These results are preserved. |
| [PRKD2 original-data recovery](reports/prkd2-public-data-recovery.md) | Original IBDverse classical/intermediate monocytes recover the missing rs112445263 association and clear marginal coverage. Model/covariance qualification remains separate. |
| [PRKD2 regulatory mechanism comparison](reports/prkd2-regulatory-mechanism.md) | rs313839 has more coherent model/peak support than rs112445263 in the fixed panel. Training-data reuse, limited comparisons, mixed modalities and strong reference LD constrain causal interpretation. |
| [IBDverse original-model audit](reports/ibdverse-model-audit.md) | RNA eligibility counts, declared software and complete gene/conditional-hit summaries recovered; actual executed settings and suitable covariance/components remain missing. Independent checks and 284 repository tests passed at that stage. |
| [ETS2 benchmark](reports/ets2-program-benchmark.md) and [GSE84161 preparation](reports/gse84161-normalization-and-mapping.md) | A measurable, context-dependent program and a reproducibly normalized pharmacology dataset. New independent-study program comparisons and GSE84161 paired treatment analyses remain unperformed. |

PFKFB3/BCL2L11 remain parked after the [coverage follow-up](reports/variant-coverage-repair.md). New work there should start from a source that resolves the recorded gaps. The previous pilot, matched comparisons, allele audits and other findings remain accessible in the [README](README.md) and [evidence log](docs/evidence-log.md).

## Execution order and decision rules

1. Keep the completed blood result, its negative comparisons and all earlier frozen science unchanged. The result constrains a research decision; it is not a clinical recommendation.
2. Run the organoid interaction only after its own source and method freeze. A weak result is valid; do not change populations, programs or thresholds to manufacture significance.
3. Keep MacroMap as an informative reserve rather than repeating generic model scans or four-donor provenance audits indefinitely.
4. Use public data and local compute first. Researcher contact, service enrollment, private-data access and laboratory arrangements require separate authorization.
5. Preserve sample/donor dependence, source values and all unavailable or contradictory results. Retain source URLs, retrieval dates, build/annotation details and transformations. Validate numerics and replay, but do not equate check volume with biological independence.

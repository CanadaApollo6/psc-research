# PSC research roadmap

Updated September 16, 2026 after the [paired organoid IL-17A analysis](reports/psc-organoid-il17.md).

**Next compute priority: fixed-program MacroMap stimulus context, followed by the qualified pediatric liver disease-control panel.** These branches are unconditional on the organoid result. The goal is discriminating evidence, not another classifier/model sweep or changed thresholds until a desired result appears.

The organoid analysis found eight NB-BH interaction candidates, six BY candidates and zero full-family Welch sensitivity discoveries. All eight directions survive paired-change donor omissions, not NB refits. Only eight donors were studied. These are model-dependent ex-vivo associations, not established PSC mechanisms. The [397-check numerical verification](reports/psc-organoid-independent-verification.json) and [byte-identical replay](reports/psc-organoid-replay-verification.json) do not establish biological replication or small-n calibration.

The [blood result](reports/psc-blood-replication.md) remains unchanged: 93 BH conjunction candidates, zero BY and no stronger disease-control conjunction support. Original frozen plans and the [four-direction roadmap](docs/computational-followup-roadmap.md) remain historical records. Current status is in the [active agenda](docs/psc-analysis-agenda.md); no researcher contact or laboratory experiment is implied.

## Current priorities

| Priority | Question or branch | Next substantive action | Limits and prerequisites |
|---|---|---|---|
| **1** | ETS2 program stimulus context | Finish independent review/executable runner, freeze MacroMap mapping, training-only matching/scaling, whole-run holdouts and one incremental-loss endpoint; then execute. | Healthy-derived macrophages, not PSC diagnosis or a direct test of ETS2 dependence. No real scores yet. |
| **2** | Pediatric liver disease-control programs | Reconcile GSE303271 gene identities for externally fixed programs, then freeze PSC–AIH and distinct ASC sensitivities. | 17 PSC / 17 ASC / 30 AIH; whole biopsy, no healthy controls or exact sample-linked clinical/stage adjustment. |
| **Hold** | Adult liver disease controls | Require corrected source expression or separately authorized raw-array reconstruction before modeling GSE159676. | Six PSC patients with two biopsies each, not 12 independent patients. Mixed-scale values and control/PBC/sarcoidosis label conflicts remain. |
| **Park** | NoPSC nuclear/spatial atlas | Reopen only for a qualified original-count/full-feature/barcode/donor export. | Modalities and possibly earlier biobank cohorts overlap; browser labels alone are not a usable count matrix. |
| **Bounded qualification** | BACH2/miR4464 genotype-associated CD4 state | Establish exact E-MTAB-14013 donor/genotype/count/annotation joins before any analysis. | Eight source PSC donors; not isogenic causality, independent replication or RNA-level proof of a translation mechanism. |
| External-data priority | PRKD2 shared genetic signal | Use qualified actual model/covariance or complete-component inputs when available. | Expression associations and metadata do not supply the [remaining inputs](data/derived/ibdverse-model-audit/remaining-required-inputs.csv). |
| Molecular/laboratory question | UBASH3A complete product | Retain the specific missing RNA/protein observations. | Coverage-limited absence and structure predictions cannot establish a complete product or its absence. |

## Completed organoid milestone and remaining limits

The [frozen protocol](docs/psc-organoid-il17-protocol.md) tests a modeled **log2 ratio of IL-17A treatment ratios**, with donor baselines and disease-by-treatment interaction. The equal-donor log-CPM mean-change contrast is a separately prespecified sensitivity, not the same estimand. Unequal within-group DEG-list lengths do not test interaction.

The merged integer export proved SCT-consistent and was not used as raw RNA. All 16 original Cell Ranger RNA libraries were recovered. The new analysis uses all 36,601 original features for the same 46,343 published retained cells. The small historical RNA-margin difference remains documented; no guessed filtering “repairs” it. The source acquisition stayed below 1 GB, without FASTQ. The [reproduction guide](docs/psc-organoid-reproduction.md) fixes all inputs, software, missingness and numerical safeguards.

Culture treatment does not randomize disease, age, procurement, history or culture selection. The eight candidates are not an authorized organoid-derived transfer program. A new transfer question would need an explicit discovery/transfer plan and a dataset capable of testing the corresponding quantity. The externally fixed liver panel proceeds regardless of this result.

## Next milestone: fixed ETS2 context calibration

[MacroMap source qualification](reports/research-path-review.md) is complete. The methods under review retain complete original-source program memberships, 18 strict non-PIC paired contrasts and a fixed nine-class whole-run held-out comparison. The original ETS2 ZIP is absent; the accepted complete pinned membership-export route is distinct from an atlas-mapped subset. Prediction uses 185 both-time lines (114 SmartSeq2 / 71 NEB), training-only matching/scaling and no model tuning. Conditional fixed-prediction run bootstrap does not estimate full refit uncertainty. Source preparation and synthetic review are not effects.

The [GSE255234 qualification](reports/ets2-study-qualification.md) resolves H5 as 5 μM heme for 6 h, alongside H10 and LPS, with four paired donors. Its historical full-universe CPM gate remains unpassed. Do not download 73.7 GB of FASTQs merely to bypass that gate. A new released-universe estimand would require its own plan, not a retrospective pass. GSE84161 anonymous donor pairing and actual array-dose provenance also remain separate unresolved inputs.

## Model capabilities and stopping rules

The [capability audit](reports/model-capabilities-and-psc-use.md) found no authenticated current AlphaGenome inference in the checked setup and confirmed only public AlphaFold DB metadata retrieval. SDK presence, historical runs, hardware, DB revisions and service entitlement are different facts. Use normal secret configuration, never chat/Git credentials. Models need exact data-gated variant/product inputs and measured controls; they do not replace missing PRKD2 covariance or an unobserved UBASH3A product.

The [finite frontier review](reports/psc-analysis-frontier-review.md) sets contrasts, limits and stopping rules. Complete executable analyses. Park a branch for a real data/access/measurement/experimental blocker rather than repeatedly auditing unchanged inputs. Exhausting this defensible current frontier is not a claim that all possible PSC science is exhausted.

## Earlier branch-specific dependencies and completed foundation

The retained sections below describe missing study inputs and experiments. They do not override the current queue or imply that contacts have occurred.

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

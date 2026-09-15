# PSC research roadmap

Updated September 15, 2026, after the [IBDverse original-model audit](reports/ibdverse-model-audit.md).

**Next compute priority: complete one qualified ETS2 program comparison in an independent public macrophage dataset.** PRKD2 now depends on specific study-data inputs. UBASH3A has a narrower public RNA-confirmation opportunity and a prepared laboratory question.

This is the current execution roadmap. The [original four-direction roadmap](docs/computational-followup-roadmap.md), [research protocol and dated backlog](docs/research-plan.md), frozen analysis plans and completed reports remain the historical record. The next analyses below have not started; documenting them does not imply that a data request was sent or a laboratory experiment arranged.

## Current priorities

| Priority | Branch | Useful work with public data and compute | What permits the next substantive result |
|---|---|---|---|
| 1 | **ETS2: program specificity and transfer** | Qualify an independent macrophage study, then compare the fixed ETS2 program with broad inflammation/stress programs using documented donor pairs and matched controls. GSE255234 is the first candidate; public raw-read reprocessing is a possible alternative to an ambiguous deposited matrix. | A complete, documented expression input and sample design, followed by a frozen dataset-specific analysis. This can test program behavior across conditions; it does not establish ETS2 causality. |
| External-data priority | **PRKD2: shared genetic signal** | Run conditional and shared-signal analyses once the original model and covariance/component inputs qualify. The public-source recovery and mechanism audits are complete within their defined scopes. | Actual RNA model/sample records and suitable study LD/covariance or complete signal components, together with the unresolved PSC sample/covariance and coefficient-scale information. |
| 2, bounded follow-up | **UBASH3A: observed RNA identity** | Inspect additional relevant CD4 sequencing libraries and independent datasets for reads linking the exact +29 junction to the complete RNA. Consider another protein search only in a suitable experiment with demonstrated ordinary UBASH3A coverage. | Informative original reads or spectra. The completed negative/low-coverage searches do not establish biological absence; RNA decay and protein production remain separate measurement needs. |
| Later, conditional | **Additional loci and prediction screens** | Revisit a specific gene or model comparison when new data close an identified coverage, context or independence gap. | A concrete new source and a prospective question. Larger prediction sweeps alone do not resolve the existing evidence gaps. |

## Next compute milestone: one qualified ETS2 comparison

The [completed benchmark](reports/ets2-program-benchmark.md) shows that the published ETS2 program can be measured in PSC liver myeloid cells, but its behavior is context dependent. It has not been validated as an ETS2-specific disease-activity score. The next question is whether its response across independent macrophage conditions differs meaningfully from broader inflammatory and stress responses.

1. **Qualify GSE255234.** The [public study record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255234) describes four donors and 16 control/LPS/heme samples. The [existing preparation](config/ets2-validation-gse255234-preparation.json) retains all 12 planned donor-paired comparisons. Resolve the deposited Ensembl-versus-described-RefSeq annotation mismatch, count-universe/filtering provenance, relevant sample QC and matched control backgrounds. Confirm H5 dose/time or retain that planned arm as unavailable under the existing rules; its name alone is not sufficient evidence.
2. **Choose a documented expression input.** Use the deposited matrix only if it meets the existing assay requirements. GEO also lists raw sequencing in SRA. If the matrix provenance remains unresolved, first verify exact run/sample mappings, public accessibility, read layout, transfer size and compute needs. A separately frozen reprocessing plan could produce our own counts with a pinned reference and complete gene universe. This would be a new documented reanalysis, not a reconstruction of unknown historical processing.
3. **Freeze the executable comparison before scoring.** Apply the [effective ETS2 validation rules](docs/ets2-validation-amendment.md): donor-level comparisons, control-only expression matching, complete arm accounting, fixed programs and broad comparators, explicit mapping coverage and deterministic settings. Preserve original and amended protocols; do not loosen a failed criterion after seeing effects.
4. **Run and report every eligible planned comparison.** Deliver a source/sample manifest, mapping and assay QC, program and comparator results, unavailable contrasts, independent numerical checks and a reproducible report. Distinguish generic stimulation context from MEK pharmacology and direct ETS2 perturbation.

**Completion boundary:** one qualified public study analyzed under fixed rules, with all results and limitations reported. If no suitable input can be established within the frozen scope, finish with the precise failed criteria and source requirements. A weak or nonspecific program response is a useful result; an unqualified dataset is not a substitute for validation.

GSE84161 is a separate pharmacology opportunity. Its [30-array normalization and mapping](reports/gse84161-normalization-and-mapping.md) are complete, and the [supplement investigation](reports/gse84161-supplement-followup.md) identifies G-432 and G-573. Anonymous donor pairing and the actual array doses remain unresolved. No program scores or paired treatment contrasts have been calculated for it. Confirmation of those two details would unlock the planned analysis without repeating the completed normalization.

The [wider candidate inventory](reports/ets2-validation-readiness.md) remains available if the first route cannot qualify. Its older compound-identity uncertainty is superseded by the GSE84161 supplement follow-up. A separate direct primary-human ETS2 perturbation remains an evidence need; context calibration or MEK response cannot stand in for it.

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

1. Qualify and complete one ETS2 public-data comparison, or document exactly why the selected input cannot qualify.
2. Pursue the focused PRKD2 and GSE84161 information requests as a separate external track when outreach is authorized. Existing ETS2 reproduction questions can accompany that track without blocking unrelated qualified work.
3. Scope any additional UBASH3A public RNA audit separately, prioritizing diagnostic read evidence and independent donors. Use the prepared packet for a laboratory feasibility discussion when a collaborator is engaged.
4. Reassess priorities after the first ETS2 result or a substantive external reply. Start no additional shared-cause fit, model sweep or expensive experimental program merely to replace an unresolved measurement with more computation.

Each analysis needs a concrete question, eligible inputs, a prospective scope and a report that retains negative, unavailable and conflicting evidence. Completed analysis, proposed work, observational associations, model predictions and experimental results remain distinct. Changes to priorities do not rewrite the original scientific outputs or frozen analysis rules.

# Sustained PSC analysis agenda

Started September 15, 2026 after the completed blood work was committed and pushed as `48c1c2f`. This is an active research queue, not a claim that later analyses have been performed. The aim is to add defensible biological evidence and useful hypotheses. A computational result is not automatically a breakthrough or a treatment recommendation.

## Execution rule

Continue while an accessible, well-defined analysis can add discriminating evidence. Before each inference, fix its sources, biological units, identity rules, estimand, testing families and failure handling. Preserve null, contradictory and unavailable results. Independent calculations and replay are required, but check volume is not biological independence. Commit and push completed milestones.

Do not keep changing datasets, programs, filters or thresholds until a desired gene becomes significant. A new question needs a new prospective plan. Do not repeat unchanged-input audits, broad model-score sweeps or coverage-limited searches without a specific new observation they can supply.

## Ranked queue

| State | Question | Evidence it can add | Next gate |
|---|---|---|---|
| **Completed and pushed** | Which PSC blood associations reproduce across Norway and Poland, with disease controls? | 93 genes meet the fixed BH conjunction, none meets BY or stronger disease-control support. Not an established PSC-specific signature. | [Result](../reports/psc-blood-replication.md); preserve unchanged. |
| **Frozen; inference starting** | Does the IL-17A response differ between PSC and non-PSC cholangiocyte-organoid donors? | A direct donor-level disease-by-treatment interaction, rather than unequal within-group DEG counts. | GSE239283: raw RNA assay provenance, exact 8 donor pairs/16 libraries, released-cell selection and feature universe; separate execution freeze. |
| **Ranked reserve; not executed** | Does a fixed ETS2 program add stimulus information beyond broad inflammation/stress programs? | Independent healthy-line stimulus context and held-out line/run behavior. | MacroMap source/design qualification is complete; freeze program mapping, reference matching, run/line holdout and endpoints. No PSC-causality claim. |
| **Unconditional next-stage question; not executed** | Do externally fixed epithelial/injury programs distinguish PSC from other autoimmune/chronic liver injury? | Disease-tissue relevance and a disease-control comparison regardless of the organoid outcome. | Qualify GSE303271/GSE159676; fix programs, units, stage/covariate limits and study reuse before tissue effects. An organoid-derived endpoint needs its own explicit discovery/transfer plan. |
| **One bounded source check; not executed** | Can the newer NoPSC atlas distinguish within-compartment response from tissue composition and cirrhosis? | A potentially stronger liver comparison with cirrhotic disease controls. | Qualify full measurements and donor labels; spatial and nuclear modalities overlap and are not independent cohorts. |
| **One bounded source check; not executed** | Does the published BACH2/miR4464 genotype association remain within fixed CD4 states? | Separate donor composition from within-state expression in the eight-donor E-MTAB-14013 study. | Exact donor/genotype/count joins; not isogenic editing, protein-level validation or a new independent study. |
| **Capability audit; no new predictions** | Can AlphaGenome or AlphaFold answer a specific variant or protein question not resolved by the data? | Context-specific model evidence, conditional on actual inputs and independent measured benchmarks. | Confirm available API/inference interfaces and permitted use. A local SDK, public structure download or expected entitlement alone is not proof of full inference access. |
| **External-input dependent** | Can PRKD2 shared-cause or UBASH3A complete-product questions be advanced? | Original model/covariance/complete-component evidence or direct RNA/protein observations. | Existing missing inputs remain. Predictions or more uninformative coverage cannot fill them. |

The [completed bounded frontier review](../reports/psc-analysis-frontier-review.md) defines five finite analyses, existing published findings, source independence limits, concrete contrasts and stopping rules. MacroMap and the externally fixed liver panel are not conditional on a positive organoid result. Preparation of the MacroMap methods is active; its real scores remain blocked until a separate freeze.

## Current organoid safeguard

No new gene-wise effects, target/program scores or disease-by-treatment fits have been run. Integer values in a merged single-cell export do **not** prove raw RNA counts. Complete checks now show that all 46,343 cell totals and positive-feature counts match the metadata SCT assay, while 46,270 cell totals differ from the RNA assay. The merged input is SCT-assay-consistent and is not eligible for raw-RNA NB inference; this does not identify the exact historical export command. A raw-count NB model must not consume transformed counts mislabeled as raw. Original per-library processed RNA matrices are now recovered, with exact joins for every published retained cell. The pre-effect protocol accepts the full original raw feature universe after final independent checks, while preserving its small historical RNA-metadata discrepancy rather than forcing a guessed filter. The bounded processed-input budget is 1 GB; no FASTQ reconstruction is authorized for this stage.

## Stop or escalate only for a substantive reason

A branch stops when its available data cannot identify the question, essential provenance remains unresolved after a bounded search, it duplicates an uninformative test, or the next discriminating evidence requires unavailable access, appropriate experiments or resources. Keep useful alternative branches active. “No further defensible accessible analysis on the current frontier” is a practical stopping statement, not a claim that all possible PSC science has been exhausted.

Use public data and existing local resources first. Do not contact researchers, enroll in services, provision paid compute or make clinical recommendations without the relevant separate authorization. Credentials and personal health information must never enter Git or reports. Expected AlphaGenome/AlphaFold access is being verified, not assumed.

# Independent ETS2 validation: protocol before dataset selection

This protocol follows the completed [ETS2 benchmark](../reports/ets2-program-benchmark.md). Its machine-readable [selection rules](../config/ets2-validation-rules.json) record the freeze time and hashes of the source definitions. No new candidate-dataset program scores or donor-level target effects have been inspected. Published study narratives and the earlier atlas and perturbation results are already known. This is a local prospective analysis protocol, not a registry preregistration or a blinded literature review.

## Question and finished result

Does the unchanged source-defined ETS2 program respond to ETS2 perturbation in a separate human macrophage experiment, and how much of its behavior is shared with general inflammation or stress?

The immediate result is a source-backed candidate inventory, an explicit readiness decision, and a fixed dataset-specific input/contrast manifest for the first eligible experiment. If the search finds only generic stimulation, other species or reused data, report that boundary. Such sources can support context calibration but cannot substitute for independent direct-ETS2 validation.

## Evidence tiers

| Tier | Eligible evidence | Maximum interpretation |
|---|---|---|
| A | Direct ETS2 loss or gain in primary human monocytes/macrophages, separate from signature derivation, with explicit biological donors and controls | Independent test of the program's response to a direct ETS2 perturbation in that experimental setting |
| B | MEK-pathway inhibition in primary human myeloid cells with donor/control mapping | Pharmacological response consistency; does not establish mediation through ETS2 |
| C | Human primary or iPSC-derived macrophage stimulation with donor/control mapping | Cell-state transfer and generic-inflammation calibration |
| D | Cell lines, mouse experiments, mixed tissues, or observational disease comparisons | Model-system or tissue context; cannot serve as the primary human perturbation test |
| R | Original source experiment, reanalysis or a dataset used to derive/select the signature | Reproduction or calibration with reuse explicitly labeled |

Missing processed measurements, unresolved replicate identities and uncontrolled assay/condition confounding are separate readiness failures. A relevant paper with inaccessible data stays in the inventory. Failure to find an eligible study in this bounded search does not prove that no such study exists.

## Dataset selection before effects

For a paired primary analysis, require at least three explicitly identified biological donors with both perturbation and appropriate matched control; repeated libraries or repeated cultures from one donor do not create new donors. Sample labels such as `rep1` are not automatically donor identifiers. Verify actual sample records, not only the publication's headline sample count.

Prefer tier A, then B for a separately labeled pharmacological test. Within a tier, require accessible gene-by-sample raw count data and an auditable sample/condition design. Prefer larger numbers of complete donor pairs, then documented control for batch/assay confounding, then a newer publication date; break remaining ties lexically by accession. iPSC differentiation clones require the originating donor and differentiation-batch structure to be recorded.

The first run is limited to one selected study. Enumerate all eligible perturbation-versus-control doses and times in that study before scoring. Treat them as related comparisons, preserve their condition labels, and do not select the strongest response afterward. If only tier C or D is usable, prepare a separately labeled calibration plan rather than relabeling it as direct validation.

Audit independence on three axes: biological donors/data reuse, use of the candidate dataset in signature derivation or cell-state selection, and research-group/publication relationship. A different accession or author list alone does not prove donor independence. An unrelated cohort already used in the source paper for choosing a cell state is distinct biological data with selection reuse; both facts must remain visible.

Only metadata, file listings, identifiers, matrix headers and published narratives may be examined before the dataset-specific freeze. Save exact URLs, retrieval times, byte counts and hashes for downloaded inputs. Record failed retrievals and unresolved metadata instead of silently substituting another source.

## Unchanged programs and mapping

Use all nine programs from the previous frozen plan: the gRNA1-down primary list; gRNA1-up, gRNA2-down and enhancer-deletion-down companions; inflammatory response, interferon-gamma response, oxidative-stress response and apoptotic signaling; and the primary list minus the mapped union of the four broad programs. Exclude ETS2 itself from every program. The primary source list has 927 genes after that exclusion.

Remap original source members to the new dataset. Do not restrict the new experiment to the 906 genes that happened to map in the liver atlas. Match exact stable Ensembl IDs after removing only numeric versions, or unique exact feature symbols using the declared source annotation. Record unmapped, ambiguous, duplicate and excluded entries. Require at least 80% of intended members and ten mapped genes for each source program. The nonoverlap sensitivity retains its parent's coverage gate and requires at least ten remaining genes.

## Quantification and matching

For raw bulk gene counts, verify dimensions, identifiers, nonnegative finite integer values, sample identity and positive all-feature library totals. Combine documented technical sequencing splits of the same biological sample before normalization. If a donor has multiple biological cultures under one condition, average their transformed expression with equal within-donor weight; do not treat them as independent donors or automatically sum different cultures. Use the original matrix's full gene universe for library totals, including zero-valued genes.

Use log2(CPM + 1) and equal-weight means over the fixed mapped genes, as in the prior benchmark. CPM is computed from raw counts; TPM, FPKM, normalized microarray intensity and already-transformed matrices are not interchangeable with raw counts. A different assay or quantity requires its own prespecified calibration procedure before numerical inspection. The current procedure does not infer missing counts from normalized values.

If a selected dataset is single-cell, aggregate raw counts by donor, condition, assay and original myeloid population before normalization, with the existing 20-cell gate and 50-cell sensitivity. Bulk samples use donor/sample QC rather than a cell-count gate. Assays and cell populations remain separate.

Construct the expression-matching reference within each selected dataset/assay/population: average each gene's log expression equally across the predeclared conditions within donor, then equally across complete donors. Sort all genes by reference abundance and stable identifier into 20 equal-count bins. Exclude the eight mapped source programs plus ETS2 from the control pool. Draw 1,000 matched gene sets, without replacement within each set, matching the target count in each bin. Use the same sets for every donor and related contrast in that stratum. Record the deterministic seed and sample-selection hashes in the dataset-specific plan and audit. Insufficient pools make matched results unavailable; do not widen bins.

The score is the program mean minus the mean of its matched-set means. Report the unadjusted mean, matched mean, mapping/detection coverage and baseline-matching error separately. These references match average expression, not biological independence, variance, gene length or co-regulation. Their percentiles are not donor-level P values.

## Fixed outputs and interpretation

For each planned donor-paired contrast and every program, report available/absent pairs, donor-median treated-minus-control score difference, full range, counts of positive/negative/zero differences, and leave-one-donor-out medians when at least three pairs remain. Keep every contrast, including weak, opposing and unavailable results. Do not treat genes, cells, doses or technical replicates as independent donors. No gene-as-replicate significance test or new clinical efficacy inference is planned.

For direct ETS2 loss, the prespecified primary direction is a lower gRNA1-down score; for gain it is higher. Report the opposite-direction companion and all broad programs even if they disagree. A compatible median, stability to donor removal and a distinction from broad responses are separate findings, not a single invented probability of validation. Small-sample directional agreement alone does not establish a universal activity score.

A pharmacological response can be consistent with the ETS2 mechanism while reflecting many other effects of MEK inhibition. A generic-stimulation response can establish contextual sensitivity but cannot establish ETS2 dependence. Even a successful direct perturbation in healthy-donor macrophages does not predict PSC clinical response.

## Audit and continuation

Pin the selected dataset, sample map, conditions, source files, exact scoring parameters and implementation version before the first new numerical run. Test identifier resolution, replicate handling and the planned contrast accounting; independently check at least one donor calculation and replay aggregate outputs. Keep donor-level intermediates ignored and publish aggregate results and provenance in the private research repository.

The focused [researcher request](ets2-lee-request-draft.md) is prepared separately. A response could close original-study provenance gaps; it does not gate the public-data search. New evidence requiring a changed method is recorded as a dated amendment, preserving this protocol and earlier results.

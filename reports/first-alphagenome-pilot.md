# First AlphaGenome pilot: mixed agreement with published PSC mechanisms

Completed 2026-09-11. Four verified variants across three PSC-associated regions were scored with AlphaGenome 0.9.0 using the `ALL_FOLDS` model selection. The protocol, variant selection and comparison rules were fixed before prediction requests. This is a feasibility and consistency check on known examples, not a discovery or an independent accuracy benchmark.

## Main result

The model ranked the published target gene first for **UBASH3A** and **PRKD2** in the prespecified cell contexts. However, the **UBASH3A expression prediction opposed the published allele-direction comparison**. The two ETS2 variants gave different rankings. Gene identification and the direction of a molecular effect must therefore be assessed separately.

| Region and variant | Change scored | Primary cell context | Published target's rank | Median signed expression score | Comparison with published direction |
|---|---|---|---|---|---|
| UBASH3A, rs1893592 | A → C | CD4-positive T cells | 1 of 30 | −0.034576 | Opposite: source implies increased expression for C versus A; model predicts decreased expression. |
| ETS2, rs2836883, GWAS lead | G → A | CD14-positive monocytes | 27 of 37 | −0.000148 | Same sign, but the predicted effect is very small and the target ranks low. |
| ETS2, rs4817988, highest fine-mapping probability | G → A | CD14-positive monocytes | 4 of 37 | −0.002687 | Not tested: the lead variant's risk allele cannot be transferred to this variant without phase evidence. |
| PRKD2, rs313839 | C → G | CD14-positive monocytes | 1 of 50 | +0.020952 | Not tested: the paper's A risk label does not match the C/G reference-database alleles. |

Scores are AlphaGenome's gene-expression log-fold-change scores for **alternate minus reference**, aggregated across the prespecified primary tracks. They are not disease-risk estimates, causal probabilities, measured expression changes, or treatment effects. No effect-size threshold or statistical significance test was defined for this small pilot. These are four variants at three loci, not four independent studies.

## What was checked before inference

- All four rsIDs and their published positions matched the build-37 records. Their build-38 positions and reference bases were verified against Ensembl release 116.
- The same bases were checked directly against AlphaGenome's GRCh38.p13 FASTA using indexed byte-range requests. No whole genome download was needed.
- rs1893592 has A/C/G alleles in the database. The pilot deliberately used the paper's A/C comparison.
- The API's live metadata confirmed CD4 T-cell, regulatory T-cell, memory T-cell and CD14-monocyte tracks. No RNA-seq track explicitly representing IL-4-stimulated macrophages was identified. Generic tracks do not establish coverage of disease-specific cell states.
- Source-direction comparisons were restricted to rs1893592 and rs2836883. The unresolved rows remained unresolved before and after inference.

## Ranking and baseline

Candidate genes had an annotated transcription start site within the fixed 1,048,576-base window. Ranking used the median absolute expression score across primary-context tracks, restricted to genes returned by the model. Missing gene coverage was reported rather than assigned a zero effect.

The nearest-transcription-start-site baseline, restricted to the same evaluable gene set, selected **RNU6-1149P** for rs1893592, **LINC02940** for rs2836883, **RPL23AP12** for rs4817988, and **PRKD2** for rs313839. This baseline includes all annotated gene types. In this example PRKD2's top model rank adds no improvement over that baseline; the UBASH3A rank is more informative, but its opposite predicted direction remains a failure of the expression-direction check.

The model did not score 10, 12, 12 and 11 of the respective annotated candidate genes. Reference annotations and model gene coverage are not identical. Rankings are within these gene sets and are not calibrated confidence values.

## Interpretation and next experiment

**The next useful work is an allele-and-transcript audit of UBASH3A and PRKD2.** For UBASH3A, the source discusses both expression and intron retention. We should trace the underlying experimental allele conventions and transcript definitions, then inspect the predicted splice-site change at rs1893592. A gene-level splicing score alone cannot establish intron-retention direction. For PRKD2, obtain an authoritative association-statistic row to resolve which allele increases PSC susceptibility. Do not repair either discrepancy by choosing the interpretation that agrees with the model.

For ETS2, the ranking difference motivates a complete credible-set comparison with suitable immune-cell contexts and prespecified comparison variants. It does not establish that rs4817988 is causal or that the model has found a new target.

Before claiming useful generalization, add matched background variants, held-out loci or other genuinely independent evidence, and appropriate cell-state coverage. `ALL_FOLDS` was used here; this is not an out-of-training-region test. Published benchmark evidence can overlap with model training sources.

## Reproduce and inspect

The exact inputs, protocol, source hashes and client commit are versioned. Full prediction tables remain in the ignored local run folder `data/predictions/20260911T223028Z/`.

- [Target summary](pilot-target-summary.csv)
- [Every evaluable gene's rank](pilot-gene-rankings.csv)
- [Target scores in all prespecified cell contexts](pilot-target-track-scores.csv)
- [Nonsigned target splicing scores](pilot-target-splice-scores.csv)
- [Run provenance and output hashes](pilot-run-provenance.json)
- [Fixed protocol](../config/alphagenome-pilot.json)
- [Verified variant inputs](../data/derived/benchmark-variants.csv)

The server's exact internal model build was not exposed by this client response; model selection, client commit, metadata hashes and request dates are recorded. A future server update may change predictions. Full output hashes identify this run.

Biological reference: [Goode et al., 2024](https://www.nature.com/articles/s41467-024-53602-w). Model implementation: [pinned official AlphaGenome client](https://github.com/google-deepmind/alphagenome/tree/aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d). Coordinate and reference conventions: [official FAQ](https://github.com/google-deepmind/alphagenome/blob/aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d/docs/source/faqs.md).

This analysis concerns public disease-susceptibility evidence. It does not establish how to treat established PSC or provide an individual clinical recommendation.

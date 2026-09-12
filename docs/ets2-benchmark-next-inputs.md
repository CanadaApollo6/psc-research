# ETS2 mechanism benchmark: available inputs

Checked 2026-09-12. This is a source-readiness assessment. No ETS2 program score, drug ranking, spatial analysis or treatment simulation has been computed in this phase.

**Later follow-up, September 12:** the [program benchmark](../reports/ets2-program-benchmark.md), [source audit](ets2-source-semantics.md) and [inhibitor-method reconstruction](ets2-inhibitor-reconstruction.md) now complete the bounded analysis below. This document preserves the earlier readiness snapshot and its original boundaries.

The [Stankey et al. primary study](https://www.nature.com/articles/s41586-024-07501-1) provides an existing experimental mechanism to reproduce before extending the PSC candidate analysis. Its [versioned public release](https://zenodo.org/records/10707942) contains code, gene sets, ranked gene statistics and four compressed spatial-count datasets. The 34,636,310-byte ZIP matches the publisher MD5 and our pinned SHA-256. See the [source manifest](../config/ets2-benchmark-sources.json) and [member inventory](../data/derived/ets2-benchmark-inventory.csv).

## What is available

| Input in the versioned ZIP | Potential role | Boundary |
|---|---|---|
| `RNA-seq/RNAseq_CRISPR/ko.g1.res_finalSYMBOL_rank.csv` and accompanying R analysis | Published ETS2-disruption signature | A ranked statistic is not a donor count matrix or fold change. Verify its precise sign and construction against the code before scoring. |
| `RNA-seq/RNAseq_overexpression/overexpression_250ng_rank.csv` and `overexpression_500ng_rank.csv` | Opposing perturbation and dose checks | The code describes limma t-statistic ranks. These are two doses in one experiment, not independent cohorts. |
| `RNA-seq/RNAseq_MEKi/PD_0.1_MEK.resSYMBOL.rank.csv` and `PD_0.5_MEK.resSYMBOL.rank.csv` | Measured intervention signatures for a later benchmark | Treatment-versus-control contrasts are specified in the code; confirm that each saved rank file has that orientation. A transcriptional response does not establish PSC clinical benefit. |
| `RNA-seq/RNAseq_biopsies/ets2_genesetsENSG_MASTER.csv` | Source-defined up/down gene sets | Preserve guide-specific and enhancer-specific columns. Select a definition before looking at its performance in the atlas. |
| `Spatial_transcriptomics/S0_NC_count_matrix.zip`, `S0_PSC_count_matrix.zip`, `S1_NC_count_matrix.zip`, `S1_PSC_count_matrix.zip` | Potential PSC spatial benchmark | Nested archives are inventoried only. Sample, donor, field-of-view and panel coverage still require inspection; names do not prove independent donors. |

The article lists the underlying disruption RNA-seq at EGA **EGAD00001011338**, overexpression RNA-seq at **EGAD00001011341**, and MEK-inhibitor RNA-seq at **EGAD00001011337**. Those are controlled-access routes. The corresponding public ZIP folders do not contain the `raw_counts.txt` inputs referenced by their scripts. A different biopsy folder does include a count file; it must not be substituted for the missing macrophage experiments. Public ranked results support a summary-level reproduction, but cannot reproduce donor-level fitting or a held-out-donor test of those experiments.

## Next bounded analysis

1. Resolve the saved statistics, gene identifiers, contrasts, donor structure and guide/dose relationships. Keep any unresolved sign as unavailable.
2. Freeze one source-defined ETS2 program and a simple scoring baseline, with expression-matched control gene sets and stress/cell-death checks. Exclude ETS2 itself from a program score if testing whether a downstream program adds information beyond ETS2 RNA; document that choice before calculating scores.
3. Evaluate the program within the liver atlas's author-defined monocyte and macrophage populations, with donor aggregation and its processing limitations retained. A generic inflammatory program and an ETS2-dependent program need separate interpretation.
4. Reproduce the expected relationship between measured disruption, overexpression and drug-response summaries. Treat this as reproduction of published evidence. It is not a new discovery, an independent prospective benchmark, or a validated simulation of patient response.

The five-gene atlas analysis supplies cell-context and coverage information for this next step. Expression of ETS2 alone does not show that the full experimental program is active. No individual-level access application or researcher request has been sent.

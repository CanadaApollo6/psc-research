# ETS2/MEK-inhibitor reproducibility brief

Prepared September 12, 2026 for possible author review. This brief concerns source provenance and computation, not a clinical claim. It contains no participant-level expression or personal health information.

**Study:** Stankey CT et al., *A disease-associated gene desert directs macrophage inflammation through ETS2*, Nature (2024), [DOI 10.1038/s41586-024-07501-1](https://doi.org/10.1038/s41586-024-07501-1).

**Code release:** [Zenodo 10707942, v1.0](https://zenodo.org/records/10707942), archive prefix `JamesLeeLab-ETS2manuscript_Stankey_CT_et_al_2024-34a9f49/`. Publisher ZIP MD5: `9d6b7ddcf377cb194bc72732d807ae80`; verified SHA-256: `f1da61f2d58b31b1ce1cc418486222a8701698adda76324a10f5f48ecc53239c`.

## Checks that reproduce

| Input and comparison | Observed result |
|---|---|
| Six columns of `ets2_genesetsENSG_MASTER.csv` versus publisher Supplementary Table S1 | Exact member-set agreement at adjusted P < 0.05 and the corresponding logFC sign |
| `ko.g1.res_finalSYMBOL_rank.csv` versus S1 gRNA1 t statistics | 11,658 unambiguous one-to-one symbol matches agree within absolute tolerance 1e-9; no sign disagreement |
| `PD_0.5_MEK.resSYMBOL.rank.csv` plus `macrophage_genesets_fgsea.csv` versus archived `fgsea_PD_0.5_mac.pathways.csv` | All nine raw ES values agree; maximum absolute error approximately 1.63e-13 |
| Archived 500-nM NES versus publisher Figure 5 source-data sheet `5j` | All nine agree; maximum absolute error approximately 2e-16 |

The ES calculation uses the original rank order, absolute-statistic hit weights with exponent one, and equal miss steps. A separate implementation gives the same result. We compared existing NES values; we did not recreate the historical fgsea random-set null, NES normalization or P values.

## Exact missing link

The released MEK script describes nine samples from three donor labels, a CPM > 0.5 filter in all nine samples, edgeR normalization, a donor-adjusted limma/voom model and `MEK_100 - ctrl` / `MEK_500 - ctrl` contrasts. It creates topTable objects and later reads the two released symbol-rank files. We did not find an active matching export or the gene-symbol mapping/collapse step connecting those objects to those filenames.

Each saved rank has 12,095 unique symbols and finite signed values. `topTable(number=12801, ...)` requests up to 12,801 rows; that request alone does not establish the number of filtered genes or explain the final row difference. Code, the publisher's treated-versus-vehicle description and the numerical pathway chain support the intended 500-nM analysis context. Exact per-gene historical export provenance is still unverified.

**Requested clarification:** complete differential-expression tables for both doses before symbol collapse, retaining Ensembl IDs, symbols, logFC, t statistics, raw and adjusted P values, with explicit contrast directions. The production code or explanation should cover expression/NA filtering, annotation/version and one-to-many mappings, duplicate-symbol resolution, tied-statistic order, export precision/serialization, and row counts before and after these steps, together with its code revision and `sessionInfo()` or package versions. These details would distinguish an exact historical export from another plausible export that gives similar pathway results. Relevant original RNA-seq accession: **EGAD00001011337**. This request does not assume its controlled-access count matrix is publicly downloadable.

## Boundaries

Saved statistics cannot uniquely recover donor expression, logFC or its standard error. The reconstruction does not establish a treatment-response model. Our independent liver-atlas follow-up is a separate observational analysis and is not offered as a replication of the original cell-culture intervention.

Detailed methods and complete pathway comparisons are retained in the [reconstruction](ets2-inhibitor-reconstruction.md) and [nine-row comparison table](../data/derived/ets2-inhibitor-reconstruction-pathways.csv).

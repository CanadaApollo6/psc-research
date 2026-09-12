# ETS2 program benchmark: methods and interpretation

The [plan](../config/ets2-program-plan.json) was frozen on 2026-09-12 at 15:05:09 UTC, before inspecting new multigene atlas scores or experimental concordance summaries. Its SHA-256 is `234560f51bc68af2b9f3258af222894e68c9e992c89ba7d2eae396673fb8c9c9`. Source definitions, gene-set sizes, selected publisher-table rows and the previous five-gene atlas results were already known. This is a focused follow-up to that analysis, not a blinded discovery study.

## Question

Does a published ETS2-dependent gene set provide a measurable, informative description of PSC liver myeloid populations? How does it compare with ETS2 RNA alone and broad inflammation or stress programs?

The study that defined the set experimentally disrupted ETS2 in cultured human macrophages. This analysis transfers the resulting gene list into observational liver data. Finding the list expressed, or finding variation among cell populations, does not demonstrate that ETS2 drives the variation in those tissues.

## Sources and programs

The [Stankey et al. study](https://www.nature.com/articles/s41586-024-07501-1) and its [version 1.0 release](https://zenodo.org/records/10707942) supply source-defined gene sets, code and ranked statistics. The [original manifest](../config/ets2-benchmark-sources.json) remains unchanged. A [separate manifest](../config/ets2-program-sources.json) pins the publisher's Supplementary Tables and Figure 5 source-data workbooks, with byte counts, hashes and retrieval provenance.

The primary program is `ETS2g1_DN`: 928 genes in the released list, including ETS2. The source audit establishes that this list exactly corresponds to Supplementary Table S1 genes with adjusted P < 0.05 and negative log fold change after gRNA1 editing. It is not the larger set obtained with a 0.1 threshold. ETS2 itself is removed before program scoring so that the score does not directly include the single-gene comparator.

| Program | Role |
|---|---|
| ETS2 gRNA1 downregulated genes | Primary ETS2-dependent list |
| ETS2 gRNA1 upregulated genes | Companion with the opposite response to editing |
| ETS2 gRNA2 downregulated genes | Related guide experiment; not an independent cohort |
| chr21 enhancer-deletion downregulated genes | Related regulatory perturbation |
| Inflammatory response | Broad biological comparison |
| Response to interferon gamma | Immune-response comparison |
| Response to oxidative stress | Stress comparison |
| Apoptotic signaling pathway | Cell-death signaling comparison |
| Primary list excluding the four comparison programs | Sensitivity to direct gene overlap |

The four biological comparisons use exact columns in the release's `GOBP_genesets.csv`. They are not negative controls known to be unrelated to ETS2. The final program subtracts their union after exact feature mapping; it removes direct overlap but cannot remove biological co-regulation.

## Liver measurements

The two versioned objects and original annotation exclusions are those in the [completed liver-cell analysis](../reports/liver-atlas-cell-context.md). This phase selects PSC donors and four unmodified author labels: `Monocyte`, `Kupffer`, `LAM-like` and `ActMac`. Whole cells and nuclei remain separate. There is no PSC-versus-control or PSC-versus-PBC expression comparison.

Each donor × assay × author population is one observation for summary purposes. Repeated libraries in that stratum are combined. A donor must have at least 20 retained cells in a population; the complete analysis is also reported at 50 cells. An absent author label is unavailable, not a measured zero.

All stored raw counts are checked. Gene counts are divided by each donor/stratum's sum across all raw features, converted to counts per million, then transformed as log2(CPM + 1). The program mean gives equal weight to its mapped genes, and donor summaries give equal weight to eligible donors. Gene sets retain measured zeros rather than selecting only their most expressed members.

Ensembl IDs are matched exactly after removing only a numeric version suffix. Entries supplied as symbols must match one unique raw feature name. Missing, ambiguous and overlapping mappings are recorded. A source program requires at least 80% of its intended genes to map and at least ten mapped genes. This annotation gate differs from observed RNA detection, which is reported separately.

## Expression-matched references

The comparison asks how each program behaves relative to genes with similar average RNA abundance in this dataset. The reference is constructed separately for cells and nuclei from the primary 20-cell donor groups. First, each gene's log expression is averaged across eligible myeloid populations within a donor; then those donor means are averaged with equal donor weight.

All raw features, including zeros, are ordered by that reference abundance and stable identifier, then divided into 20 equally sized bins. The background pool excludes the union of the eight source programs and ETS2. For each program, 1,000 random gene sets match its gene count in every bin, sampling without replacement within a set. Random sets can share genes with each other. The seed and exact sampling order are frozen. The same sampled sets are used across donors and at the 50-cell sensitivity threshold.

If a bin has too few eligible background genes, matched results are unavailable. Bins are not widened to obtain a preferred result. The unadjusted program mean remains available separately. Matching adequacy and baseline abundance differences are reported.

The adjusted score is the program's mean log expression minus the average mean for its 1,000 matched sets. A negative value means lower expression relative to this chosen background; it does not mean ETS2 is inhibited. A positive value is not a calibrated measure of pathway activation. Matching within the atlas also means that the overall score cannot be treated as a comparison with an external healthy baseline.

The fraction of random sets below a program is a conditional calibration percentile. It is not a population P value: genes can be co-regulated, and the random sets are not new donors. Average-expression matching does not match variance, gene length, detection or cell function.

## Comparisons and specificity

Three contrasts were selected before scoring: activated macrophages minus Kupffer cells (primary), LAM-like macrophages minus Kupffer cells, and monocytes minus Kupffer cells. Only identical donors and assays are paired. At least three donor pairs are required; summaries retain the median, full range and positive-pair count. Leave-one-donor-out medians require at least three pairs to remain and are not confidence intervals.

Within each population, the primary score is compared with ETS2 RNA and the four broad biological scores using descriptive Spearman correlations. At least four donors are required, and constant vectors yield an unavailable correlation. Small donor counts make these correlations unstable. They neither establish incremental predictive accuracy nor prove that ETS2 has an effect independent of inflammation.

All programs, populations, thresholds and unavailable comparisons are reported. Related signatures and overlapping cell/nucleus donors are not independent replications. The source tissue represents advanced PSC; the results cannot directly characterize uncomplicated disease or predict an individual's treatment response.

The author-defined populations are themselves transcriptomic annotations. A difference between an activated-macrophage label and another population can partly reflect the expression features used to label those cells. We have not established an annotation procedure independent of the ETS2 program; this is another reason to interpret cell-state contrasts as descriptive context rather than validation of a causal mechanism.

## Experimental summary benchmark

Publisher Supplementary Tables S1 and S2 are the primary source for ETS2 gRNA1, gRNA2, enhancer deletion and the two overexpression doses. The complete source audit keeps discrepancies in saved rank files separate from those final results. A saved vector is oriented only through numerical agreement with a publisher table or an explicit code path producing that saved file; expected biology alone cannot determine its sign.

The frozen comparisons report shared unique gene counts, signed-statistic correlations, log-fold-change correlations when available, and direction agreement. They do not calculate gene-as-donor P values. Effects on the gRNA1-derived program in gRNA1 itself are circular source reproduction. Guide and dose comparisons share a study and may share participants. Public summaries cannot reproduce donor-level differential-expression fitting without the missing count matrices.

MEK-inhibitor comparisons remain directional only where the saved statistics can be oriented from source evidence. None of these calculations simulates a treatment response, selects a drug, or establishes clinical benefit in PSC.

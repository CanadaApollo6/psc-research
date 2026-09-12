# GSE84161: prospective array calibration framework

This framework follows the [study investigation](../reports/gse84161-study-investigation.md) and preserves the earlier [validation rules](../config/ets2-validation-amended-rules.json). It defines an array-specific method before program effects are inspected. It is **not an executable run freeze**: donor correspondence, annotation crosswalks, software versions and final sample/QC records must be pinned in a separate run manifest before scoring.

## Question and evidence limit

Does the source-defined ETS2-dependent program change under the study-reported MEK inhibitor in LPS-stimulated primary human monocytes, and how does that response compare with TPL2 inhibition and broader inflammatory/stress programs?

This is tier B comparative pharmacology. These cells were cultured briefly after blood isolation; they are not the original ETS2-edited inflammatory macrophages or PSC liver cells. A response does not establish ETS2 mediation, direct binding, a clinically useful dose or treatment efficacy. The compound identity and target-engagement evidence are taken only from verified source records.

## Inputs and biological units

Use all 30 original monocyte CEL files from GSE84161, keeping their GSM identifiers and exact source sample labels. The source has six conditions across five numbered replicate groups. Retain those numbers as **candidate blocks**, not independently confirmed donor identifiers. The paper's donor-adjusted model and five-donor design support pairing, but a source-backed block-to-donor map remains a requirement for donor-paired summaries. Do not infer it from clustering, correlations or treatment responses.

Do not append the related GSE84153 neutrophil RNA-seq samples to the monocyte cohort. The superseries repeats the monocyte records and is not an independent replication. Do not count repeated records or technical file copies as new samples.

## Fresh array processing

Use the official Bioconductor `affy` RMA implementation on all 30 arrays together, with full chip probe definitions, background correction and quantile normalization enabled, `bgversion=2`, and no expression-selected probe subset. Pin R, Bioconductor, `affy`, the HG-U133 Plus 2 CDF and all dependencies in the execution record. RMA outputs log2 expression; do not apply CPM, another logarithm or donor-residual subtraction. [Official RMA documentation](https://www.bioconductor.org/packages/release/bioc/manuals/affy/man/affy.pdf).

Pin the exact GPL570 annotation snapshot and its dates. Keep probe sets with exactly one unique, nonblank Entrez Gene identifier; exclude ambiguous and unmapped probes from gene summaries with counts and reasons. For every retained gene, take the median of all its eligible probe-set log2 RMA measurements in each array. The fixed mapping and aggregation give each gene one value regardless of probe count. Do not choose the most variable, largest-effect or most significant probe. Record the number of contributing probes and the limits of combining probes that may measure different transcripts.

This defines a new calibration. The original deposited table selected high-variance probes and has unresolved historical annotation/processing details. It is retained as source provenance; agreement with that table is not assumed and is not the rule used to select our probes.

## Gene programs and annotation

Use the same nine programs and original Ensembl members as the previous benchmark; exclude ETS2 itself. Obtain a versioned Ensembl-to-Entrez crosswalk and pin its bytes before mapping. Require unambiguous reciprocal one-to-one gene correspondence for score membership; preserve missing, retired, many-to-one and one-to-many cases. Do not rescue a failed mapping with an expression-selected alias. Exact program coverage and all losses must be reported.

Retain the 80% mapped-member and ten-gene minimum gates. The primary source set has 927 members after ETS2 exclusion; do not inherit the liver atlas's smaller intersection. The nonoverlap program removes mapped broad-comparator overlap and retains its parent's coverage gate. All genes with valid measurements remain eligible for the reference universe, including low signal; do not filter on treatment effects or copy the paper's differential-expression threshold.

## Fixed contrasts and matching references

The following seven contrasts account for all conditions. The corresponding source sample pairs are finalized only after donor blocks are confirmed. Each paired summary requires at least three complete biological donors.

| Contrast | Role | Reference used for both sides |
|---|---|---|
| MEK inhibitor + LPS minus vehicle + LPS | Primary pharmacological test | LPS vehicle controls |
| TPL2 inhibitor + LPS minus vehicle + LPS | Related upstream inhibitor | LPS vehicle controls |
| MEK inhibitor minus unstimulated vehicle | Basal drug response | Unstimulated vehicle controls |
| TPL2 inhibitor minus unstimulated vehicle | Basal upstream-inhibitor response | Unstimulated vehicle controls |
| Vehicle + LPS minus unstimulated vehicle | Stimulation-context calibration | Unstimulated vehicle controls, applied to both sides |
| TPL2 inhibitor + LPS minus MEK inhibitor + LPS | Comparative pharmacology under stimulation | LPS vehicle controls |
| TPL2 inhibitor minus MEK inhibitor | Comparative pharmacology at baseline | Unstimulated vehicle controls |

The LPS-versus-basal contrast is computed using a single basal reference for both arrays; it is not a subtraction of independently calibrated basal and LPS scores. The LPS arrays therefore have distinct context-specific score instances where needed. Do not interpret differences between absolute scores generated under different references. No interaction contrast or cross-assay absolute-score comparison is planned.

For each of the two control backgrounds, average each gene's log2 RMA across the confirmed control donors with equal donor weight. A shared control enters once. Order genes by reference abundance with lexicographic stable Entrez identifier as the tie-breaker; split into 20 equal-count bins. Exclude the union of all eight original mapped programs plus ETS2 from the matching pool. Draw 1,000 matched sets without replacement within a set, preserving the bin membership counts. Reuse the same draws across all arrays scored against that control background. Do not widen insufficient bins.

Retain the previous seed algorithm, with dataset key `GSE84161` and provisional stratum keys `monocyte_array_basal_vehicle` and `monocyte_array_lps_vehicle`; the exact implemented strings and software versions must be fixed in the run manifest. The score is the mean program log2 RMA minus the mean of the matched-set means. Report raw program mean, matched background, mapping and matching error separately. These controls match expression abundance, not biological specificity or independent donor variation.

## QC and reporting

All CEL identities, header dimensions, byte hashes and sample links must pass before processing. The full RMA matrix must have the expected probe and sample dimensions and finite values. Preserve source QC and exclusions. Use standard array diagnostics to flag problems, but do not invent an exclusion threshold after inspecting an ETS2 response; any new exclusion rule requires a recorded amendment before dependent effects are inspected. Retain all planned contrasts with unavailable statuses where needed.

The execution manifest must name the diagnostic implementations, reference criteria and any exclusion thresholds before program effects are inspected. Those choices are still pending. A diagnostic flag is not an automatic exclusion under an unspecified rule.

For each eligible program and contrast, report donor-median differences, full range, positive/negative/zero counts and leave-one-donor-out medians where three pairs remain. Do not treat genes or probes as biological replicates. Small donor counts and donor/assay uncertainty preclude a treatment-response claim. If the donor map remains unresolved, do not publish these paired endpoints under an assumption disguised as confirmation.

## Before execution

The [preparation manifest](../config/gse84161-array-preparation.json) records the current readiness state. The next execution freeze must include the confirmed sample/donor map, source-backed compound labels, complete raw-array inventory, CDF and annotation hashes, source-to-array gene mapping, exact normalization/matching parameters, QC decisions and runtime versions. Preparing RMA inputs can proceed independently of researcher correspondence; paired biological interpretation remains subject to the donor gate.

# PSC BACH2: prospective follow-up decision review

**16 September 2026 UTC. Design review only. No expression analysis is authorized or executed.** Existing cached publication, qualification and feature/cell-annotation evidence were used. No new HTTP, matrix, gene/donor measurement, ADT value, effect estimate or score was read or computed. The closed ECM/source files and prior reports remain unchanged.

## Recommendation: one narrow RNA check, then close

**One finite donor-level follow-up has modest but real value as a check of a published RNA observation.** Recommend only the single BACH2 endpoint below, subject to a separate root implementation/analysis freeze. Do **not** continue this branch as a test of protein translation, causal genotype effects, or independent validation.

The reason to continue is specific, not merely that count files now exist. The [Poch/Bahn et al. paper](https://doi.org/10.1016/j.xcrm.2024.101620) reports no BACH2 RNA difference in the **same eight-person CITE/scRNA cohort**. The cached Figure S2B caption explicitly says:

> “BACH2 mRNA expression levels in CD4+ (C0, C1) and CD8+ TN (C3) between SNP-carriers and non-carriers.”

A single donor-equal comparison in the released **C0/C1 naïve CD4 compartment** can put an effect estimate, uncertainty and donor-influence assessment around that population-matched observation. The CD4 choice follows the paper's miR4464/BACH2/iTREG focus, before any deposited BACH2 quantity is inspected. It is **not** an all-CD3 endpoint; CD8 C3 is outside this follow-up.

The caption does not identify whether the original display pooled cells, averaged donors, separated C0/C1, or used a particular test. The proposed pooled-population estimator is therefore **not an exact reproduction of Figure S2B's historical computation**. Methods describe donor-covariate logistic regression for cluster-marker detection versus other cells. That does not identify the BACH2 genotype comparison's statistical unit. Neither a valid donor-level null test nor pseudoreplication should be assumed from that paragraph.

This is an original-data robustness/precision check, not a novel-gene search. If the project needs an independent biological validation, causal allele assignment or a translation endpoint rather than this limited check, **PARK**: the released data cannot supply those answers. A wider gene/state analysis is not a substitute.

## What is qualified, and what remains outside reach

The [completed cohort qualification](psc-bach2-cohort-qualification.md) supplies a common **36,601-feature RNA axis**, **55,460 fixed retained cells**, and eight accession-local donors: four source `SNP` and four `noSNP`, all male PSC cases. It establishes a **released RNA-feature source-count-format quantity**, not untouched original UMIs or validity of a count likelihood. All 38 deposited antibody targets are readable; none is BACH2. The paper has separate western-blot and experimental protein evidence, but this release does not provide a direct BACH2-protein feature. Mature miR4464 and translation are not measured by this 5′ polyA gene-expression endpoint.

The source MEX headers declare `cellranger-6.1.1`; Methods/IDF name Cellranger 3.0.2 and GRCh38 release 93. Identical axes and repeated export headers do not reconstruct the executed reference, export command, historical RNA/SCT layer or missing margins. The Methods' 56,209-cell count differs from the 55,460-cell figure/annotation target. The proposed question uses the **released fixed target**, without repairing either history or adding source-only cells.

Figure 3C concerns four later re-clustered naïve activation states. Their memberships are not deposited. Its caption refers to re-clustering C0, while Results use plural “clusters”; this scope ambiguity is preserved. C0 and C1 must not become substitute SC0–SC3 labels. The separate healthy-donor miR4464/mock mRNA result in Figure 5C is not this eight-donor genotype comparison.

## The one prospective candidate

### Fixed membership and independent units

Use E-MTAB-14013 and the already qualified sample/donor/genotype/barcode joins. Starting from the complete fixed 55,460-row annotation, include only these exact `paper_clusters` values:

| Exact source label | Required short label / Seurat value | Fixed retained cells |
|---|---|---:|
| `C0: CD4+ TN RTE` | `C0` / `0` | 9,320 |
| `C1: CD4+ TN mature` | `C1` / `1` | 7,339 |

The union has **16,659 cells**. All eight donors are represented in each label, including at the already reported 50-cell metadata threshold. That is feasibility evidence, not an RNA-detection gate. Do not introduce a new cell or gene abundance filter, recluster, infer missing activation states, equalize state sizes, downsample cells, or replace donors.

The cell key retains sample identity. Remove only the exact matching `sample_name + '_'` prefix from that sample's annotation barcode. Never join globally by bare barcode. All 3,855 source-only cells remain outside the fixed retained set; they are not replacements for the separate historical cell-count discrepancy.

The sole feature is the exact source triple **`ENSG00000112182` / `BACH2` / `Gene Expression`**. An annotation-only lookup found one identical match at full source-feature position 12,314 in each archive. This is source identity, not independent current-reference or allele reconciliation. **BACH2 abundance, detection and sparsity remain uninspected.** No target set, ADT proxy or additional gene is included.

Replication is **eight donors, not 16,659 cells**, with one unpaired four-versus-four source-group contrast. Use the source `SNP`/`noSNP` labels without inventing genotype dosage, risk nucleotides, strand or REF/ALT coding. Do not equate rs56258221 with other project BACH2 anchors.

### Released-quantity estimand

Pool the complete selected C0/C1 cells **within each donor**. For donor *d*, define:

- **B_d:** the sum of released BACH2 RNA-feature units in those cells.
- **R_d:** the sum over **all 36,601 released RNA features in those same cells**, including implicit-zero rows.
- **Y_d = log2(1 + 10⁶ × B_d / R_d).**

The single estimand is **Δ = mean(Y_d in SNP) − mean(Y_d in noSNP)**, with equal donor weights. RNA and ADT units must never share a denominator. Do not reconstruct SCTransform, apply the historical gene/QC filters again, normalize over selected genes only, or certify these units as original UMIs.

This is a donor-average difference in **log2(CPM+1) of the released RNA-feature quantity**. It is not an absolute RNA amount, transcription rate, ordinary raw-RNA fold change, or protein effect. The +1 convention is fixed before values; low-abundance compression is a limitation, not a reason to tune it later.

**Pooling deliberately retains the observed C0/C1 mixture and RNA-content weighting.** Equal weighting of donors does not equalize cells or states. A difference could arise from C0/C1 redistribution, residual activation mixtures inside either label, or changes in other RNA that alter the denominator. This endpoint cannot separate composition from a cell-intrinsic effect. Do not claim that it completes the broader composition-versus-within-state question in the earlier frontier review.

### One complete inference family

The family is **one gene × one fixed pooled population × one source-group contrast**. No separate C0, C1, CD8, ADT, program, trajectory or clinical-response tests are included.

If separately approved, report Δ, a **nominal, model-based two-sided 95% Welch t interval**, and the two-sided working-model P value. Let `v_g = s_g²/4`, using donor sample variances. Use `SE = sqrt(v_SNP + v_noSNP)` and the actual Welch–Satterthwaite degrees of freedom:

```text
ν = (v_SNP + v_noSNP)² / (v_SNP²/3 + v_noSNP²/3)
```

Do not use an infinite-df normal interval. No multiplicity adjustment is needed for this single declared test; that does not create a genome-wide discovery or an independent confirmation. The review is prospective relative to new deposited-value access, **not blinded to the published narrative**.

The working assumptions are independent source donors and approximately Gaussian donor-level transformed quantities within groups. Four per group cannot establish calibration, representative sampling or clinical exchangeability. Do not switch tests after a normality check, substitute a count-likelihood/NB model, bootstrap cells, or interpret genotype-label permutations as a randomized experiment.

A zero B_d gives Y_d=0 when R_d>0: retain it and flag it, with no detection filter. A missing donor, invalid RNA denominator, changed axis/identity or incomplete fixed membership makes the **whole fixed endpoint unavailable**, rather than reducing the four/four cohort. If the estimated SE is zero/nonfinite or df is invalid, preserve any defined point estimate but leave interval/P unavailable; never force a perfectly certain interval or P=0/P=1.

### Fixed diagnostics, not extra tests

1. Check the frozen source/annotation pins, exact eight-donor joins, complete RNA axis and RNA-only denominator.
2. After authorization, retain anonymous donor-level selected-cell counts, BACH2/source-RNA totals, BACH2-positive-cell counts/fractions and Y_d. Positive released counts are not proof of original UMI detection. Do not discard low or zero observations.
3. Show C0/C1 **metadata cell fractions per donor** only as mixture context. Do not add a genotype-composition test, perform adjustment with them, or claim that cell fractions alone explain the pooled RNA effect.
4. Report **all eight leave-one-whole-donor-out Δ point estimates**, their range and sign changes. Do not report leaveout P values/intervals, delete an influential donor, treat overlapping leaveouts as replications, or use sign stability as an extra discovery gate.

These diagnostics must not become a menu for choosing a favorable normalization, subgroup or estimator. The claim reviewer initially suggested an unlogged RNA-fraction contrast. This review instead fixes the **mean log2(CPM+1) difference** above, before value access. Those are different estimands. The provisional unlogged suggestion is superseded, not an alternate primary or sensitivity endpoint. Root must freeze this exact choice before execution.

## What an answer could constrain

| Prospective outcome | Defensible conclusion | Not established |
|---|---|---|
| A clear donor-level shift on this fixed released scale, with consistent influence diagnostics | An unqualified reading of similar **aggregate naïve-CD4 RNA** needs qualification for this estimand. The donor-level association is not explained merely by giving a donor more cells in the group comparison. | Exact failure of the historical figure/test; within-state regulation; rs56258221 causality; a protein or translation effect. |
| Interval includes zero and allows substantial shifts | This cohort is imprecise for this donor-level RNA contrast. | RNA equality, no genotype biology, or support for a translation-only mechanism. |
| A relatively narrow interval around zero | Bounds on differences on this particular released log-CPM scale, conditional on the working model. | Biological equivalence without an independently justified margin, or unchanged RNA/protein in every cell/state. |
| Sparse/degenerate measurement or one-donor dependence | Limited assay support or fragile source-cohort evidence. | Permission to switch states, genes, tests or donors to rescue the result. |

No equivalence margin is supplied by this review. Do not derive one from the observed interval. Even an RNA difference would not refute the paper's separate reporter/transfection evidence; stable RNA would not itself prove miRNA-mediated translation control.

## Clinical and genetic limits remain decisive

Table S2 already shows group differences in IBD labels, treatment exposure and severity measures. Three SNP-group records carry `PSC-assoc. colitis`, one source record lists Vedolizumab, and severity ranges differ. `/` remains an undefined source literal, not proof of no disease or treatment. There is no verified patient-row → archive key; birth years are not sampling ages. No individual clinical values may be assigned by source row order.

Age, ancestry, collection/run/batch, medication timing, clinical selection and their genotype balance remain unresolved. Eight cases cannot support an unrestricted confounder adjustment. Sorting, cryopreservation and sample-dependent historical retention can also affect the observed population. Case-only recruitment and LD mean this is not an isolated allele intervention. The paper itself notes linkage with rs72928038 and possible miR4464 targets other than BACH2. A new donor-level RNA interval does not remove those issues.

The same original eight donors are **not independent validation**. E-MTAB-14103 has no qualified genotype/participant bridge and cannot be promoted to replication. More cells, leaveouts, repeated export headers or agreement between software implementations do not add independent biological samples.

## Freeze and stop rule

**No implementation or real-value access follows automatically from this recommendation.** Root must first approve this exact limited goal and freeze source/cell/gene pins, code and package versions, output schema, estimand, missingness rules and diagnostics. A future implementation can stream the existing cached sources into bounded donor aggregates; no new download, whole-transcriptome search or reclustering is needed.

If approved: **one extraction, one primary contrast, the fixed diagnostics, then close**, regardless of direction, P value, interval width, sparsity or influence. Any failed fixed input/all-eight-donor gate yields unavailable/PARK, not a repair, donor exclusion or alternate population. A weak/null result triggers no extra genes, BACH2/tolerance programs, ADT proxies or new source search. A positive result likewise triggers no automatic expansion or causal/protein claim. If this limited original-cohort RNA check is not worth the project's attention, stop at the completed qualification rather than manufacture a broader analysis.

## Provenance and review validation

- Existing [metadata qualification](psc-bach2-input-qualification.md), [first-archive qualification](psc-bach2-count-qualification.md), [complete-cohort qualification](psc-bach2-cohort-qualification.md), and their three aggregate JSONs were reused. Source integrity statistics were not used to select a gene effect or population result.
- Cached paper: **PMC11293351.1 / PMID38901430**, *Intergenic risk variant rs56258221 skews the fate of naive CD4+ T cells via miR4464-BACH2 interplay in primary sclerosing cholangitis*; DOI above. Only design/methods and qualitative claim identity were assessed, not plotted values or numerical genotype–gene effects.
- Figure S2 caption: cached `supplement-layout.txt`, lines 49–59; SHA-256 `f943d3ead2f74929170821eb6c4ef0e2dbe7c9415cbd1c205b45be1f200b3fbe`. Its already pinned PDF is SHA-256 `5d37298b5a7efa1b0708f25fd2842c6eb1b0b129f4abfe3467cfb466de7d11ba`, converted previously with `pdftotext -layout` 26.08.0. This review read the caption only, not expression-panel values.
- Cell annotation: `scrna_meta.tsv`, SHA-256 `0ba598b3edd30dc7f2ba76f0e7779725e60cce50daa01a0e2e04faebcec7c26c`. The deterministic label filter defines the proposed population without any expression threshold.
- Annotation-only BACH2 lookup used native `.venv` and cached three-column feature files. Each source file has SHA-256 `fb11f5aa4ca7f366ab3ff0ceb22579a6d9bfd386ce66babc7519c280853b3e1e`. No matrix or gene/donor quantity was opened.
- The independent claim reviewer checked the main text and cached caption, distinguished S2B/3C from the separate transfection endpoint, and supported only the same narrow population-matched candidate. It did not execute a method or certify a new biological result.

Source pins, the annotation lookup receipt, claim audit and exact **proposed—not approved** candidate specification are in `work/psc-bach2-followup-review/`. Only this new report and owned work artifacts were written. No production code/importer changed, so no model run, test-suite rerun or catalog rebuild was warranted for this design-only task.

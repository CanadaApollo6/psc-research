# Inputs that would make the next colocalisation test informative

Prepared September 12, 2026 after the [completed regional assessment](../reports/colocalisation.md). This is a follow-up chosen after the results, not an amendment to the completed analysis. No researcher contact has been made.

**Completed follow-up:** the work specified in sections 1 and 2 has now been assessed in the [variant-coverage report](../reports/variant-coverage-repair.md). Two PFKFB3 deletions are exactly recoverable; the four-base deletion remains absent. All four BCL2L11 variants are measured in the separately selected larger OneK1K cohort, but their associations are weak and regional RNA coverage still fails. Exact DICE filter causes remain unknown. The original specification below is retained; the [updated source questions](coverage-repair-source-questions.md) identify what is still needed.

## 1. Match the important PFKFB3 deletions across sources

The next bounded compute task is to establish exact allele identity and coverage for these GRCh38 RNA variants:

| RNA source variant | Published rsID aliases observed in the credible-set file | Why it matters |
|---|---|---|
| `chr10_6161637_TAGAG_T` | rs149556258; rs71390193 | 47.55% of BLUEPRINT component-1 weight; 7.20% of DICE component-1 weight. |
| `chr10_6161575_TAC_T` | rs138248884; rs1491512823; rs745635008 | 27.58% of DICE component-1 weight. |
| `chr10_6170685_AGTATGGAGGTGGGGACAGG_A` | rs1384613200; rs140499372; rs3834900 | 1.74% of BLUEPRINT component-1 weight; 6.91% of DICE component-1 weight. |

The [original credible-set rows](../data/derived/coloc-published-credible-sets.csv.gz) preserve all aliases; the [full-vector audit](../reports/coloc-summary.json) distinguishes component weights from marginal PIPs. Do not collapse these variants to a single base, match only their rsIDs, or liftover just an indel's anchor coordinate.

Work to complete: obtain the local reference sequences on both builds; normalize each full allele against its source reference; verify a reciprocal mapping of the affected sequence; search the complete bounded GWAS and reference-genotype rows; preserve alternative representations and all missing cases. If that produces an adequate common variant universe, freeze a new analysis before recomputing any posterior. Retain every currently reported SNP result and both coverage failures.

Finished when each important deletion has a verified exact match or a documented source gap, and the full signal-weight coverage is recalculated. A failed match is not a null association. Better coverage would still leave the GWAS reference-LD and platform-covariance limits.

## 2. Resolve the BCL2L11 missing-variant mechanism

The [missing-variant table](../reports/coloc-missing-gwas-variants.csv) enumerates the 20 largest omitted GWAS SNP weights for every gene/context. Start with these GRCh38 allele keys in QTD000439, QTD000444 and QTD000454:

| Source rsID | Unordered allele key; not REF/ALT | Eligible GWAS weight absent from each of these DICE comparisons |
|---|---|---:|
| rs13405741 | `chr2:111155479:C:T` | 19.96% |
| rs72837816 | `chr2:111159291:A:C` | 16.14% |
| rs72836352 | `chr2:111137297:C:T` | 10.39% |
| rs144569746 | `chr2:111150990:C:T` | 9.10% |

Inspect source genotype/imputation and allele-specific filtering metadata, and determine whether original or less-filtered aggregate QTL statistics exist for these exact alleles. Do not manufacture missing RNA effects from LD proxies. If they are not available, assess a separately defined CD4 eQTL cohort with adequate variant coverage and sample size before reading target effects. Record cohort/donor overlap explicitly. The existing three nine-carrier DICE results do not substitute for a strong, well-resolved molecular signal.

Finished when the reason for the missing records is sourced, or a fixed alternative dataset supplies sufficient regional evidence. If neither is obtainable, keep BCL2L11 unresolved and record that limit without spending more AlphaGenome calls on the same evidence gap.

## 3. Ask for the exact original statistical inputs if outreach is authorized

The existing [unsent author questions](finemap-provenance-questions.md) now include the platform-count discovery. The useful request is for build- and allele-coded original GWAS association inputs, per-variant sample scope, original signed in-sample LD or covariance, and matching run manifests. The 11,386 log count equals the Omni subset total; ask whether that explains the archived run instead of assuming it does.

For the eQTL Catalogue maintainers, the precise questions are whether the specified BCL2L11 records were excluded by genotype/imputation filters and which sample metadata correctly identify QTD000464/469. Any response is new source evidence, retained alongside the original records. Sending these questions requires separate authorization.

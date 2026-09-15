# IBDverse original model and data audit

Public-source audit, 15 September 2026. No new association, fine-mapping fit or shared-cause posterior.

The audit recovers the **reported model version, RNA-eligible donor counts and original PRKD2 gene/conditional-hit summaries**. These narrow the next request substantially. A new PSC–PRKD2 shared-signal calculation remains unqualified because the exact executed model records and adequate RNA/PSC covariance, or complete original RNA signal components, are still missing.

## What was recovered

| Evidence | Classical blood monocytes | Intermediate blood monocytes |
|---|---:|---:|
| Original context | `dMean__Myeloid_0_blood_all` | `dMean__Myeloid_1_blood_all` |
| Public RNA-eligible individual labels, at least five cells | **94** | **93** |
| Cells in the released analysis object | 45,193 | 49,333 |
| Previously inferred regression residual degrees of freedom | 69 | 66 |
| Expression PCs implied **if** those RNA counts are the model N and the only other covariates are five genotype PCs | **18** | **20** |
| Original selected PRKD2 lead, GRCh38 REF→ALT | 19:46,765,721 T→A | 19:46,700,002 G→A |
| Lead beta per ALT (SE) | +0.4250392 (0.07871816) | +0.95814425 (0.124381125) |
| Lead nominal P in gene-summary export | 8.9217467×10⁻⁷ | 8.9941011×10⁻¹¹ |
| Original gene FDR, `qval` | **0.00897880** | **6.94302×10⁻⁶** |
| PRKD2 rows in complete conditional-hit file | **1**, rank 1 | **1**, rank 1 |
| Author's global LD-clump label for this lead | 19:46,765,721 T→A | **19:46,718,300 C→G (rs313839)** |

The [full source/model table](../data/derived/ibdverse-model-audit/PRKD2-source-model-summary.csv) retains original numeric strings. The four extracted source rows preserve the separate permutation P, beta-approximation P, FDR `qval`, fitted `true_df` and blank `qvals` fields. The new lead-summary P values differ from the preceding nominal exports by about 5–7×10⁻⁷ in relative terms. The intermediate lead SE also differs: 0.124381125 versus 0.12438112, a 5×10⁻⁹ absolute difference. Both original versions remain in the [field comparison](../data/derived/ibdverse-model-audit/lead-source-field-comparison.csv); none is overwritten or harmonized.

**94/93 are directly observed RNA eligibility counts, not verified genotype-model sample counts.** The public analysis object lacks `Genotyping_ID`, and its exact identity with the historical model input has not been established. We joined its sample identifiers to the complete published sample/individual map, with no unmatched cells, and applied the preselected five-cell threshold. The conditional PC arithmetic is `94−2−5−69=18` and `93−2−5−66=20`; these selected PC values were not found in an executed-run record.

One selected independent hit per context means the source procedure retained one hit. It does not prove exactly one biological causal effect. The intermediate lead's membership in the rs313839 clump is useful context for the earlier mechanism results; it does not establish a shared PSC cause. Global clumping combines lead variants across annotations and does not supply per-context signed correlations.

## Model provenance

The [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13441992/) describes donor-level pseudobulk mean log-normalized expression, inverse-normal transformation, at least five cells per individual, at least 30 individuals per annotation and primary gene expression in at least 20% of individuals. It reports TOPMed-imputed genotypes, ALT dosage, a TSS ±1-Mb window, MAF filtering, five genotype PCs and optimization of expression PCs by the number of FDR-significant eGenes.

The [Reporting Summary, page 1](https://pmc-oa-opendata.s3.amazonaws.com/PMC13441992.1/41586_2026_10627_MOESM2_ESM.pdf) explicitly reports **TensorQTL 1.0.9**. Page 3 confirms five genotype PCs plus expression PCs. The pinned [TensorQTL v1.0.9 source](https://github.com/broadinstitute/tensorqtl/blob/8dfd0b634bb7d66e82d190114a353a73d9f101fa/tensorqtl/core.py) defines regression residual df as `N−2−number_of_covariates`, supporting the interpretation of our earlier 69/66 consistency result. Its allele-count implementation thresholds dosages and casts sums to integers, so dividing `ma_count` by twice the allele frequency is not a justified exact sample-count recovery.

The [study configuration](https://github.com/andersonlab/IBDVerse-sc-eQTL-code/blob/de8596d9d1ab3ed65636f7861127c291aaddade2/configs/qtlmap-multi_tissue_all.config) searches 2–100 expression PCs in steps of two, sets five genotype PCs and leaves extra covariate files blank. The separate `nr_expression_pcs=5` setting belongs to the **disabled SAIGE branch**, not the active TensorQTL model. The linked optimization script selects the most eGenes and, for ties, the fewest PCs.

The study README reports QTLight v1.4.0. No exact v1.4.0 release/tag was listed in the inspected public repository; the named `v1.4` branch resolves to an October 2024 snapshot. This is not an authenticated executed-run revision. That snapshot's container recipe installs TensorQTL without pinning a revision, and its mapping script does not contain the later conditional-hit export. The study's dated container URL, code/configuration snapshots and reported software version therefore **do not replace the run manifest or container digest**. The [repository version table](../data/derived/ibdverse-model-audit/code-repository-versions.csv) preserves these distinctions.

The source field `true_df` is 58.12637/57.293793 in the gene-summary files. It is a fitted quantity in the permutation beta-approximation procedure, not the directly counted participant N or the 69/66 nominal regression residual df. We do not substitute it for either.

## Release and metadata accountability

- [S-BSST2922](https://www.ebi.ac.uk/biostudies/studies/S-BSST2922) lists six files. Its complete 393-GB base-eQTL ZIP inventory contains **6,047 members**: chromosome nominal tables, gene-summary tables and independent-hit tables. It contains no `Covariates.tsv`, PC-optimization record, genotype/phenotype mapping file, execution log or LD matrix. The previous complete directory is reused unchanged.
- The four selected small members were read completely with exact HTTP range identity, local/central filename, decompression EOF, size and CRC checks. They contain 16,209 and 16,037 gene-summary rows, and 2,198 and 2,122 independent-hit rows. All four contain exactly one PRKD2 row. Only **3,806,532 compressed/header bytes** were transferred for these members.
- [E-MTAB-16986](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-16986) lists eight files. Complete small metadata tables contain 732 samples and 421 individuals; blood represents 95 distinct individuals with Crohn's disease before the cell-type/model exclusions. The live project site's label says 733 samples; that label is preserved as a discrepancy with the 732-row metadata and Reporting Summary, not silently used as an analysis denominator.
- The object linked as **Analysis anndata** on the [project site](https://www.ibdverse.info/) contains 2,196,874 observations and ten observation columns. Only sample ID, predicted cell label and tissue were decoded. The published sample map resolves all observations; aggregate counts are in the [metadata receipt](ibdverse-model-audit-h5ad-metadata.json). Expression matrices were not decoded or analyzed. No individual identifiers enter the derived tables.
- The released `LD_clumped_eQTLs.txt` has 195,321 membership rows and three columns, exactly matching the decoded GitHub copy. Its **16 PRKD2 rows form four global clump labels**. This is membership, not a signed LD matrix or complete component weight table. The clumping code accesses internal study genotype files; those paths are not public downloads.
- The 83-MB `base_coloc.zip` has three members, `IBD.gz`, `UC.gz` and `CD.gz`. Their bounded nested-gzip headers contain colocalisation summary columns and no N, covariate or covariance fields. Only headers were inspected; no assertion about PRKD2's presence/absence in the full streams is made. These are IBD/CD/UC results and cannot substitute for a PSC comparison.
- The supplementary fine-mapping methods concern RNF14, not PRKD2. Supplementary Tables 2 and 3 contain no PRKD2 row and do not supply the missing model counts. The project's genotype link leads to an [EGA access-controlled dataset page](https://ega-archive.org/datasets/EGAD00001016069); its public page does not establish a downloadable model-matched genotype panel. No access request was submitted.

The [14-file public inventory](../data/derived/ibdverse-model-audit/public-file-inventory.csv), [base-eQTL member-type counts](../data/derived/ibdverse-model-audit/base-eqtl-member-types.csv), [base-coloc inventory](../data/derived/ibdverse-model-audit/base-coloc-inventory.csv) and acquisition ledgers specify the finite search boundary. Public files elsewhere may exist; this audit does not claim a universal absence.

## Browser discrepancy

With PRKD2, Blood, Myeloid, Cell type and displayed q≤0.05 selected, the [public browser](https://data.ibdverse.info/) returns one classical-monocyte lead. Its [CSV snapshot](../data/derived/ibdverse-model-audit/browser-filtered-export-as-published.csv) reports `qval=0.003`. That matches the archived permutation P **0.002997002997** rounded to three significant digits, whereas the archived gene FDR is **0.008978799297**. We use the explicit archived `qval` for gene FDR and retain the browser field as published. The browser's global clump selection also means a missing intermediate row is not evidence that the original intermediate-monocyte association is absent.

Initial EBI file requests returned 403 and the portal fallback returned 404. After a delayed retry the direct public host responded normally; the requested small metadata and four original members were recovered. All failed and successful receipts are retained. The browser's UI download did not yield a usable download result; its visible public download link subsequently returned the 409-byte CSV through the bounded acquisition helper.

## Exact remaining inputs

1. **For each fixed monocyte context:** the actual `Covariates.tsv` sample header and covariate names/counts, final genotype/phenotype inclusion and QC counts, and `optim_pcs.txt` or `optimise_nPCs-FDR0pt05.txt`. These would confirm whether N=94/93 and expression PCs=18/20 were actually used. The source optimizer publishes these files, but they are omitted from the inspected release. Individual identifier disclosure is unnecessary for a first count/column confirmation.
2. **RNA signal separation:** study-matched signed dosage LD/covariance with exact GRCh38 variant, REF/ALT and row order, donor and covariate provenance; or complete original per-signal posterior/BF vectors with conditioning variants and model parameters. A clump list or one lead per signal is insufficient.
3. **Execution provenance:** exact QTLight/TensorQTL revisions, container digest and final run configuration/command, including dosage, missingness and duplicate-ID handling.
4. **PSC covariance and sample provenance:** per-variant platform N and cross-platform overlap/covariance, or adequate study-matched LD with justified sample assumptions, plus documented coefficient/SE scale. The previously documented platform totals—11,386, 3,504 and 14,890—do not reconstruct that covariance.

The [machine-readable missing-input list](../data/derived/ibdverse-model-audit/remaining-required-inputs.csv) and [reproduction guide](../docs/ibdverse-model-audit-reproduction.md) make this a bounded data/provenance request. The earlier missing-allele recovery, marginal coverage and regulatory predictions remain unchanged. **No new H4 posterior is reported.**

## Verification

All **284 repository tests pass**. Independent verification checks 233 distinct source/receipt files, 158 Git blobs against four pinned repository trees, the complete base-eQTL directory, all four selected members and the donor-count aggregation. Offline replay reproduces **16 generated files exactly**, with the additional coloc inventory independently reconstructed. The catalog rebuild is unchanged; all 1,254 preceding files outside the two navigation documents remain byte-identical, and those two documents receive insertions only.

See the [independent verification](ibdverse-model-audit-verification.json), [offline replay](ibdverse-model-audit-replay-verification.json) and [stage record](ibdverse-model-audit-stage-verification.json) for the evidence and runtime. These verify the audit, not a causal biological mechanism.

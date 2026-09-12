# PFKFB3 deletions and missing BCL2L11 variants

Completed September 12, 2026. This follows the [regional colocalisation assessment](colocalisation.md); its frozen results remain unchanged.

**Two PFKFB3 deletions can be recovered, but the strongest RNA signals still lack adequate coverage. All four missing BCL2L11 variants are measured in a larger, separate cohort, where their RNA associations are weak. Neither finding establishes a shared PSC–RNA mechanism.**

| Question | Finding | Research decision |
|---|---|---|
| Were the PFKFB3 deletions lost because of variant formatting? | Two match exact sequences in the original PSC data and public reference genotypes. The four-base deletion is absent from both checked sources, including equivalent repeat representations. | Preserve the two recovered records. Obtain the four-base deletion's actual PSC statistics before another dominant-signal test. |
| Were the BCL2L11 records missing because we selected the wrong allele or gene? | All four exact variants are present for BCL2L11 in BLUEPRINT. Their positions are absent across every gene in all nine selected DICE archives. | This localises the DICE gap to release-level variant availability; the precise filtering reason remains unknown. |
| Does a larger cohort strengthen BCL2L11? | The reannotated OneK1K release contains all four variants in all three selected CD4 groups. All twelve nominal P values are 0.117–0.952. | The tested expression hypothesis gains no support here. Keep other cell states and molecular mechanisms unresolved. |

**PFKFB3: exact sequence identity, including repeat representations**

Every deletion's complete reference allele and 100 bases on each side maps uniquely, contiguously and reciprocally between GRCh38 and GRCh37. The mapped DNA agrees between reference builds. Independent `bcftools norm` checks agree with our normalization for all six build/variant combinations. Full alternate haplotypes, rather than rsIDs alone, establish identity. [Exact identity audit](../data/derived/coverage-repair-exact-identities.json), [pinned reference sources](../config/coverage-repair-sources.json).

| Deleted bases | GRCh38 position and REF → ALT | GRCh37 position | Exact original PSC row | PSC discovery P |
|---|---|---:|---|---:|
| 4 | 6,161,637 `TAGAG → T` | 6,203,600 | Absent | Unavailable |
| 2 | 6,161,575 `TAC → T` | 6,203,538 | rs138248884 | 0.00282325 |
| 19 | 6,170,685 `AGTATGGAGGTGGGGACAGG → A` | 6,212,648 | rs140499372 | 0.0614214 |

The search audits all 7,639 original PFKFB3-region GWAS records for potentially overlapping representations and normalizes the ten indel records within the new reference block. The four-base deletion has equivalent minimal placements at GRCh37 6,203,600–6,203,617; no original GWAS record with the appropriate length change overlaps that tract. This is a documented source gap, not a measured null association. The two recovered deletions also match explicit REF/ALT records in the bounded 1000 Genomes phase-3 query. Their EUR reference-allele frequencies are 0.2296 and 0.2465, close to the source GWAS risk-allele frequencies of 0.22 and 0.23. Original GWAS allele_0/allele_1 columns are not treated as REF/ALT. [Bounded indel audit](../data/derived/coverage-repair-bounded-gwas-indels.json), [reference query provenance](../config/coverage-repair-reference-queries.json), [Ji et al. discovery study](https://www.nature.com/articles/ng.3745).

Current Ensembl responses merge rs149556258 into rs71390193 and rs140499372 into rs3834900. They describe longer, multiallelic repeat alleles, but normalization reproduces the exact target deletions in both builds. The audit preserves the queried identifiers, returned identifiers and complete allele strings. These identifier changes do not supply missing GWAS statistics.

Using the complete published molecular-component Bayes-factor vectors, the two recovered records would provide this additional signal support:

| Published RNA component | Previous coverage | Exact additions | Recoverable coverage | Fixed 90% RNA gate |
|---|---:|---:|---:|---|
| BLUEPRINT QTD000031, component 1 | 50.47% | 1.76 percentage points | **52.23%** | Fails |
| DICE Treg QTD000469, component 1 | 54.06% | 34.49 percentage points | **88.55%** | Fails |
| BLUEPRINT QTD000031, component 2 | 99.9423% | Negligible | 99.9423% | Already passed |

The missing four-base deletion alone carries 47.55% of BLUEPRINT component 1 and 7.20% of DICE component 1. These are fractions of full component evidence, not disease-causality probabilities. The denominator retains every published variant, including other unmeasured variants. This is a support audit: we have not inserted indels into an old GWAS fit or calculated a new colocalisation posterior. The previous secondary-component result remains intact. [Coverage table](../data/derived/coverage-repair-pfkfb3-signal-coverage.csv), [per-deletion component weights](../data/derived/coverage-repair-pfkfb3-deletion-weights.csv).

**BCL2L11: what is missing, and what is known about filtering**

The four exact GRCh38 alleles are rs13405741 `T>C`, rs72837816 `A>C`, rs72836352 `C>T` and rs144569746 `C>T`. All four also pass reciprocal sequence and public-reference checks. rs72837816 is globally multiallelic; its `A>G` alternative is not substituted for the PSC-specific `A>C` record.

Ten successful nominal-archive queries searched every gene and allele at the four positions. All 36 DICE variant/context combinations have **no row at the position**, while BLUEPRINT supplies all four BCL2L11 effects and 23 rows for other genes. Those four BLUEPRINT effects exactly match the earlier gene queries; they are not newly discovered measurements. BLUEPRINT nominal P values are 0.00361–0.0123, with 32–33 minor-allele carriers. [Forty-cell missingness audit](../data/derived/coverage-repair-bcl2l11-missingness.csv), [BLUEPRINT values](../data/derived/coverage-repair-bcl2l11-blueprint-effects.csv), [complete query records](../config/coverage-repair-catalogue-queries.json).

| Source checked | What it establishes | What it cannot establish |
|---|---|---|
| Original DICE paper | Original association testing excluded variants with MAF below 5% or missingness above 5%, after its array/imputation workflow. | Whether any of these four variants failed a particular original filter. |
| Catalogue release notes and genimpute v22.01.1 | Array studies were re-imputed using the GRCh38 30× reference. The Catalogue profile sets R2 >0.4; the workflow also filters post-imputation MAF >1%, with pre-imputation HWE and missingness rules. | The actual per-variant filter history or run overrides for these DICE records. |
| Public DICE QC report | Sample-level genotype/expression and population QC information. | A variant-level exclusion ledger for these four alleles. |
| Original “unfiltered” TFH/TH17/TH2 downloads | Public aggregate files exist; their headers have no individual genotype columns. | Their target rows were not extracted in this follow-up. |

Sources: [Schmiedel et al. original methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289654/), [Catalogue release notes](https://www.ebi.ac.uk/eqtl/Release_notes/), [versioned Catalogue profile](https://github.com/eQTL-Catalogue/genimpute/blob/v22.01.1/conf/eqtl_catalogue.config), [versioned filtering implementation](https://github.com/eQTL-Catalogue/genimpute/blob/v22.01.1/modules/postimpute_QC.nf), [DICE QC report](https://www.ebi.ac.uk/eqtl/static/QC_reports/Schmiedel_2018_QC_report.html), [original downloads](https://dice-database.org/downloads).

The three original DICE archives total 10.17 GB. Header probes show ordinary gzip rather than indexed BGZF, and their `.tbi` requests return 404. We did not stream them in full or infer target absence from a header. The exact original-study measurements remain unexamined. General filtering rules and absence from the Catalogue are insufficient to assign a particular cause. [Access record](../config/coverage-repair-original-dice-access.json).

**A separately selected larger cohort**

The follow-up rule selected a separate CD4 cohort before target inspection, preferring the reannotated OneK1K release. We froze all three relevant CD4 helper/regulatory groups and the same BCL2L11 TSS ±1 Mb interval. The original and reannotated OneK1K releases reuse the same cohort; they were not counted as two replications. OneK1K is a separate named cohort from DICE/BLUEPRINT, while its cell groups share participants. [Cohort plan](../config/coverage-repair-alternative-plan.json), [OneK1K study resource](https://onek1k.org/), [Yazar et al. (2022)](https://doi.org/10.1126/science.abf3041).

| Reannotated CD4 group | Dataset | Source sample count* | Target minor-allele carriers | P-value range across four targets | Eligible PSC evidence covered | Eligible RNA evidence covered |
|---|---|---:|---:|---:|---:|---:|
| Regulatory T cells | QTD000889 | 528 | 82–85 | 0.117–0.144 | 99.9990% | 72.74% |
| Tcm / naive helper T cells | QTD000891 | 997 | 175–177 | 0.792–0.889 | 99.9990% | 65.34% |
| Tem / effector helper T cells | QTD000892 | 653 | 105–106 | 0.715–0.952 | 99.9990% | 70.63% |

\*The r8 metadata and nominal AN/2 fields agree on 528, 997 and 653. The Catalogue study page lists 981 OneK1K donors, so the count of 997 needs clarification. We preserve these source values without claiming 997 verified independent donors. The coverage calculation uses beta/SE with fixed expression SD=1 and does not require an assumed exact cohort N. Imputation R2 is unreported for these rows. [Source metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/65995fe96c8c6da1e750b51bac3b4c2f5ba9260f/data_tables/dataset_metadata_r8.tsv), [complete target effects](../data/derived/coverage-repair-alternative-target-effects.csv).

The 16,401 retained nominal rows represent 5,009, 5,044 and 5,003 distinct variants before exclusions and alias handling. The primary SNP comparison retains 3,790, 3,815 and 3,786 RNA variants; 2,586, 2,594 and 2,578 overlap the eligible GWAS. It excludes 993–1,000 non-SNV variants per group, 223–229 nonshared palindromic SNPs lacking the frozen frequency check, and one low-MAF SNP in QTD000892. Coverage describes this eligible SNP universe; it is not full molecular-signal coverage including indels. Exact exclusions remain available in the [variant audit](../data/derived/coverage-repair-alternative-variant-audit.csv.gz).

All three comparisons fail the fixed 90% RNA-coverage threshold, despite recovering essentially all eligible PSC evidence. None has a published BCL2L11 credible set. No new colocalisation posterior is reported. The twelve target effects are weak, and the strongest eligible regional nominal P values are 0.000395–0.00129; these do not establish a well-resolved shared molecular signal. This does not rule out BCL2L11 in other cell states, splicing, protein regulation or another mechanism. [Regional coverage and missing-variant tables](../data/derived/coverage-repair-alternative-coverage.json).

**Reproducibility and the next useful input**

The r8 files use a new Parquet location, so the initial TSV paths returned 404. The first metadata-only reader also found no gene-name min/max statistics. Both attempts are preserved. A separate access amendment, frozen before target effects, selects row groups using position bounds and then chromosome/gene identifier columns. Just 29.68 MB of nominal-archive byte ranges were read; every range is hashed, with remote identity checks. The three complete small credible-set files were also inspected. [Query provenance](../config/coverage-repair-parquet-queries.json), [reproduction guide](../docs/coverage-repair-reproduction.md).

All 69 Python tests pass. Exact identities, coverage and target tables replay byte for byte, and the three Parquet extractions replay from the ignored range cache. All 235 earlier colocalisation artifacts retain their recorded hashes. No AlphaGenome request, researcher contact or clinical inference was made.

The most useful next inputs are the **actual PSC statistics for `chr10:6161637 TAGAG>T`, with study allele definitions and matching LD/covariance**, and **the DICE variant-level filtering ledger for the four BCL2L11 alleles**. OneK1K's donor-count discrepancy also deserves a source explanation. [Precise unsent requests](../docs/coverage-repair-source-questions.md) are ready. Another shared-signal analysis should be frozen only when its common variant universe and molecular evidence are adequate; the current results do not justify promoting either gene as a treatment target.

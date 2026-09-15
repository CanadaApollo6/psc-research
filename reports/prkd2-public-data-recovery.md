# PRKD2 public-data recovery

September 15, 2026. This is a follow-up to the [completed signal analysis](prkd2-signal-analysis.md). It searches for missing data under a separately frozen scope; it does not replace the earlier statistical results.

## Decision

The original **IBDverse** data recover **rs112445263–PRKD2**, the association missing from the reprocessed monocyte files. Classical monocytes show P=4.748×10⁻⁵ and intermediate monocytes P=1.153×10⁻⁹, with the working PSC risk allele C associated with lower RNA. Both contexts now pass the unchanged **marginal evidence-coverage** requirement: about 90.70% of PSC evidence and 97.11%/99.9919% of RNA evidence are represented. These are useful candidates for a further signal analysis, not a demonstrated shared PSC cause.

The original DICE download also recovers the association, with P=0.224, but lacks an explicit standard error. No complete original-source component vectors or adequate study-matched LD/covariance have been recovered for a new multi-signal comparison. Exact source sample/model information also remains necessary. **No new shared-cause posterior is reported.**

The original BLUEPRINT file does not contain that association, and the harmonised PSC copy does not recover the four targeted RNA indels. All **33 accessible Catalogue contexts** lack the key PSC variant in their PRKD2 data. Nine other metadata-listed contexts returned unavailable-file responses; their original author sources are assessed separately below. These are correlated contexts and reannotations, not 42 independent cohorts.

All four preselected original IBDverse blood-monocyte members were read completely. Non-classical monocytes fail RNA coverage; the S100A8/9 group's apparent numerical coverage is undermined by an extreme finite-sample distribution limitation described below. The four contexts share one study family.

## What was fixed before the new association queries

The [recovery plan](../config/prkd2-recovery-plan.json) was frozen at 16:52:01 UTC, SHA-256 `8f212c4d24e87c62648f6802edafc36288f833c35eafd9d42efc9136b3aa708f`. It retains PRKD2, ENSG00000105287, and the original GRCh38 TSS ±1-Mb interval: one-based positions **45,717,127 through 47,717,127**. Targets were all previously missing associations contributing at least 1% of a full marginal evidence vector, plus the existing rs313839 identity/direction control. Selection uses the preceding analysis, not effects in the newly inspected sources.

The six [verified targets](../data/derived/prkd2-recovery-targets.csv) are:

| Role | GRCh38 REF→ALT | GRCh37 REF→ALT |
|---|---|---|
| Missing PSC SNP, rs112445263 | 19:46,700,099 C→A | 19:47,203,356 C→A |
| Missing BLUEPRINT RNA insertion | 19:46,681,839 C→CT | 19:47,185,096 C→CT |
| Missing BLUEPRINT RNA insertion | 19:46,685,258 C→CA | 19:47,188,515 C→CA |
| Missing BLUEPRINT RNA deletion | 19:46,681,255 CT→C | 19:47,184,512 CT→C |
| Missing DICE RNA insertion | 19:46,717,081 T→TCGG | 19:47,220,338 T→TCGG |
| Identity/direction control, rs313839 | 19:46,718,300 C→G | 19:47,221,557 C→G |

The 42 eligible monocyte gene-expression contexts come from the pinned Catalogue metadata: 28 release-7 contexts, including the two already examined, and 14 release-8 contexts. Blood, stimulus and secondary lung contexts retain their labels; macrophages and monocyte-derived dendritic cells are excluded. The live metadata repository still matched the cached commit `65995fe96c8c6da1e750b51bac3b4c2f5ba9260f`. The [full context table](../data/derived/prkd2-recovery-cohort-screen.csv) and [252-row target table](../data/derived/prkd2-recovery-target-screen.csv) include failures and multiple probe measurements.

The unchanged model requirements include exact allele identities, at least 100 shared variants, at least 90% of the full evidence in **each** trait, and suitable complete signal vectors or LD/sample-covariance inputs. A source-specific coefficient cannot fill a hole in another pipeline's component vector.

## Original BLUEPRINT and DICE

The [DICE download page](https://dice-database.org/downloads) links the unfiltered classical-monocyte VCF. The [BLUEPRINT archive](https://ftp.ebi.ac.uk/pub/databases/blueprint/blueprint_Epivar/qtl_as/QTL_RESULTS/) provides the complete nominal monocyte summary file. Both complete archives were downloaded and read through decompression EOF. A separate gzip/line-filter/field-parser verification independently reproduces the retained rows; see the [archive verification](prkd2-recovery-archive-verification.json).

| Source | Complete compressed archive | Total source rows | PRKD2 rows retained | rs112445263–PRKD2 |
|---|---:|---:|---:|---|
| Original DICE | 3,289,307,806 bytes | 215,708,754 | 3,060 | Present: P=0.224, beta=+0.22 per A |
| Original BLUEPRINT | 1,016,268,157 bytes | 59,468,244 | 3,099 | Absent after complete full-allele comparison |

DICE has 89 rows for rs112445263 across all genes, one of which is PRKD2. The original DICE export also contains the targeted T→TCGG insertion, with P=0.2129 and beta=+0.30. The three targeted BLUEPRINT indels are absent from both original PRKD2 vectors. These findings concern these specific released analyses, not whether the participants carried the variants.

Original DICE reports rs313839 at P=0.1005, beta=+0.32 per G; original BLUEPRINT reports P=3.078×10⁻²⁵, beta=+1.133, SE=0.1092 per G. The reprocessed DICE value previously reproduced was P=2.289×10⁻¹⁰. The source-specific differences in significance are substantial and remain unresolved. The control's direction agrees with the working assignment of risk C lowering RNA, but this does not make the estimates interchangeable or independently replicate the same analysis.

The original [DICE methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6289654/) and [BLUEPRINT methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5119954/) use their own expression processing and association models. DICE's rounded VCF `Statistic` is not silently divided into beta to invent an SE. BLUEPRINT's monocyte RNA analysis uses 194 donors in the original publication; the reprocessed Catalogue context lists 191. Source sample descriptions do not establish an exact per-variant covariance matrix.

The [original-source overlap audit](prkd2-recovery-analysis.json) preserves all source gene associations in its denominators, including unshared, unmapped and out-of-window records. BLUEPRINT has one exactly duplicated source row; it remains in the retained raw extract and counts once in evidence calculations.

| Source | Distinct shared full edits | Full PSC marginal evidence represented | Original RNA marginal evidence represented | New multi-signal comparison |
|---|---:|---:|---:|---|
| Original DICE | 2,501 | 99.8996% | Unavailable: no explicit complete SE vector | Unqualified |
| Original BLUEPRINT | 2,283 | 61.6979% | 99.99998%, diagnostic calculation | Unqualified; PSC coverage fails |

These percentages are evidence-weight coverage, not percentages of causal probability or validation accuracy. PSC uses the preceding full-vector marginal BF convention, including platform-specific N and case fraction. The BLUEPRINT RNA diagnostic uses the published beta/SE and a quantitative-effect prior SD of 0.15 in source units. It is not a colocalization posterior. No original-source full signal components or matched covariance were recovered.

## Catalogue screen and possible reference-panel explanation

All 26 newly queried release-7 contexts completed. Five release-8 contexts had accessible Parquet files and completed regional PRKD2 extraction; nine returned HTTP 404. Combined with the two existing complete regions, **33 contexts were measured and all lack rs112445263–PRKD2**. For release 7, the exact SNP position was queried across all genes; the primary SNP is absent from those point streams entirely. This bounds the represented PSC marginal weight below about 79.64%, even before other missing variants are considered.

For secondary indels, point absence does not exclude a shifted equivalent representation elsewhere. That limitation does not affect the single-base primary SNP. The target table preserves source-gene/probe distinctions and does not average distinct measurements.

The separately frozen [reference audit](../data/derived/prkd2-recovery-reference-panel-audit.csv) checks the six targets in the public 1000 Genomes high-coverage **20201028_3202_phased** chromosome-19 file. It contains the four exact indels and rs313839 but no site at rs112445263's coordinate. The Catalogue's inspected genotype-processing configuration references this reference-panel family. This is **consistent with a reference-panel coverage gap**, but the exact historical internal reference checksum and variant-level QC records are unavailable. It does not prove the Catalogue filtering cause or invalidate the PSC association. No individual reference genotypes were retained for this new audit.

## Original author sources for the unavailable Catalogue contexts

The [author-source plan](../config/prkd2-recovery-author-plan.json) limits follow-up to the same nine contexts and freezes IBDverse cell-type mapping before inspecting its associations. An HTTP 404 from the Catalogue is not treated as absence from the authors' data.

**IBDverse:** the [project site](https://www.ibdverse.info/) releases a 393,071,674,087-byte ZIP containing separate chromosome/context files. Its cached central directory identifies 6,047 members. The pinned author annotation maps Myeloid_0 to classical monocytes, Myeloid_1 to intermediate monocytes, Myeloid_3 to non-classical monocytes and Myeloid_7 to S100A8/9 monocytes. Only their four chromosome-19 members are selected; dendritic cells and the aggregate myeloid group are excluded. The [publication](https://www.nature.com/articles/s41586-026-10627-z) specifies GRCh38, TOPMed imputation, ALT dosage and nominal TensorQTL associations. The study includes Crohn disease and healthy participants; these are not PSC case–control expression measurements.

The selected members total **561,794,284 compressed bytes**, **17,390,643 source rows** and **14,649 retained PRKD2 associations**. Each complete member passes CRC, size and decompression EOF checks. A separate standard-ZIP reader and field parser reproduce every retained source row; see the [independent verification](prkd2-recovery-ibdverse-verification.json). The [target table](../data/derived/prkd2-recovery-ibdverse-targets.csv) and [full-evidence audit](prkd2-recovery-ibdverse-analysis.json) preserve all four outcomes.

| Original blood-monocyte context | PRKD2 rows | rs112445263 beta per A (SE) | Nominal P | PSC marginal evidence represented | RNA marginal evidence represented | Interpretation |
|---|---:|---:|---:|---:|---:|---|
| Classical | 3,634 | +0.580192 (0.133639) | 4.748×10⁻⁵ | 90.7035% | 97.1135% | Clears marginal coverage; further model inputs required |
| Intermediate | 3,637 | +0.869322 (0.122734) | 1.153×10⁻⁹ | 90.7035% | 99.9919% | Clears marginal coverage; further model inputs required |
| Non-classical | 3,634 | +0.437078 (0.157520) | 0.007011 | 90.7035% | 74.8380% | Fails RNA coverage |
| S100A8/9 | 3,744 | +0.401249 (0.373790) | 0.395396 | 90.7052% | 100%, numerical only | Severe distribution limitation; not qualified |

These are nominal source P values, not gene-level FDRs. Betas are in the original inverse-normal expression units. The marginal RNA coverage calculation uses the same Gaussian approximate-BF convention and 0.15 prior SD as the other new original-source diagnostics; its model assumptions need assessment rather than automatic acceptance.

The S100A8/9 export contains beta/SE ratios as large as about 47 alongside P≈0.000445. The separately [frozen format audit](../config/prkd2-recovery-ibdverse-format-plan.json) therefore compares every PRKD2 row in all four contexts against Student-t distributions with integer residual degrees of freedom 1–500. The best matching values are **69, 66, 73 and 2**, respectively. Across all rows, the largest absolute log₁₀-P discrepancy at the matching distribution is below 7.6×10⁻⁷. Independent R and Python calculations agree. This is an inference about the export's distribution convention, **not a recovered participant count or covariate list**.

With a distribution consistent with just two residual degrees of freedom, S100A8/9's large beta/SE ratios cannot be interpreted as ordinary normal Z scores. Its apparent 100% Gaussian-BF overlap is not evidence of a well-determined shared signal. The [format audit](prkd2-recovery-ibdverse-format-audit.json) retains the original P values, effects, standard errors and numerical coverage; it adds a limitation rather than changing the source or selecting a favorable result. Classical and intermediate monocytes remain the useful candidates, with finite-sample/model assumptions still explicit.

**Nassiri 2025:** the [author repository](https://github.com/isarnassiri/Supplementary_data_Genetic_determinants_of_monocyte_splicing/tree/694ee67e9829a0056c7290e0e8c5f107d7941da1) exposes three gene-level files labelled non-significant. All three were downloaded: 69,271 untreated, 25,822 IFN-γ and 53,862 LPS rows. Their schema contains forward/backward conditional P values, ranks and selected-hit flags, with no explicit allele pair. They contain 46, 48 and zero PRKD2 rows, respectively, and none for the primary SNP. These are conditional exports, not a complete nominal association vector; their omissions cannot establish null associations. The [94 retained PRKD2 rows](../data/derived/prkd2-recovery-nassiri-conditional-rows.csv) remain separate from the marginal source statistics.

**Natri 2024:** the [GSE227136 source directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227136/suppl/) includes public LIMIX beta/SE results and SNP identity information. The relevant immune TAR.GZ is **11,026,774,251 bytes**, with no supplied random-access index, exceeding the frozen 300 MB per-context transfer budget. The separate compressed SNP-information table is about 974 MB. This stage verifies source availability and conventions but does not scan the two secondary lung-monocyte contexts. Their target status remains **unmeasured**, not absent. The README labels the assessed allele as the major/A1 allele; it must not be assumed to be the genomic reference. The separate mashR `lfsr` output is not a nominal P value.

## Harmonised PSC copy

The entire [GCST004030 harmonised archive](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004030/harmonised/) contains 7,866,424 rows and passes the publisher's MD5 `4a9995f393bdf1dd4d09ca9129800b1a`. Both build controls agree with the verified GRCh38 alleles; no coordinate-metadata contradiction was found. The `base_pair_location` column in this harmonised release is also already converted, so it must not be treated as the original GRCh37 coordinate.

The six [target outcomes](../data/derived/prkd2-recovery-harmonised-targets.csv) recover the two known SNPs and **none of the four RNA indels**. The entire union of both build intervals was retained before identity matching. There are 4,668 fully verified rows with harmonised coordinates in the fixed GRCh38 window, compared with 4,727 in the earlier original-source evidence envelope. The 64 retained union rows without a harmonised variant ID are also examined in both source-allele orientations; none supplies a target match, and ambiguity/reference failures remain explicit. Identifier changes and source-envelope differences are preserved in the [secondary audit](prkd2-recovery-secondary-audit.json); this copy is not used to shrink the original denominator.

For rs112445263, the source C-effect OR 1.318 becomes A-effect OR 0.7587253414. For rs313839, C-effect OR 1.322 becomes G-effect OR 0.7564296520. Both reciprocal checks pass. These are allele representations of the same PSC study, not new biological replication; the published label conflicts remain recorded.

## Verification and remaining requirements

The [reproduction guide](../docs/prkd2-recovery-reproduction.md) records source URLs, complete-archive and byte-range hashes, extraction/replay commands, the corrected auxiliary control query and the boundaries of portable verification. Original source rows and all earlier PRKD2 analyses are preserved. This stage emits no new H4 posterior and does not treat unavailable inputs as zero evidence for a shared cause.

All **266 repository tests pass**. Independent R coverage calculations agree with Python within 3.7×10⁻¹⁵, and the distribution-convention checks agree within 2.3×10⁻¹³. All **27 calculated outputs replay exactly** in an isolated work copy. The [final verification record](prkd2-recovery-stage-verification.json) confirms preservation of 993 prior tracked files, with four intentional navigation/evidence updates.

The immediate useful follow-up is now specific: obtain **complete original IBDverse classical/intermediate PRKD2 signal vectors or suitable study-matched LD/covariance**, together with the exact donor counts, genotype/expression covariates and variant-QC provenance used in those exports. Those contexts already contain the missing SNP and clear the marginal coverage gate. The original PSC platform-specific N/covariance limitations also remain. Catalogue reanalysis sample sizes are not silently assigned to the original exports. The inspected ZIP contains nominal, gene-summary and independent-hit tables; an independent-hit table is not a full component vector.

The original DICE row further narrows the source-provenance question but cannot repair a different reprocessed model. A laboratory mechanism test, personal-genome interpretation and clinical intervention remain outside these public-data results.

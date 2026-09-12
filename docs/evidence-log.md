# Evidence and decisions

## 2026-09-11 — Scope and current results

The primary independent project is regulatory variant interpretation. A liver-atlas metadata audit is available as a secondary resource. The first four-variant AlphaGenome run is now complete; its results are below and in the [pilot report](../reports/first-alphagenome-pilot.md). No liver-atlas expression-matrix analyses, new gene discoveries, treatment effects, or personal risk estimates have been produced.

The current import contains a signal summary table, molecular colocalisation rows and sequencing-library metadata. The machine-readable counts are in [data-audit.json](../reports/data-audit.json).

## 2026-09-11 — Source inconsistencies to preserve

Source: [Goode et al. (2024)](https://www.nature.com/articles/s41467-024-53602-w), Table 1, Supplementary Data 2, and the corresponding Results text.

- The Results text describes credible sets ranging from 1 to 62 variants, but Table 1 reports **122** for the second MST1 signal. The import preserves 122; its correctness is unresolved.
- Table 1 spells a candidate-gene label **BCL211**. This is retained verbatim and is not treated as a validated gene identifier.
- The AP003774.1/CCDC88B discussion describes rs663743 as a single-variant resolution in one passage but also describes two signals with two-member sets. Table 1 contains the latter pattern. Do not silently resolve this conflict or use its prose as an exact variant manifest.
- Some figure captions, coordinates and probabilities differ from corresponding prose/table entries. In particular, Table 1 gives UBASH3A's highest PSC probability as 0.62, while the Results text rounds it to 61%. Use the explicit source field and do not equate PSC and molecular-QTL probabilities.
- Supplementary Data 2 contains a tissue-state typo, `Unstimulatedimulatedimulated`. The imported field is marked `as_published`; a future normalized field can record a reviewed correction.

These observations concern published source consistency. They do not establish that the authors' underlying analyses are incorrect. The pipeline does not infer complete credible sets from set sizes or the highest-posterior SNP.

## 2026-09-11 — GEO metadata audit

Public metadata for [GSE243977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243977) and [GSE247128](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247128) contain 45 library records: 26 control-style titles, 15 PSC and 4 PBC.

Title prefixes yield 21 control, 10 PSC and 3 PBC donor labels. These are provisional labels; repeated libraries are not independent donors. The [paper](https://pubmed.ncbi.nlm.nih.gov/38199298/) describes 24 control livers, so these two subseries alone do not establish the full donor mapping for every published analysis.

The automated check flags 19 records: ten title/description chemistry discrepancies and nine nucleus/cell wording discrepancies. Generic descriptions may explain these. Preserve the raw labels, inspect fuller public annotations, and seek clarification if an analysis depends on unresolved mappings. This does not block the separate genetics project.

## 2026-09-11 — Research interpretation

- Benchmark genes are existing findings. Reproducing them is a method check.
- AlphaGenome predicts molecular effects; its scores are not PSC causal probabilities.
- Model training overlap and correlated evidence can inflate apparent validation.
- Cell-state coverage matters. Absence of a predicted effect in an unsuitable track is not evidence of no biological effect.
- Genetic susceptibility, established-disease progression, drug efficacy and clinical safety are different claims.
- Public professional contacts in the directory have not agreed to collaborate; no outreach was sent during repository setup.

## 2026-09-11 — Verified variants and completed model pilot

Twenty-seven pinned reference responses document build-37/build-38 variants, Ensembl release 116, nearby genes and AlphaGenome's own reference bases. Four selected substitutions passed coordinate and reference checks. rs1893592 is multiallelic; C was selected because the published comparison is A/C. rs313839 has C/G alleles, so the A risk label in Supplementary Data 2 cannot support a directional comparison. rs4817988 has no directly tabulated risk allele in the selected source, and the lead rs2836883's direction was not transferred.

The fixed `ALL_FOLDS` feasibility run ranked UBASH3A 1/30 but predicted the opposite expression direction; ETS2 27/37 for rs2836883 and 4/37 for rs4817988; and PRKD2 1/50 with source direction excluded. Reported ranks use the same evaluable gene universe as the nearest-TSS baseline. This is not a held-out or matched-control validation. Missing immune stimulation states and model/source overlap limit interpretation.

The next scientific priority is to trace allele conventions and transcript-level evidence for the disagreements, not to reinterpret the source to favor the model.

## 2026-09-11 — Completed allele/transcript follow-up

See the [full audit](../reports/mechanism-audit.md). The earlier model predictions and prespecified comparisons remain unchanged.

- UBASH3A expression directions disagree in the primary literature. Newman 2017 reports elevated exon/junction coverage and intron retention; Ge and Concannon 2018 report reduced total mRNA with no significant absolute change in their assayed intron-9 transcript, increasing its ratio to total RNA. Todd 2018 explicitly discussed the conflicting results. Different endpoints and contexts must be kept separate.
- Ensembl transcript annotation places rs1893592 at +3 of the same intron in five transcripts. Exon/intron numbering differs between transcript models; the physical donor boundary is shared.
- Our additional prediction's v0.1 endpoint had a **one-base indexing error**. The original config is archived. Official model code labels the last exonic base for positive-strand donor tracks, whereas the client transcript helper reports the first intronic base. The corrected index and erroneous original reading are both published; no prediction was repeated. This is an exploratory, corrected analysis.
- At that boundary A→C reduces predicted donor usage in both CD4 assay tracks and increases predicted usage of a donor 29 nt downstream. The original maximum donor score corresponds to that downstream gain. The 29-nt splice form was measured previously by Mucaki et al. 2020, but their one-individual-per-genotype lymphoblastoid comparison does not agree with the model's direction. Existence of a splice form is not validation of an allele effect.
- The exact original GWAS row labels **rs313839-C** as PSC risk (OR 1.322; control frequency 0.84). Forward-strand C frequency is 0.837 in 1000 Genomes EUR, supporting C as the working risk allele for this palindromic C/G SNP. Goode's thesis explicitly labels G as risk, and the 2024 supplement labels A. Preserve both conflicts; direct molecular-QTL coefficient alignment remains outstanding.
- The original PRKD2-region lead is rs60652743-A. Carryover of that A into the supplement is a possible explanation, not a demonstrated correction. rs313829, shown in the regional figure, maps to chromosome 7 and cannot substitute for rs313839 on chromosome 19.
- The targeted original GWAS rows are discovery-stage results, not the full meta-analysis or the later fine-mapping dataset. Their P-values and sample scope must not be silently substituted.

No outreach, personal-genome analysis, intervention recommendation or clinical efficacy claim resulted from this work.

## 2026-09-11 — Measured RNA and external cohort checks

See the [full QTL follow-up](../reports/qtl-followup.md). Eight indexed archive queries preserve 121 association rows. Thirty-nine public source files pin metadata, headers, indexes and supporting documents. No further AlphaGenome request was made.

- Direct PRKD2 rs313839 C>G coefficients are positive in reprocessed BLUEPRINT monocytes (beta 0.908576, P 1.7097e-23, n=191) and DICE classical monocytes (0.921304, P 2.28923e-10, n=91). ALT is explicitly the effect allele. With working risk C, this supports risk-associated lower RNA and agrees with the model. BLUEPRINT reuses the original cohort; DICE is separate. This does not retroactively turn the frozen pilot exclusion into a prespecified success.
- UBASH3A rs1893592 A>C coefficients are positive in DICE resting CD4 (1.24972, P 1.21399e-23, n=88) and 4-hour CD3/CD28-activated CD4 (0.710505, P 3.66393e-06, n=89), opposite the model's negative gene-level score. The two DICE conditions share donors.
- No UBASH3A splice trait is assigned in the retrieved DICE phenotype metadata. Its filtered splice indexes omit chromosome 21. Neither observation is a measured null effect. Catalogue connected-component files retain selected tag traits rather than every phenotype.
- Lepik_2017 blood contains the exact +29 boundary after documented coordinate conversion, with beta 0.523483 per C (P 1.98464e-15, n=471). The source cluster is labeled minus while the assigned gene is plus. A genome-wide metadata check finds 115,561/116,015 unique single-gene phenotypes disagree in this dataset (99.61%). Preserve both labels; the result supports boundary/direction agreement provisionally, with unresolved strand provenance and mixed-cell context.
- TwinsUK has a positive association for a different exon-spanning junction, not the +29 junction. Its sample count is not evidence of independent unrelated donors. Both the positive association and the missing +29 measurement remain reported.
- Supplementary Table 2 of AlphaGenome supplies 18 relevant track metadata rows, with ENCODE accessions. Donor-level overlap with QTL cohorts remains unresolved, and ALL_FOLDS does not hold out these genomic regions.

The original literature labels, initial wrong donor index, corrected model results and opposing prior splice measurements remain in the record. Researcher questions were updated but not sent.

## 2026-09-11 — Strand mechanism and direct public-junction check

See the [new report](../reports/ubash3a-junction-followup.md). Source retrieval continued into September 12 UTC. No new AlphaGenome request or researcher contact occurred.

- Ensembl reference DNA and both independent Snaptron junction records support positive-strand GT–AG boundaries for the canonical and +29-nt forms. Lepik's minus suffix is preserved. The published workflow maps its reverse-stranded branch to RF in HISAT2 and to `-s 2` in regtools, while the declared regtools 0.6.0 recipe defines RF as `-s 1`. A source-derived paired-read truth table demonstrates how this would reverse labels; actual historical parameters/container contents remain unverified. This is a plausible explanation, not a confirmed correction.
- Before target count/genotype-group inspection, a frozen plan specified the two junctions, donor aggregation, European primary cohort, 10-read threshold and 20-read sensitivity check. GEUVADIS was already used by Ji et al.; this analysis is exploratory reuse, not a new cohort or held-out model validation.
- The two complete Snaptron rows matched the indexed archive. Direct ID queries confirmed all 667 sample/run mappings against recount3 and ENA. They represent 464 donor aliases, with 447 matched public 1000 Genomes genotypes. The target VCF record has no rsID; its GRCh37 position and A/C alleles are explicitly validated and the missing ID preserved.
- Of 360 genotyped European donors, only 18 have at least 10 combined junction reads: 12 AA, 6 AC and 0 CC. The 20-read threshold retains five. The intended population-adjusted HC3 interval is undefined with singleton population groups; both analyses remain not estimable. Sparse coverage is not a null biological effect. Nine +29 supporting reads were observed across eight of all 464 RNA donors; this cannot establish a stable allele effect.
- The Lepik cluster lead file records seven tested traits but selects the distinct exon-spanning endpoint. Its coverage archive contains a UBASH3A base-coverage PDF, not a donor-by-junction matrix. Public maintainer correspondence from 2023 explains the absence of full splicing files and offers a targeted-transfer route in that case. A short unsent email and six-row endpoint specification now request the exact blood/CD4 comparison and historical strand settings.
- Thirty-five source files and the separate single-variant genotype query are pinned. Full public sample-level joins stay ignored; aggregate counts, missingness and model status are versioned. The next independent compute task is complete credible-set inspection and a controlled extension, with the UBASH3A measurement gap kept explicit.

## 2026-09-12 — Complete fine-mapping archive and fixed comparison rules

See the [archive audit and comparison report](../reports/psc-finemap-and-comparison.md). The 16,716,205-byte fine-mapping archive from [Zenodo 11143561](https://zenodo.org/records/11143561) matches the publisher MD5 and our pinned SHA-256. All 71,083 variant rows are imported: 46,248 GWAS rows across 18 regions and 24,835 molecular-QTL rows across eight datasets. This counts dataset rows, not independent variants or replications. The separate 7.64 GB expression-QTL archive is not downloaded.

- All 25 nonempty configuration tables reproduce their matching SNP-table marginal probabilities within four-decimal rounding precision. Preserve both source PIPs and separately calculated configuration marginals. `-inf` source log Bayes factors remain `-inf`; a missing configuration is not zero probability.
- All 26 logs identify FINEMAP v1.1, while Methods state v1.3. Sixteen GWAS configurations contain multivariant models despite a one-causal-variant maximum in their logs. The posterior count-of-causal-variants distributions disagree in 23/25 comparable datasets, including every molecular-QTL dataset with a configuration table. One PRKD2-region monocyte configuration is empty; `_2CV` is a separate source, not a replacement.
- No explicit credible-set, original association-input or LD files are in the release. We cannot reproduce exact per-signal published sets or certify run provenance. Marginal PIPs must not be renormalized into a credible set for multicausal regions. GWAS logs' 11,386-individual field must not be silently equated with the paper figures' full sample scope.
- UBASH3A and SIK2 alone have complete singleton-only configurations. Conditional reconstructions yield 95% sets of 7 and 790 variants, respectively, and 99% sets of 103 and 1,803. These are separately labeled reconstructions, not claims to recover the paper's sets. Table 1's five-member UBASH3A result remains unchanged; its exact construction rule is unresolved.
- The Table 1 crosswalk preserves literal rsIDs, semicolon aliases, provisional coordinate matches and unresolved identifiers. CD28 signal 2's rs231799 has 0.17 in the table versus 0.0004 in the archive. BCL211/BCL2L11, IL2-IL21/IL21 and ETS2/PSMG1 are region links rather than validated equivalences between gene names. GPR35, HDAC7 and SIK2 occur in the archive outside the Table 1 region list.
- Before new model scores, a deterministic rule selected IL2RA, BACH2 and BCL2L11. All 12 records with archived PIP >=0.1 are accounted for. Nine literal rsIDs have two-build database mappings; six are reference-verified biallelic SNVs. Five become candidate anchors and rs7750559 remains a reserve. rs4147359 and rs72837816 are multiallelic, rs68021656 is an indel, and three anonymous/array identifiers remain unresolved. Do not replace these with guessed SNP alleles.
- The fixed protocol requires two low-PIP comparison variants per matched anchor with EUR frequency, mutation/context, TSS-distance and LD checks, an identical per-region window/gene universe, fixed CD4 RNA-seq tracks and region-level descriptive contrasts. Comparator identities and matching covariates remain unprepared; `ready_for_inference` is false. Incomplete LD coverage of ambiguous high-PIP records must be disclosed. ALL_FOLDS and uncertain source probabilities do not provide independent accuracy labels.
- Forty-three public sources are pinned; no AlphaGenome prediction or researcher contact occurred. Unsent provenance questions are available. A reply would improve confidence but is not required to prepare the limited prospective comparison. Existing model predictions and clinical interpretation remain unchanged.

## 2026-09-12 — Completed matched comparison, including unmatched candidates

See the [complete report](../reports/matched-comparison.md). The original protocol and preparation freeze remain unchanged. The final nine-input manifest was locked at 04:33:28 UTC; the nine AlphaGenome requests ran successfully from 04:34:00 to 04:34:08 UTC. This is a local prospective freeze, not a registry preregistration or held-out model benchmark.

- Every literal rsID with low source PIP in the selected regions was queried on both builds: 6,117 records, with 535 meeting the initial annotation/distance criteria. A further 1,257 non-rsID source records did not meet the fixed identifier requirement. Canonical replacements returned by Ensembl were preserved in raw responses but not silently substituted; 92 original rsIDs lacked a usable direct mapping. These are documented mapping exclusions, not incomplete network requests.
- All 541 targets (535 potential comparisons plus six eligible high-PIP SNVs) match both Ensembl references and the model reference. Exact public 1000 Genomes phase 3 genotypes exist for 540; BACH2 rs9362750 G>A is absent. Frequencies and dosage LD use 503 EUR donors. The six excluded high-PIP identifiers/alleles remain outside the LD screen, as declared before inference.
- Three anchors obtained fixed comparison pairs: rs7923054 with rs7919234/rs111845808; rs72837826 with rs56366063/rs35718992; rs13405741 with rs4848397/rs10176092. All assigned comparisons use 503 jointly called donors and satisfy the fixed frequency, context, distance and LD gates. Each pair's mutual r² is below 0.015. No comparison variant is reused.
- Both BACH2 anchors failed matching. Of 126 initial candidates, 124 fail the frequency requirement for each anchor and one lacks an exact genotype. The remaining candidate qualifies for rs7750271 only, leaving one eligible comparator there and zero for rs72928086. Neither can form a pair. Criteria were not relaxed, and neither anchor received a prediction.
- The final gene universes were fixed before scoring from Ensembl 116 and the model's GENCODE v46 annotation: 30 IL2RA genes, 16 BACH2 genes and 23 BCL2L11 genes. The two merged CD4 RNA tracks were pinned and live metadata checked before prediction. All nine requested variants have all fixed genes and both tracks: 228 variant–gene pairs and 456 signed scores. No missing output was converted to zero or used to shrink a universe.
- Primary anchor contrasts are −0.0007176101 (IL2RA rs7923054), +0.0008407832 (BCL2L11 rs72837826) and −0.0003911853 (BCL2L11 rs13405741). The BCL2L11 regional median is +0.0002247990; the equal-weight median across two evaluable regions is −0.0002464056. BACH2 is unavailable. There is no consistent model-score separation in this small comparison; correlated anchors are not independent replications.
- Maximum-scoring genes are noncoding RNAs or processed pseudogenes; none matches its region's name. IL2RA itself ranks second for rs7923054. BCL2L11 ranks 17th for rs72837826 with mixed signs across its two tracks, and 14th for rs13405741. The raw scorer is natural-log ALT-minus-REF with a pseudocount; no risk direction or percent change in measured expression is inferred. Underlying expression means are not retained, so expression-level/pseudocount sensitivity remains untested.
- Eighty-five matching sources, a separate metadata snapshot and three indexed genotype subsets are pinned. Full donor-level inputs, credentials and full predictions remain ignored. All 42 tests pass. No outreach or new clinical finding resulted. The next independent compute phase will account for measured RNA evidence across all nine variants and all 228 fixed gene pairs; it is a follow-up chosen after this model run, not a revised primary endpoint.

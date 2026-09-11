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

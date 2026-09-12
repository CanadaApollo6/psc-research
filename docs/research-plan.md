# First computational project

Updated 2026-09-12. Status: the pilot, allele/transcript audit, measured-RNA follow-up, direct public-junction analysis and [complete fine-mapping archive audit](../reports/psc-finemap-and-comparison.md) are complete. The archive contains 18 GWAS regions and eight molecular-QTL datasets, with unresolved log/configuration provenance. A [prospective comparison protocol](../config/psc-controlled-comparison.json) and five reference-verified candidate substitutions are fixed. Comparator matching is next; no predictions have been requested for the new regions. The earlier UBASH3A public-junction analysis remains inconclusive because of sparse coverage.

## Question and finished result

Can regulatory-sequence predictions help prioritize **variant → gene → cell type → molecular effect** hypotheses at PSC-associated regions?

The first finished result is a reproducible feasibility pilot on **UBASH3A, ETS2, and PRKD2**, including selected variants, genes, tracks, parameters and disagreements. A follow-up resolves transcript numbering, documents a corrected donor-track index and defines a +29-nt UBASH3A splice hypothesis. Negative or inconclusive results are retained.

A prospective comparison at IL2RA, BACH2 and BCL2L11 is now defined before viewing their prediction scores. It will compare model effect magnitudes for statistically prioritized variants with matched low-PIP variants. This is not an independent accuracy benchmark: archive probabilities are uncertain labels, causal genes are not established at these regions, and the fixed ALL_FOLDS model does not hold out these genomic regions. A larger computer run alone is not evidence of a better biological result.

## Work we can do independently

1. **Prepare genetic inputs.** Start with the published Table 1 summaries and Supplementary Data 2. Resolve each selected rsID against a versioned variant database, verify build 37 versus build 38, check the reference sequence, and record both alleles and the risk-allele direction. Preserve all alternative mappings and exclude ambiguous cases with a reason. A lead variant and the highest-posterior variant can differ.
2. **Define the benchmark.** Treat the existing gene links as prior evidence, not discoveries. Enumerate candidate genes using a fixed rule, inspect the model's available tissues and immune-cell states, and declare expected outputs before querying. Do not silently substitute whole liver for a missing immune-cell or stimulation context.
3. **Run a bounded AlphaGenome pilot.** Use precomputed Atlas scores where appropriate; use hosted reference-versus-alternate predictions for additional molecular detail. Record client/model versions, tracks, score definitions, sequence windows, and missing coverage. Keep risk-versus-protective effects distinct from alternate-versus-reference effects. No local GPU is required for hosted inference.
4. **Compare against evidence and simple baselines.** Evaluate gene ranking against a nearest-gene baseline, check direction when allele alignment permits it, and include prespecified matched comparison variants. Related variants from one locus are not independent tests. Report all three loci, not only successful examples. Investigate sensitivity to justified track and window choices without selecting the best-looking run as the main result.
5. **Assess what is new.** Check for prior experiments and other variant-prioritization studies. Track overlap between model training data and benchmark evidence; rediscovering training-associated signals is not independent validation. Use a genuinely separate dataset or experiment where available, and mark independence as unknown where it cannot be established.
6. **Extend only with interpretable evidence.** Obtain complete credible-set memberships for the next loci if available; the imported summary table does not provide them. Keep fine-mapping probabilities, model scores, and external evidence in separate columns. Do not multiply them into a claimed causal probability without a justified statistical model.

## Starting evidence

- [Goode et al. (2024)](https://www.nature.com/articles/s41467-024-53602-w): PSC fine-mapping and molecular-QTL evidence. Some underlying individual-level data are controlled-access. Published tables and the linked summary-statistic releases can be used independently.
- [Study-generated fine-mapping and QTL release](https://zenodo.org/records/11143561): the 16.7 MB fine-mapping archive is downloaded, checked and imported. The separate 7.64 GB expression-QTL archive remains undownloaded.
- [Stimulated macrophage data](https://zenodo.org/records/7967759): candidate source for cell-state context; not downloaded here.
- [PSC liver atlas](https://pubmed.ncbi.nlm.nih.gov/38199298/): potential cell-expression context. Current work covers metadata only. Expression differences in advanced-disease tissue would not establish a causal variant effect or directly describe early disease.
- [AlphaGenome official client](https://github.com/google-deepmind/alphagenome) and [Atlas announcement](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/): model access and scope.

## Backlog and acceptance criteria

| Order | Work item | Finished when |
|---|---|---|
| 1 | Normalize the three regions' inputs | Completed for four selected variants. Coordinates and bases verified; two source-direction comparisons remain unresolved and were excluded. |
| 2 | Freeze pilot and inspect track coverage | Completed for the feasibility run. A matched comparison-variant set remains necessary before a controlled benchmark. |
| 3 | Connect AlphaGenome and run the pilot | Completed: four successful requests with pinned client and run provenance. |
| 4 | Produce pilot report | Completed for all three regions, including the UBASH3A direction disagreement and ETS2 variant-dependent ranking. |
| 5 | Audit source alleles and transcripts | Completed: source-specific RNA endpoints, shared transcript boundary, one additional splice request, original GWAS rows and allele-frequency checks. C is the working PRKD2 risk allele; contradictory source labels and molecular-QTL alignment remain explicit. |
| 6 | Align PRKD2 molecular-QTL coefficients | Completed for reprocessed BLUEPRINT and separate DICE monocytes: G increases expression in both, agreeing with risk C lowering expression and with the model. Exact original-paper labels remain an author clarification; historical exclusions remain unchanged. |
| 7 | Test UBASH3A splice hypothesis independently | Public summaries and direct GEUVADIS reanalysis completed. Only 18/360 genotyped European donors meet the fixed junction-read threshold; no CC donor is retained. Record the comparison as inconclusive. Motifs establish compatible positive-strand boundaries; historical workflow settings still need confirmation. A six-endpoint blood/CD4 data request and short email are drafted, not sent. This is not independent model validation. |
| 8 | Audit the complete fine-mapping archive and define the comparison | Completed: all 71,083 variant rows imported; all 25 nonempty SNP/config pairs agree within rounding, but 23 corresponding posterior-count distributions disagree with logs. Original per-signal credible sets cannot be validated from this release; singleton reconstructions remain separate. Five candidate substitutions are verified on both builds and the model reference, with every high-PIP exclusion recorded. Selection and matching/analysis rules are fixed before scores. |
| 9 | Prepare matched comparison inputs | Next independent compute phase: verify EUR MAF and pairwise LD, annotate the low-PIP pool, assign two comparison variants per eligible anchor using the fixed protocol, and pin one common gene universe, exact tracks and interval per region. Preserve the current freeze record; store new matching sources in a separate manifest. Record unmatched anchors and incomplete high-PIP LD-screen coverage. Finished when the exact matching table and final input manifest are locked before inference. |
| 10 | Run and report the fixed prospective comparison | Run only complete matched groups: up to 15 variant requests for the five current anchors, within the protocol's 18-request cap. Report every comparison and exclusion with regions as the reporting unit. Retain signed molecular scores separately from magnitude contrasts; do not claim held-out accuracy, a causal posterior or a treatment effect. |

## What this could contribute

A useful result could identify a previously underexplored regulatory mechanism, improve the ranking of candidates for an experiment, or establish where a model fails on PSC-relevant biology. Publication or novelty is not guaranteed. A mechanism affecting disease susceptibility may not control established disease or make a safe drug target. Those questions require further evidence, including experiments and clinical studies.

Researcher outreach can help review the result and propose experiments. It does not hold up the public-data work above. No personal genome or patient records are required for this project.

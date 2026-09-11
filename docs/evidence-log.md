# Evidence and decisions

## 2026-09-11 — Scope and current results

The primary independent project is regulatory variant interpretation. A liver-atlas metadata audit is available as a secondary resource. No model predictions, expression-matrix analyses, new gene discoveries, treatment effects, or personal risk estimates have been produced.

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

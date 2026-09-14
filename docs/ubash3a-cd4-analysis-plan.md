# UBASH3A CD4 experiment: proposed analysis and decision rules

September 14, 2026. **Design only: no experimental samples or results.** The RNA identities and comparisons below are specified now. The complete confirmatory plan is not frozen: the laboratory must qualify the method and complete the fields in the accompanying design record. The four-donor pilot estimates feasibility and variability; it does not count as independent confirmation of decisions made using its data.

## Measurement identity

Use the assay-target table's exact one-based exon coordinates and retain complete raw molecule evidence. Keep transcript structure, splice event, allele, cell state, biological donor, independent culture, RT preparation, library and molecule-family identifiers separate. An RNA must have all relevant junctions linked on the same original molecule; missing upstream or terminal information gives an unresolved or partial identity. Alignments must be checked against GRCh38 with alternate structures included. Preserve genuine novel structures separately from reference assignments.

Local +29 reads are a useful screening endpoint but cannot be pooled into a full-transcript decay curve. Do not infer the normal RNA's rs1893592 genotype from RNA at this excised position. Separately obtained donor genotype may describe the sample, but neither allele-specific decay nor a genotype association is a primary endpoint here. Repeated technical measurements, cells and UMIs are not independent donors.

## Estimands

For each qualified mature RNA j and donor d, estimate its decay constant **k(j,d,active)** and **k(j,d,suppressed)** in inverse hours. A single-exponential model has labeled abundance L(t) = L(0) exp(−kt), after validated recovery and cell-number corrections. Inspect its fit and the full time course. If residual labeling, precursor maturation or division is material, use the prequalified expanded kinetic model or report apparent turnover instead of claiming a degradation constant. Half-life is ln(2)/k only when k is positive and estimable.

Define the intervention-associated decay component:

**N(j,d) = k(j,d,active) − k(j,d,suppressed).**

The primary donor contrast is **D(d) = N(C29,d) − N(B29,d)**. Positive D means the canonical-derived +29 RNA has a larger suppression-associated reduction in decay than the boundary RNA. Estimate the donor-population mean D and its two-sided 95% interval, propagating curve uncertainty and preserving donor pairing. This contrast measures differential intervention dependence; attribute it specifically to NMD only with successful pathway controls and an independent perturbation.

Also report the active-condition half-life comparison of B29 and C29. Positive D alone does not prove that B29 has the longer untreated half-life, because their non-NMD turnover can differ. A broad claim of preferential baseline survival requires evidence for that longer baseline half-life as well.

Predeclared secondary measurements are the equivalent U29−B29 intervention contrast, the independently measured B0/C0/U0 curves, and log fold changes in abundance within each identified RNA. U29 shares B's reference upstream exon structure, whereas C29 contains an additional upstream exon. The U comparison informs this structural ambiguity, but cannot replace a missing primary C comparison. T29, the already NMD-annotated ENST00000635325.1 family and other observed structures are exploratory. A causal claim about one downstream junction requires a matched sequence intervention; these natural-transcript comparisons alone do not isolate that feature.

## Calibration and confirmatory freeze

Before confirmatory donor data, record the exact molecular-counting/normalization method, age-label classifier or enrichment recovery, limit of blank/detection/quantification, cross-assignment matrix, minimum informative molecule count, RNA/viability/purity thresholds, perturbation acceptance rule, time schedule, kinetic model, handling of non-exponential curves, meaningful positive margin δ in inverse hours, donor-level uncertainty method, desired power/precision, and sample-size calculation. The pilot may refine these fields; every refinement must be recorded and frozen before new confirmation samples are measured.

Use an assay-appropriate likelihood for the actual measurement unit. Calibrated molecule counts, digital-PCR concentrations and label probabilities have different errors. They must not be fed interchangeably to an unweighted log-linear fit. Paired donor estimates and technical uncertainty must be retained. The final method should be checked by simulation and calibration controls once the measurement format is known; no result-generating analysis program is claimed to be validated now.

Freeze one primary cell state and one primary contrast. Prefer the resting state if it qualifies; if feasibility requires the predeclared activated state, record that choice and the restricted estimand before confirmation. Do not select a state, dose or time window for the largest target effect. Orthogonal perturbations and secondary endpoints require separate uncertainty reporting and an explicit multiplicity plan before confirmatory claims.

## Interpretations

| Observation | Permitted conclusion |
|---|---|
| Only a local +29 junction is resolved | Exact local splice event detected; complete RNA unidentified |
| T29 is verified but B29 is not | Exploratory alternative-ending RNA supported; original boundary hypothesis untested |
| B29 is not detected with qualified recovery | Not detected in the tested state within the assay's stated sensitivity; no universal absence claim |
| Primary comparator missing, intervention controls fail, or decay curves lack information | Primary contrast not estimable; explain the failed requirement |
| RNA abundance rises, but no qualified decay measurement exists | Abundance response; cannot distinguish altered production/splicing from decay |
| D has a sufficiently precise positive interval, but only one perturbation was used | Differential dependence on that intervention; NMD mechanism remains to be confirmed |
| Primary interval supports a positive D, independent NMD perturbation agrees, and controls pass | Support for differential NMD-dependent decay in the tested donor population/state |
| D interval includes zero and also meaningful positive effects | Inconclusive; a nonsignificant result does not establish resistance or equality |
| Upper interval bound is below the predeclared positive δ with valid controls and adequate precision | Evidence against an advantage of that specified size in this setting; not proof of exact equality |
| Altered peptide is detected | Translation evidence for compatible RNAs; source-RNA identity and protein persistence require separate evidence |

No zero-count pseudocount should manufacture a half-life or effect. Retain missing data and failed samples explicitly. A zero in a validated measurement model can contribute a detection bound; it cannot be silently removed or replaced by a small positive value. Technical failures and biological non-detection need separate codes.

## Records and delivery

The supplied CSV templates are empty schemas, not example observations. `sample-manifest-template.csv` tracks sampling and process identity; `molecule-evidence-template.csv` joins raw molecules to full structures; `decay-measurements-template.csv` records abundance, measurement unit, recovery, age-label information and quality state. Numeric fields are blank until measured. Keep actual donor genotype and sample-level data in authorized study storage, not in this public-input planning package.

Request raw sequence/signal files where applicable, sample and library metadata, molecule-level evidence, complete end/junction calls including excluded reads, calibration and pathway-control measurements, time-course measurements, software versions, parameter files and file hashes. An abundance spreadsheet without linkage and calibration cannot resolve the original failure of identifiability. Report the feasibility and confirmation phases separately, including every outcome that prevents a conclusion.

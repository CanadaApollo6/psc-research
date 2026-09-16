# BACH2: fixed naïve-CD4 released-RNA observation check

**Complete.** The one frozen endpoint, full independent source/native audit, exact scientific replay and full aggregate publication have passed their separate checks. This fixed branch is closed; no further scientific analysis is planned.

## Result

The SNP-minus-noSNP difference in mean donor log2(1 + CPM) is **−0.4441001**, with nominal actual-df Welch **95% CI [−1.0942751, +0.2060748]** and two-sided **P = 0.1432691**. This is an imprecise negative point estimate, not a nominal P < 0.05 discovery or a demonstration of equivalent RNA levels.

All eight leave-one-donor-out **points** remain negative, ranging from **−0.6315878 to −0.2975364**, with zero strict sign reversals. This shows direction stability to the specified single omissions in this sample. The omissions are overlapping sensitivities, not independent replication or proof of a population difference. No omission CI or P was computed.

| Source group | Donors | Mean donor log2(1 + CPM) | Sample SD |
|---|---:|---:|---:|
| SNP | 4 | 5.26602480 | 0.416728747 |
| noSNP | 4 | 5.71012493 | 0.313682619 |

| Fixed contrast quantity | Result |
|---|---:|
| SNP − noSNP | −0.444100136 |
| Standard error | 0.260796680 |
| Actual Satterthwaite degrees of freedom | 5.57342779 |
| Nominal pointwise 95% Welch interval | [−1.09427507, +0.206074803] |
| Nominal two-sided P | 0.143269126 |
| Inference status | `available_nominal_Welch` |

Displayed decimals are rounded for readability. Computational precision is not biological precision. The full source aggregate retains the original numerical strings; this contrast is **not a raw RNA fold change, absolute RNA amount, transcription rate or protein measurement**. There is one prespecified gene/population/contrast, not a genome-wide or project-wide error-control claim.

## What was fixed before values

The [prospective claim appraisal](psc-bach2-followup-design-review.md) and [protocol](../docs/psc-bach2-naive-rna-protocol.md) selected one limited same-cohort check. The [root freeze](../config/psc-bach2-naive-rna-execution-freeze.json) was committed and pushed in **`a1b0a6e1d880f9b97c49bdbc0358aef7db408a47` before extraction**. No source, gene, population, normalization, statistical method, numerical criterion, donor exclusion or clinical join changed after values were accessed.

- Use only the eight qualified [E-MTAB-14013](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14013) donor units from [Poch/Bahn et al.](https://doi.org/10.1016/j.xcrm.2024.101620): four source `SNP` and four `noSNP`, all male PSC cases. The labels do not independently verify allele, build, ancestry or a causal nucleotide effect. Groups are unpaired; **donors, not cells, are replicates**.
- Start with the same 55,460 retained annotation rows. Select the exact C0/0 `CD4+ TN RTE` and C1/1 `CD4+ TN mature` union: **9,320 + 7,339 = 16,659 cells**, both states in every donor. Do not add the 3,855 source-only cells or guess historical filtering.
- Join sample plus barcode, stripping only that sample's prefix. Bare barcodes collide across samples. Preserve original source labels and annotation crosswalks; no new clustering or activation-state reconstruction is performed.
- Use the unique source feature **`ENSG00000112182 / BACH2 / Gene Expression`**, at **full source feature row 12,314, one-based** in all eight axes. This is not a genomic coordinate. Cached Ensembl 116 reciprocal identity checking does not reconstruct the executed source annotation or variant alleles.
- Pool selected cells within each donor. Let B be their released BACH2 RNA units and R their released units over **all 36,601 RNA rows**, excluding **38 ADT rows**. Calculate `Y = log2(1 + 1000000 * B/R)` and give the four donors in each group equal weight. Do not average cell-level logs, depth-weight donors, downsample or balance C0/C1 fractions.

The [metadata](psc-bach2-input-qualification.md), [first archive](psc-bach2-count-qualification.md) and [full cohort](psc-bach2-cohort-qualification.md) qualifications preserve source URLs, retrieval hashes, feature axes and unresolved history. Released source-count-format units are **not certified untouched original UMIs**. Cellranger 6.1.1 MEX headers versus the Methods/IDF 3.0.2 and reference history remain unresolved. The original-UMI/count-likelihood gate is unpassed; no NB or count likelihood was used.

Figure S2B's broad naïve-CD4 compartment motivated this check. Historical normalization, pooling and testing units were insufficiently specified for exact figure-method replication. Published results were available during design. This is not untouched independent validation, the Figure 3 activation-subcluster reconstruction, the Figure 5C healthy-donor transfection experiment, or an all-CD3 analysis. E-MTAB-14103 contributes no genotype-linked replication cohort.

## Complete fixed diagnostics and gates

All eight donors, fixed source/membership axes and positive RNA denominators passed. Zero B was permitted and would have been retained; **zero of the eight donors had B = 0**. No abundance or detection filter was applied. Across all selected cells, **2,408 of 16,659** had positive released BACH2 counts. This is cohort-wide QC, not an additional genotype-group detection endpoint, a cell-level inferential sample size or evidence that undetected cells lack protein.

All eight omission points were defined: **8 negative, 0 positive, 0 zero, 0 strict reversals** relative to the negative full point. No donor was removed from the primary comparison. No omission P/CI, alternative gene/state, permutation or normality-based test switch was added. Donor B, R, Y, source joins, state fractions and individual omissions remain private; only the predeclared aggregate summaries are published.

The working Welch model assumes independent donors and approximately Gaussian transformed relative RNA within groups. Four donors per group cannot demonstrate that calibration, representative sampling, clinical exchangeability or random genotype allocation. The interval includes negative and positive differences on this specific score scale. No independently justified equivalence margin was specified, so neither P > 0.05 nor the omission results support equivalence or biological negligibility.

## Verification and replay

Root repeated **167 pre-effect native tests**: 80 author, 46 independent-method and 41 separate execution-auditor tests, with zero failures/errors/skips. Guarded actual default author/auditor checks attempted no real quantity opens or result writes. The complete source selection, 55 annotation/receipt hashes and native runtime—including **2,362 installed NumPy/SciPy file hashes**—were checked. Author, independent reviewer and root accepted final V7/candidate-v4. A final report-only bookkeeping update was explicitly linked to the old/new hashes; accepted operative code, test receipts and handoff were not silently replaced.

The one real primary completed the exact four-file native inventory, with completion published last. Its completion SHA is **`850ff604b819fb0aae660fabc9a42248d511a5a334b07764afdb46d3dd961c4f`**. The separately sealed auditor authenticated all completed payload hashes before numerical parsing, independently streamed **all eight complete MatrixMarket sources**, reconstructed every selected-cell B/R/detection quantity and checked all donor Y, group means/SD, sole contrast/SE/df/CI/P/status, all eight omission points and the entire public projection/QC/provenance. The audit passed. All **102 outer metadata/code/protocol/runtime/source-receipt/native-authorization seals** and all **16 actual matrix/archive hashes** match the frozen inputs.

The independent reference uses 4,000-digit direct logarithms/exact rational point calculations and centered moments. Nonzero quantities require relative agreement within 1e-10, with exact zero/sign/constancy/availability and no absolute epsilon; standardized CI bounds require relative 1e-9. The tail reference uses a positive 100-digit Decimal incomplete-beta series with float log-gamma normalization. Critical values use SciPy `betainc` inversion with `brentq`; `stdtr/stdtrit` only corroborate native availability. Shared installed libraries remain a limit to software independence. These computational checks do not establish biological precision, small-sample calibration or clinical validity.

Mechanical replay changes only the authorized fresh destination and its prescribed provenance. **`measurements.private.json` is byte-identical**. In `public-aggregate.json`, the only changed leaf is **`provenance.authorization_sha256`**; every scientific value, status, QC field and other provenance digest is exact. Native authorization differs only by `output_directory`; completion changes exactly the resulting destination/auth/public-payload bindings. No broad provenance subtree was excluded. The replay and independent audit are computational checks, not independent biological observations.

The earlier method review's access to a closed cohort-qualification JSON containing historical donor technical totals remains disclosed. Those totals were not extracted, compared, selected on or analyzed. The reviewer also reported authorized opaque hashing of eight compressed archives without decoding members. Production preparation did neither, and no selected-population BACH2 quantity was measured during method preparation. The final metadata comparator uses individual integrity/acquisition receipts; this does not erase the historical access.

## Complete aggregate package and reproduction

- [Exact public aggregate](../data/derived/psc-bach2-naive-rna/public-aggregate.json), [one-row endpoint](../data/derived/psc-bach2-naive-rna/endpoint.csv), [two group rows with sample SD](../data/derived/psc-bach2-naive-rna/groups.csv), and [all fixed influence/QC fields](../data/derived/psc-bach2-naive-rna/influence-qc.json).
- [Independent full-source audit](../data/derived/psc-bach2-naive-rna/aggregate-verification.json) and [scientific replay](../data/derived/psc-bach2-naive-rna/mechanical-replay-verification.json): exact approved receipt copies, not new audits by the reporter.
- [PNG](figures/psc-bach2-naive-rna.png) / [SVG](figures/psc-bach2-naive-rna.svg), [manifest](../data/derived/psc-bach2-naive-rna/manifest.json), [completion](../data/derived/psc-bach2-naive-rna/completion.json), and [root presentation verification](psc-bach2-naive-rna-presentation-replay-verification.json).
- [Independent aggregate interpretation review](psc-bach2-naive-rna-results-review.md), [reporter design](psc-bach2-naive-rna-reporting-design.md), and [exact restoration/reproduction guide](../docs/psc-bach2-naive-rna-reproduction.md). A clone alone lacks required ignored inputs. Never overwrite existing runs or silently repin missing sources/runtime.

Root repeated all **27 reporter tests**, checked all **49 CSV scalar fields** against source values, verified the full JSON/influence/QC projection, and inspected the actual PNG. A separate source-only reporter review found no blocker. Root presentation replay checks **24 actual hashes**: all **ten payloads plus manifest** are byte-identical; only the completion's origin preview directory and consequent hash change. This is distinct from the native scientific replay described above. After separate root approval, **all 12 published files** match the approved preview bytes. No scientific statistics were recalculated for presentation.

Aggregate n=4 output is not formal differential privacy or a guarantee against inferential disclosure. The full native repository run executed **1,202 tests**, with **zero failures and the same seven historical errors** from missing ETS2/GSE84161 caches and R runtime. This is not a clean whole-repository suite. Frozen science was not altered to repair those historical dependencies.

## Interpretation and finite close

The observed point estimates lower transformed relative BACH2 RNA in the source SNP group, but the model-based interval remains imprecise and crosses zero. The same-cohort check therefore neither establishes a genotype-associated population shift nor establishes comparable/equivalent RNA levels. It cannot measure translation efficiency, protein abundance, mature miR4464, RNA stability or a specific cell-intrinsic regulatory mechanism. No named BACH2 ADT endpoint was deposited; this is not evidence of biological protein absence.

Pooling retains C0/C1 state mixture and RNA-content weighting. Clinical IBD/treatment/severity imbalance, ancestry/LD, batch, retention and annotation limitations remain possible explanations; this analysis does not identify which operated. There is no verified archive-to-clinical-row bridge. A slash was not converted to no disease/treatment, birth year was not treated as sampling age and no clinical covariate was imputed. Study donor labels support an independence assumption, not an independent identity/kinship investigation.

An RNA shift could qualify a broad RNA-similarity reading, but this uncertain relative-RNA observation cannot prove or refute the separate protein/translation experiments. No causal allele, disease mechanism, treatment benefit or clinical recommendation follows. The one fixed endpoint and all its diagnostics close here regardless of direction, P, precision or influence. No result triggers extra genes, states, cohorts, normalization choices or the [parked reserves](psc-remaining-frontier-appraisal.md).

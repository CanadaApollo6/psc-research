# GSE84161: a useful MEK comparison with source details to resolve

Investigated September 12, 2026. **We recovered and verified all 30 original monocyte array files, linked them to the original paper, and prepared a separate array analysis framework. We have not calculated an ETS2-program response in this dataset.**

The study is Senger et al., *The kinase TPL2 activates ERK and p38 signaling to promote neutrophilic inflammation*, Science Signaling (2017), [PMID 28420753](https://pubmed.ncbi.nlm.nih.gov/28420753/). Its [GSE84162 umbrella accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84162) contains [GSE84161 monocyte arrays](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84161) and GSE84153 neutrophil RNA sequencing. The July 2019 GEO release date is distinct from the paper's April 2017 publication.

## What the study adds to our question

The authors compared inhibition of MEK with inhibition of TPL2, an upstream signaling kinase. They reported ERK-pathway modulation in monocytes and a cell-dependent TPL2/p38 branch in neutrophils. These are published findings, not results of our reanalysis. [Primary article text](https://www.researchgate.net/publication/316231923_The_kinase_TPL2_activates_ERK_and_p38_signaling_to_promote_neutrophilic_inflammation).

For our ETS2 work, the useful feature is the six-arm monocyte design: vehicle, MEK inhibitor and TPL2 inhibitor, each with and without LPS stimulation. It lets us ask whether the fixed ETS2-dependent gene program responds alongside, or differently from, broader inflammation and stress programs. Comparing two interventions is informative, but neither directly manipulates ETS2. Any later response would remain comparative pharmacology in briefly cultured blood monocytes, without establishing ETS2 mediation or a treatment effect in PSC liver.

## What we verified

| Item | Finding | Consequence |
|---|---|---|
| Biological design | GEO states five donors; all six arms appear in each of five numbered replicate groups. | Thirty arrays represent a candidate five-donor design, not 30 independent people. |
| Sample correspondence | The public monocyte records have replicate numbers, with no explicit donor key; the source analysis includes donor adjustment. | Pairing is strongly supported but remains unconfirmed under our prospective donor rule. Only anonymous group correspondence is needed. |
| Original files | All 30 advertised CEL files are present in a 146,739,200-byte archive. Every gzip stream passes integrity checks; all headers identify the expected chip and dimensions. | Original array preparation is feasible. File integrity does not establish biological or array-signal quality. |
| Repeated records | The 30 monocyte records and ordered feature identifiers agree between the subseries and superseries. | The superseries contributes no additional monocyte replication. Its 20 neutrophil records remain separate. |
| Deposited expression table | Each array has the same 19,944 unique probe-set identifiers. GEO describes RMA followed by highest-variance probe selection per Entrez gene. | Those rows are an already selected source table; we will use original CEL files for our new calibration. |
| Annotation | The full GPL570 snapshot has 54,675 probes. Of the deposited 19,944, 267 have no Entrez mapping and 414 map to multiple IDs in this snapshot. | Historical annotation or selection details are needed for exact reproduction; the current mapping cannot be silently substituted for the original one. |

The metadata and identifier checks are reproducible in the [metadata audit](gse84161-metadata-audit.json), [sample map](../data/derived/gse84161-samples.csv), and [raw-array inventory](../data/derived/gse84161-raw-arrays.csv). We decompressed raw bytes for checksums and headers; no probe-intensity values were interpreted.

The annotation details matter: GPL570 was last updated December 14, 2020, while its individual annotation rows are dated October 6, 2014. These dates do not identify the authors' historical annotation version. Restricting the present full platform to probes with one unique Entrez mapping yields 41,834 probes representing 20,486 genes, before the planned Ensembl crosswalk. That is an annotation inventory, not measured program coverage. The deposited table is labeled RMA signal; donor residualization is described for heatmaps, and should not automatically be attributed to the deposited values.

## A concrete source discrepancy

GEO gives **2.8 µM TPL2 inhibitor and 0.31 µM MEK inhibitor**. The paper's Figure 2D caption gives **926 nM and 103 nM**, respectively—approximately three times lower. Neither accessible record resolves which concentrations generated these particular arrays. The article also leaves the chemicals labeled TPL2i and MEKi; cited inhibitor papers do not uniquely identify them. [GEO design](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84161), [original figure caption](https://www.researchgate.net/publication/316231923_The_kinase_TPL2_activates_ERK_and_p38_signaling_to_promote_neutrophilic_inflammation).

We preserve both dose statements and assign no guessed chemical identity. This limits dose and compound attribution, while still allowing preparation under the source's treatment labels. The publication and failed supplement retrievals are recorded in the [provenance review](../docs/gse84161-publication-evidence.md) and [literature manifest](../config/gse84161-literature-sources.json). No usable supplement or original analysis code was obtained in this focused search.

## The next analysis is now defined

The [array framework](../docs/gse84161-array-protocol.md) specifies fresh RMA normalization across all 30 CEL files; a pinned, unambiguous probe-to-gene mapping; and the median of all eligible probes for each gene. It carries forward the same nine programs, expression-matched controls, mapping thresholds and donor-level reporting rules. It does not reuse the source's highest-variance selection.

Seven contrasts cover stimulated drug responses, basal drug responses, stimulation itself, and TPL2-versus-MEK comparisons. These correspond to 35 candidate within-block comparisons if donor correspondence is confirmed; they are neither 35 independent donors nor completed results. Both sides of each contrast use the same control background. The [preparation record](../config/gse84161-array-preparation.json) pins the exact candidate sample pairs and separates completed checks from pending execution requirements.

On our side, the next work is to prepare the array software and a versioned Ensembl–Entrez crosswalk, then freeze the executable method before inspecting program effects. The existing local R 4.4.3 runtime works, but its array packages are absent. No new software was installed in this investigation. Exact donor-block confirmation remains necessary for the planned paired biological summaries. The compound and dose questions limit mechanistic attribution separately; they do not hold up software or annotation preparation.

We prepared [four specific source questions](../docs/gse84161-source-questions.md) covering donor groups, compounds/doses, historical annotation and deposited-value processing. They are unsent. The prior ETS2 results and both validation-rule snapshots remain unchanged. This investigation strengthens the next experiment we can perform with public data; it does not change Gary's care or establish a new treatment option.

## Reproduce the audit

The [study source manifest](../config/gse84161-study-sources.json) pins the two SOFT caches and raw archive, including retrieval precision and local checksums. Large source files stay ignored. From the repository root with the recorded caches present:

```bash
.venv/bin/python scripts/audit_gse84161.py
.venv/bin/python scripts/audit_gse84161_raw.py
.venv/bin/python -m unittest discover -s tests -v
```

The [verification record](gse84161-study-verification.json) records actual checks and preserved historical artifacts. These commands audit metadata, identifiers and byte integrity; they do not run RMA or calculate program scores.

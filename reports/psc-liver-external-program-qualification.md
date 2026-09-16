# External liver program source qualification

**16 September 2026 UTC. Source-only; bounded pass closed.** No liver expression, normalization, matching, scores, effects or hypothesis tests were used. No organoid-derived membership was selected. Frozen organoid files were not changed.

## Decision

| Requested source | Status | Complete membership |
|---|---|---:|
| Published human epithelial **IL17A-alone response** | **HOLD — not qualified** | Unavailable, not an empty list |
| Externally fixed fibrosis/ECM comparator | **Qualified as an unsigned ECM pathway** | 321 official MSigDB symbols |

The requested epithelial response endpoint is not ready. An ECM pathway does **not** fill that gap. No automatic pathway replacement, small hand-picked list, UP-only assumption, or reuse of the eight organoid interaction candidates is authorized. The root must decide whether to retain HOLD or authorize a later, separately bounded changed source endpoint. No additional search is part of this closed pass. This source blocker need not delay MacroMap.

## Qualified comparator: complete official ECM membership

- **Program:** `REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION`.
- **Identifiers:** MSigDB `M610`; exact source `R-HSA-1474244`; collection `C2:CP:REACTOME`.
- **Version:** MSigDB **2026.1.Hs**, linked to **Reactome release 95** in the official page's version history. The JSON response filename independently contains `.v2026.1.Hs.json`.
- **Organism:** `Homo sapiens`.
- **Source:** the complete [official gene-set page](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION.html) and its [linked JSON](https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION&fileType=json), both retrieved without login on 2026-09-16.
- **Definition:** exactly the **321 unique entries in the official JSON `geneSymbols` array**, in its original order. This is a versioned MSigDB gene-set projection, not a claim to have independently recovered the whole Reactome molecular-participant universe.

The complete HTML declares **340 original `Human_Ensembl_Gene_ID` identifiers mapped to 321 genes**. All 340 rows were recovered. Their 321 unique symbols exactly equal the JSON members. All original Ensembl IDs and row indexes remain attached to the corresponding symbol in the public CSV. Published NCBI Gene IDs remain strings. No current-reference, alias, locus or tissue reconciliation was performed.

There are **12 many-to-one symbol groups** containing 31 original rows, hence 19 extra source IDs beyond the 321 official members: ACTG1(2), CAST(2), COL11A2(5), CTSB(2), DDR1(5), ELANE(2), ITGA9(2), MMP11(2), PRSS1(2), PRSS2(2), TNXB(3), TRAPPC4(2). These are source-published mappings, not independent genes or locally chosen merges. No repeated original ID, missing original ID/symbol/NCBI ID, or conflicting source assignment was found. This does not establish unambiguous mapping to a current reference or to liver features.

### Directions, weights and limits

Every member is **`UNSIGNED`**, with **membership weight 1 per official symbol**. In the CSV, `published_direction=UNSIGNED` is our status label for the absence of published perturbation-response directions, not a literal publisher direction or inferred effect sign. Multiple source IDs do not multiply a member's weight. This weight is a uniform membership convention, **not** a published positive treatment-response coefficient. No expression-based or disease-learned weights were estimated. A score formula and any tissue mapping/missingness policy still need a separate root method freeze.

This is non-perturbation pathway membership. “Constitutive” here means that it is not defined by a cytokine perturbation; it does not mean ubiquitous or stable expression. It is **not a directional fibrosis-response signature, fibrosis-specific program, PSC-specific program, epithelial program, or IL17A response signature**. Treatment identity, dose, duration, experimental donor count and differential-expression threshold are not applicable to this pathway source and were not invented. The source does not declare an Ensembl annotation release or genome build. Allele conventions do not apply.

The original denominator for future mapping is **321 official members**, not 340 source IDs or a successfully mapped subset. Keep unresolved members explicit. Missing mappings must not become zero expression or gene absence. No liver coverage, activation, fibrosis effect or clinical conclusion has been calculated here.

## Epithelial IL17A response: specific missing source inputs

The pass identified **Nograles et al. (2008)**, *Th17 cytokines interleukin (IL)-17 and IL-22 modulate distinct inflammatory and keratinocyte-response pathways*, PMID **18684158**, PMCID **PMC2724264**, DOI [10.1111/j.1365-2133.2008.08769.x](https://doi.org/10.1111/j.1365-2133.2008.08769.x). The captured Europe PMC bibliographic response and abstract describe keratinocyte cytokine-response profiling in a human-skin context. They are not a complete perturbation signature or a verification of exact IL17A-alone culture conditions.

The missing inputs are:

1. The **complete published member list and all its published directions**, not selected genes mentioned in an abstract.
2. Explicit IL17A-alone treatment identity, comparator, epithelial culture model, dose, duration and experimental replication.
3. Original identifier namespace/version and the published criteria defining that complete list.

No missing dose, donor count, downregulated arm or gene membership was inferred. A skin keratinocyte response, even if subsequently qualified, would not be a cholangiocyte response or an independent liver disease-control cohort. The original published perturbation would remain discovery evidence, not biological replication of this project's liver question.

### Observed access results

- `NOGRALES_IL17_UP` and the final candidate-name lookup `CHIRICOZZI_IL17A_RESPONSE` each returned the official MSigDB **Gene Set Not Found** page. These were source-name hypotheses, not accepted programs. This does not prove that those studies have no published lists or that any historical set was definitely deprecated.
- The Nograles PMC article URL returned **HTTP 200 with a reCAPTCHA challenge**, not article text.
- The Europe PMC article page returned **HTTP 403 with a Cloudflare challenge**. Its public `fullTextXML` endpoint returned **HTTP 404 with a zero-byte body**.
- An ordinary MSigDB keyword-search URL redirected to a login page. No login, registration, credentials, challenge solving or access-control bypass was used.
- Two initial acquisitions failed after HTTP at a local output-path error, before their bodies and status could be saved. No outside-owned-area file was created. The path calculation was corrected and the owned directory is now checked before HTTP. The exact failed response sizes/status and any redirect history remain unavailable; they were not reconstructed from later responses.

Thus the result is a **source-access/completeness blocker for this bounded pass**, not evidence that epithelial IL17A responses do not exist. No official IL17 pathway has been qualified as a substitute either. The unavailable endpoint has `member_count: null`, no directions or weights, and no rows in the membership CSV. It must not be treated as a zero-member/zero-score program.

## Acquisition and reproducibility

Known caches were checked first. The existing ETS2/eight-comparator sources did not supply an IL17 epithelial or ECM definition. The new pass used **10 main logical acquisitions plus 2 ECM acquisitions**, with **4 recorded redirect hops**. Initial failed redirect histories are unknown. These are logical acquisition limits, not an exact count of physical HTTP hops.

- **330,514 retained response-body bytes**, including challenges, error pages and searches; **200,901** belong to the two complete ECM sources.
- **1,150,648-byte conservative accounting charge**, including the full **820,000-byte cap reserve** for the two lost initial bodies and **134 declared redirect-body bytes** not retained.
- Ceiling: **5,000,000 new public metadata/membership bytes**. The accounting charge is below that ceiling. This is **not an exact wire-transfer audit**: the initial failed responses and redirect histories were lost. Retained captures reached EOF; no truncated body is used as membership evidence.
- No new expression data, private API, service enrollment, contact, dependency installation or Git operation occurred in this task. Acquisition is closed.

The [source configuration](../config/psc-liver-external-program-sources.json) pins versions, URLs, retrieval times, source byte counts/hashes, directions, weights, projection rules, limits and the IL17 HOLD. The [membership CSV](../data/derived/psc-liver-external-program-membership.csv) has 321 rows and preserves all 340 source IDs. Its SHA-256 is **`8cf6573a97d01cc8e83558634bf62b6d48d9be015e67f50d15868071c18b7202`**; size **125,834 bytes**.

| Complete membership source | Retained bytes | SHA-256 |
|---|---:|---|
| Official HTML, retrieved 2026-09-16T02:16:29.722831+00:00 | 197,895 | `0bce8c070ad8e24961cd22fabd57da8cf26161f7ece5d3da49715812ac1895a6` |
| Official JSON, retrieved 2026-09-16T02:17:25.700864+00:00 | 3,006 | `2046e084f1c36c74cf107c8448d19bb4f47bb2bf37b45b2388fd3d1b837b6540` |

Ignored full request receipts, source excerpts, original 340-row table, projection audit and failed-capture log remain under `work/psc-liver-external-programs/` and `data/raw/psc-liver-external-programs/`. Original source bytes remain unchanged; all parsing and semantic access judgments are separate.

### Offline checks

```text
.venv/bin/python scripts/qualify_psc_liver_external_programs.py --build
.venv/bin/python -m unittest discover -s tests -p test_qualify_psc_liver_external_programs.py -v
.venv/bin/python scripts/qualify_psc_liver_external_programs.py --verify
```

The qualifier has **no network, liver mapping, expression or scoring mode**. It fails on incomplete rows, changed source/version pins, missing/extra/repeated symbols, repeated source IDs, conflicting source mappings, invented directions, or an IL17 empty-list/substitution attempt. `--build` writes only the owned membership CSV. `--verify` checks pins and exact CSV reconstruction read-only. All **27 synthetic tests pass**. The full repository suite ran **664 tests with zero failures and seven errors**, all at the same unrelated missing historical ETS2/GSE84161 cache or R-runtime paths. No historical inputs were repaired. The complete log and error labels are retained in the source config and ignored work directory. The independent final source audit passed **47/47 checks with no issues**. It compared all 321 CSV symbols/order and all 340 original IDs/indexes/NCBI assignments to the raw official sources, plus **4,494 per-member constant/pin comparisons**. It used its own streaming HTML parser, without importing or running the producer, and made no network request. It also reconciled all acquisition charges and the IL17 HOLD. This is source/serialization verification, not biological validation. The audit did not rerun or certify the full test suite. Snapshot hashes in its receipt can precede the appended validation summaries in this report/config; source bytes and membership remain unchanged.

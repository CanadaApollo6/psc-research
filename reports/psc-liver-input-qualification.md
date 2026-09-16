# PSC liver input qualification

## Decision

This is input qualification only. No gene contrasts, target/program scores, models or biological effects were computed. The organoid handoff remains closed and unchanged. The root analysis owns any later panel and analysis freeze.

| Resource | Qualified interface | Decision |
|---|---|---|
| GSE303271 | Complete 41,096-label by 64-sample submitted count matrix; exact GEO/library joins | Count input usable with source-annotation and clinical-metadata limits |
| GSE159676 | Complete 17,046-feature by 33-array matrix; historical platform join; paired-biopsy structure | **Hold expression analysis:** mixed-scale source literals and unresolved cohort labels |
| NoPSC liver atlas | Bounded browser label/schema metadata only | **Park:** no qualified original-count/feature/barcode/donor export |

## GSE303271: pediatric whole-biopsy bulk RNA

- The full raw-count download is 3,021,812 bytes. Its 2,630,144 values are finite, nonnegative, exact integer literals. All 41,096 submitted labels are unique; all 64 columns join exactly to GEO titles and explicit library names. There are 64 distinct GSM, BioSample and SRA accessions, no repeated column identifiers, and no exactly duplicated numeric columns. These observations alone do not prove independent patients.
- Original groups are **17 PSC, 17 ASC and 30 AIH**. The paper describes 64 patients, including 34 labelled PSC, while the GEO series uses PSC n=34. The sample-level deposit splits these into PSC and ASC. Retain ASC separately: no explicit ASC recoding rule was recovered.
- Tissue is bulk RNA from cryopreserved whole-biopsy liver tissue obtained during clinically indicated diagnostic biopsies. This is not a whole organ or sorted portal fibroblasts, macrophages or cholangiocytes. No paired-region scheme is stated for this bulk deposit; the paper's separate two-patient spatial assay is not extra bulk replication.
- Keep all 998 all-zero rows in the qualified released universe. There are 998,756 zero cells overall. Per-library count sums range from 40,155,767 to 350,237,961; median 143,139,304.5. These are sums of deposited assignments, not independently verified read totals. No library or feature was dropped, rescaled or normalized.
- GEO names **GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4**. The paper's DESeq2 1.26 normalization/variance stabilization is downstream analysis, not the downloaded raw-count quantity. The GTF/provider/release/checksum, assignment/strandedness/multimapping options and executed command are unavailable in inspected sources. The source read-length/depth wordings also differ; neither is silently corrected.
- The 41,096 row labels are preserved as submitted. This is the complete **released** feature universe, not a reconstructed annotation universe or a current validated gene-ID mapping. No aliases were silently updated.
- Exact sample-linked age, sex, fibrosis/stage, treatment, IBD status and assay batch were not recovered. The main paper mentions clinical measurements, but linked figure-supplement/workbook metadata did not provide a verified sample crosswalk. Uninspected figure sheets are not evidence that such information cannot exist. Cohort summaries, library-number order and submission dates were not assigned as individual covariates.

**Feasible later comparisons:** PSC versus AIH; ASC versus AIH; PSC versus ASC, with the original labels fixed. These are pediatric disease-control comparisons. There is no healthy-control group, no stage-adjusted interface without additional exact covariates, and no cell-specific causal interpretation. Sample independence remains supported by the paper's 64-patient statement rather than a separate person/visit crosswalk. No comparison was run.

## GSE159676: adult liver-tissue arrays

### Samples, donors and source conflicts

The complete 2,252,130-byte series matrix contains 33 exact GSM-matched records: **12 PSC, 7 NASH, 3 Primary biliary cirrhosis, 3 Autoimmune hepatitis, 1 Haemochromatosis, 1 Alcohol related and 6 Liver tissue healthy**. These are source labels, not a corrected diagnosis table.

1. **PSC is six patients with two biopsies per patient, not 12 independent patients.** The public HAL author proof states this explicitly (physical PDF p44, printed p42, Supplementary Methods; repeated on physical p101). GEO titles form six complete numeric a/b blocks. The title-based blocks are consistent with that design, but a separate GSM-to-person key and the anatomical meaning/date of a/b were not recovered. Do not label them as named regions or independent replication.
2. Both Methods and GEO prose describe **five** normal/control patients, but the deposit has **six** control-style records. There is no recovered explanation for the extra record, duplicated donor, or historical exclusion. All six records remain preserved; no control was dropped.
3. The HAL proof reports **PBC 2 plus sarcoidosis 1**; GEO instead labels **Primary biliary cirrhosis 3** and has no sarcoidosis label. This is a source-version conflict, not permission to relabel any sample. HAL is an author proof, not silently represented as the final publisher supplement.
4. Microarray RNA comes from frozen liver tissue obtained from explants or diagnostic biopsies. The human array Methods do not describe sorting or microdissection. They are distinct from the paper's portal-fibroblast experiments. The control tissue is tumor-free liver from patients with colorectal cancer metastases, not healthy-volunteer biopsies. Patient-specific procurement, region, fibrosis/stage, age, sex and batch are not linked to the deposited arrays. The separate French RT-qPCR cohort's fibrosis strata must not be transferred to them.
5. The Norwegian PSC biobank is named. Overlap with the later Norwegian atlas remains unresolved, not proven absent.

### Numerical integrity: analysis hold

The GEO pre-table metadata declares **Quantile normalized** for all 33 arrays. The HAL Methods name `oligo`, `rma` and `hugene10sttranscriptcluster.db` (physical PDF p46, printed p44, Supplementary Methods). Exact versions, executed options, explicit deposited log base and export/selection history are not supplied. Conventional RMA defaults are not an executed run record. Do not log-transform this release again by assumption.

All 562,518 source values are finite and positive, but **4,829 integer literals are greater than 100**, spanning **3,912 features and all 33 arrays**. Their range is **671 to 2,958,313**. The other values range from **1.170617 to 12.92883**. The same 4,829 values exceed 20; the >100 flag is a numerical-integrity diagnostic, not a gene or sample exclusion rule. No target list was used.

All 660 available cells in the earlier cached GEO quick preview match the complete matrix numerically, including nine large literals. This traces the issue to released source values rather than this parser; it does not identify the historical cause or establish that the authors' underlying analysis was wrong. **No decimal insertion, division, imputation, extra log, sample exclusion or feature exclusion was applied.** A source-corrected processed release or a separately authorized raw-array reconstruction is needed before expression use. Gene effects cannot be rescued by simply treating these values as counts.

### Features and annotation

- All 17,046 submitted feature IDs are unique and exactly match GPL6244's complete **33,297-row** GEO annotation (dated **Aug 09 2016**). No numeric columns are exactly duplicated. Identical numeric feature vectors can occur; they were not assumed to be duplicate IDs or collapsed.
- The platform is `[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]`. Platform metadata records NetAffx build 35 in July 2016; that is an annotation release label, not a claimed genome build. No coordinates were inferred or lifted.
- Historical Entrez mapping has **16,806 one-ID features**, **115 multi-ID features**, and **125 without an Entrez ID**. The one-ID set represents **16,020 unique Entrez IDs**; **591 IDs have more than one measured feature**. No mapping ambiguity or multiple-feature aggregation was silently resolved.
- The historical GEO annotation is not claimed to be the exact `hugene10sttranscriptcluster.db` snapshot used by the authors. Current gene-ID reconciliation and the source rule selecting 17,046 of 33,297 platform rows remain unavailable. Missing platform measurements were not filled in.

**Conditional later comparisons:** PSC versus NASH is a disease-control interface with six paired-biopsy PSC blocks, not n=12 independent PSC. PSC versus control-style tissue additionally needs the five/six control reconciliation. Comparisons using PBC/other labels need their source-version conflict resolved. Every expression comparison remains on hold for numerical integrity and requires a frozen donor/block and feature-aggregation rule. No comparison was run.

## NoPSC atlas: one bounded linked-export check

The normal public browser serves search-label lists, annotation/display bundles and five labelled differential-expression spreadsheet links. Those spreadsheet payloads and all individual-gene endpoints were **not** requested. No original-count bulk inventory, stable feature/barcode axis or donor-to-observation join was established.

| Browser metadata | Nuclear | Spatial |
|---|---:|---:|
| Unique search-menu labels | 25,539 | 17,319 |
| Annotation-prefix observations | 72,990 | 45,998 |
| Complete categorical fields inspected | 2 | 6 |
| Declared whole annotation-bundle bytes | 1,314,400 | 1,749,037 |
| Retained prefix bytes | 146,317 | 276,568 |

The label lists include `Total RNA count` and `Total RNA features`, so their lengths are **not** verified gene-universe counts. The annotation streams were closed before the first continuous payload. Nuclear prefixes contain disease and cell-class labels; spatial prefixes contain disease, subtype, stage, location and spatial-class labels. No donor-ID or barcode field occurs in these prefixes, but uninspected tails are not evidence of absence. Observation row counts are not independent donor counts.

The paper reports **12 PSC + 4 cirrhotic-control nuclear donors**, overlapping **23 PSC + 7 cirrhotic-control spatial donors**. Modalities are not replication. The browser's measurement quantity is not established by the article's SCTransform description. Norwegian-biobank overlap with GSE159676 remains unresolved.

**Reopen only** for a normal public original-count release with assay/units, complete feature IDs, a stable barcode axis, anonymous donor mapping and file URLs/sizes. Report that inventory and obtain root authorization before bulk transfer. No usable bulk inventory was found, so no bulk request is pending. Exact endpoint URLs, sizes and partial-schema limits are in the [machine-readable result](../data/derived/psc-liver-input-qualification.json).

## Provenance, boundaries and reproduction

- Complete inputs, failed/size-refused requests, reused cache receipts and SHA-256 pins: `data/raw/psc-liver-qualification/`. Source joins, unchanged-value numeric arrays, full feature/annotation tables and integrity flags stay in ignored `work/psc-liver-qualification/`.
- New bodies read/retained total **23,735,752 bytes**. Charging the full advertised bodies of the two stopped atlas-prefix GETs gives **26,376,304 bytes**, below the 80,000,000-byte overall cap and each worker allocation. Oversize requests were stopped from headers without body reads. No FASTQ, CEL/raw-array archive, atlas matrix, researcher contact, account or environment change was used.
- Matrix gzip CRC/end markers, all numerical cells, sample joins, original labels, feature uniqueness and mapping ambiguity were checked. Synthetic tests cover truncation, malformed joins, duplicate IDs, missing/fractional/nonfinite values, float-rounding loss and bounded acquisition. An independent stdlib audit rechecks complete matrices and metadata without importing the qualifier. The source-level audit is separate from inference.
- Only lossless decompression, TSV/CSV parsing, exact metadata/annotation joins, documented a/b title parsing, and float64 copies for QC were applied. Source files/literals remain unchanged. Public outputs contain aggregate qualification only. No RNA normalization, gene aggregation, gene contrast, target/program score, model, effect or clinical recommendation was produced.
- Python/package versions, source URLs/retrieval times, Methods locations, all assay/annotation caveats, code hashes and ignored-artifact pins are in the [JSON qualification](../data/derived/psc-liver-input-qualification.json). All four closed versioned organoid pins and 67 preceding source/ignored-artifact pins remain unchanged; no organoid file was written.

Offline commands:

```sh
.venv/bin/python scripts/qualify_psc_liver_inputs.py --qualify
.venv/bin/python scripts/qualify_psc_liver_inputs.py --verify
.venv/bin/python -m unittest discover -s tests -p test_psc_liver_inputs.py -v
```

Acquisition is a separate, explicit `--acquire` action. Qualification and verification make no network requests. Historical unrelated repository tests may still require their own absent ETS2/GSE84161 caches or R runtime; this task does not alter those closed branches.

### Main primary sources

- [GSE303271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303271); [pediatric paper](https://doi.org/10.1172/jci.insight.199226), Methods/RNA-Seq acquisition and Data availability.
- [GSE159676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159676); [adult article](https://doi.org/10.1002/hep.32456); [public HAL author proof](https://hal.science/hal-03667547/file/HEP-21-2383.R1_Proof_hi.pdf), physical pp44 and46.
- [GPL6244 historical annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz).
- [NoPSC nuclear browser](https://riim-data.no/NoPSC_Liver_Atlas/SINGLE_NUC/), [spatial browser](https://riim-data.no/NoPSC_Liver_Atlas/SPATIAL/) and [atlas article](https://doi.org/10.1097/HEP.0000000000001432).

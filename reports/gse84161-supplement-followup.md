# GSE84161 supplement follow-up

**Investigated:** 2026-09-14. **Result:** both study inhibitor codes are resolved; anonymous donor correspondence and the concentrations used for the deposited arrays remain unresolved. This is a public-source provenance update, with no new program scores or treatment comparisons.

| Question | Current status | Evidence |
|---|---|---|
| Which TPL2 inhibitor? | **G-432** | Original supplement, Fig. S3A/E, PDF page 7. |
| Which MEK inhibitor? | **G-573** | Original supplement, Fig. S3E, PDF page 7, and Fig. S4A, page 9. |
| Do replicate numbers 1-5 each identify one donor across all six arms? | **Not explicitly confirmed** | Five donors and five complete groups support this interpretation, but no sample-to-donor key was located in the inspected material. |
| Which inhibitor concentrations produced the 30 CEL files? | **Conflicting source statements remain** | GEO and the article's array caption differ by approximately threefold; the supplement does not reconcile them. |

## Original supplement recovered

The publisher's actual public link is [aah4273_sm.pdf](https://www.science.org/doi/suppl/10.1126/scisignal.aah4273/suppl_file/aah4273_sm.pdf), with a lowercase filename. A normal browser download succeeded after the previous requests had returned HTML. The saved file is a valid **28-page PDF, 5,302,678 bytes**, SHA-256 `c477aab305df79ce53214491924baeb44333b290f6b0f0e1e277cab277841bd5`. The [new source manifest](../config/gse84161-supplement-sources.json) records the download receipt, exact source locations and failed attempts separately.

Fig. S3E explicitly associates TPL2i with G-432 and MEKi with G-573. Its caption continues on page 8 and describes potency measurements for the study's inhibitors. Fig. S4A independently names G-573. These direct labels supersede the compound-code uncertainty in the preserved [September 12 publication audit](../docs/gse84161-publication-evidence.md). They do not establish a vendor catalog number, batch or CAS identifier.

Two additional primary sources provide structure references:

- **G-573:** Hatzivassiliou et al. (2013), the paper cited by the original supplement, displays G573 in Supplementary Fig. 1A, **PDF page 22** of its [36-page supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature12441/MediaObjects/41586_2013_BFnature12441_MOESM60_ESM.pdf). The figure also depicts other MEK inhibitors separately. The corresponding reference records, [PubChem CID 46926364](https://pubchem.ncbi.nlm.nih.gov/compound/46926364) and [RCSB component 573](https://www.rcsb.org/ligand/573), explicitly cross-link and share InChIKey `ZEZHPEIEEFTILY-QMMMGPOBSA-N`. This verifies agreement between reference records; the drawing does not specify stereochemistry or establish the 2017 study batch.
- **G-432:** Wang et al. (2023) publishes a structure labeled G432 in [Figure 1 - figure supplement 2A](https://pmc.ncbi.nlm.nih.gov/articles/PMC10411973/#fig1s2). This later primary source supplies a structural reference. No new machine-readable structure or CAS assignment was made here.

## Array concentrations remain unsettled

The original statements are retained separately. Values below are all displayed in nM; GEO's original values are **2.8 µM TPL2i and 0.31 µM MEKi**.

| Source and experimental context | TPL2i / G-432 | MEKi / G-573 |
|---|---:|---:|
| GEO GSE84161 overall design; deposited monocyte array study | 2,800 nM | 310 nM |
| Senger et al. main article, Fig. 2D array caption; prior primary-text audit | 926 nM | 103 nM |
| Supplement Fig. S2E, pages 5-6; maximum doses shown in a 20-minute LPS phosphorylation assay | 926 nM | 103 nM |
| Supplement Fig. S5A, page 10; six-hour LPS monocyte cytokine assay | 926 nM | 103 nM |

Fig. S2E uses 100 ng/ml LPS; Fig. S5A uses 10 ng/ml. The latter shares the array study's nominal stimulation duration and LPS concentration, but measures cytokines rather than identifying the array samples. These supporting experiments strengthen the evidence that the lower concentrations were deliberately used somewhere in the study. **They do not determine which concentrations generated the deposited CEL files.** No dose correction was located in the inspected PDF.

The GEO statement comes from the checksum-verified **September 12 SOFT snapshot**, rechecked for metadata on September 14; it is not a newly downloaded GEO version. The main Fig. 2D statement is inherited from the [dated primary-paper text audit](../docs/gse84161-publication-evidence.md), whose author-upload host and retrieval limitations remain documented. No source value has been overwritten or selected as the definitive array dose.

Fig. S9, page 16, also provides supporting monocyte pathway-engagement evidence: ERK phosphorylation decreases with both inhibitors without corresponding p38 suppression in the shown assay. Its caption specifies at least 30 minutes of pretreatment, 100 ng/ml LPS for 20 minutes, and results representative of five experiments. This does not demonstrate target engagement in the same 30 array samples at six hours.

## Anonymous donor correspondence remains open

All **30** monocyte sample metadata records were rechecked. Each description contains a distinct source label, **SAM636535 through SAM636564**. These 30 labels are not reused across conditions and therefore do not supply a five-donor key. The records contain no sample-relation field providing a BioSample or donor cross-reference.

Replicate numbers 1-5 still each contain the six expected conditions, and the stated study design uses five donors with donor adjustment in the analysis. This supports five candidate blocks but does not explicitly certify donor correspondence. The [new metadata audit table](../data/derived/gse84161-supplement-donor-audit.csv) preserves every original title, description, SAM label and replicate number; its `confirmed_donor_group` column is empty. The existing sample manifest is unchanged. Donor references in the neutrophil proteomics supplement are not transferred to these monocyte arrays.

**The separate table archive was not obtained.** The publisher lists [aah4273_tables_s1_to_s9.zip](https://www.science.org/doi/suppl/10.1126/scisignal.aah4273/suppl_file/aah4273_tables_s1_to_s9.zip), displayed as 666.59 KB. Normal browser download attempts did not save it, and direct requests returned blocked navigation or HTTP 403 HTML. The PDF's page 27 legend describes Table S2 as gene-level fold-change and significance summaries, but the actual spreadsheet sheets were not inspected. A donor key or additional dose detail within that archive has therefore **not** been ruled out.

## Consequence and remaining questions

The chemical-code gap is closed. Before executing the planned paired program comparisons, the two remaining study-metadata questions are:

1. Does each replicate label 1-5 identify the same donor across all six treatment/stimulation arms, or is a different anonymous GSM-to-donor mapping required?
2. Were the deposited arrays exposed to 926/103 nM or 2,800/310 nM TPL2i/MEKi, or does concentration vary among samples?

The unopened table archive remains a possible public source for these answers. The [earlier source questions](../docs/gse84161-source-questions.md) are preserved as an unsent historical draft; their compound-identity question is now answered at the study-code level. Historical analysis-code provenance is unchanged and is separate from the fresh RMA procedure already completed.

The array locks, matrices, mapping gates and earlier results remain unchanged. This update does not turn MEK pharmacology into a direct ETS2 perturbation or establish a clinical effect. No researchers were contacted.

## Verification and reproducibility

The [verification record](gse84161-supplement-verification.json) checks all 11 usable source caches against their byte counts and hashes, PDF page counts and file signatures, the 30-row metadata extraction against the preserved sample table, the empty confirmed-donor fields, the G-573 reference-record cross-link and unit conversions. It separately records **79 unchanged previously pinned artifacts** and four intentional updates to repository navigation/evidence documents. Prior reports, analysis code and execution locks are preserved.

This was a documentation and source-metadata investigation. No import logic or analysis implementation changed, so the previous 144-test result is not presented as a newly executed test run. Raw articles and supporting media remain in ignored local storage; only small metadata, findings and provenance records are versioned.

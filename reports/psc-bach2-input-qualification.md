# PSC BACH2 input qualification

Date: 16 September 2026 UTC. **Metadata qualification only; no expression effects.**

## Decision

**E-MTAB-14013 has a usable source donor–genotype–cell-annotation interface. Its count and feature payloads are not yet qualified.** All 55,460 retained metadata rows join to eight distinct accession-local individuals, with four source `SNP` and four `noSNP` donors. All are male PSC cases. The paper identifies these groups as homozygous carriers and non-carriers of **rs56258221**.

The eight filtered archives are publicly reachable by HEAD and total **325,315,169 bytes**. No archive body was downloaded. Approval is required before that next step. **E-MTAB-14103 remains blocked as a genetic replication:** its 12 liver samples have ages and sex but no deposited genotype field or verified participant bridge.

The new clinical check matters: the scRNA genotype groups differ in published IBD labels, treatment exposure and severity measures. This is an observational reanalysis of the original study, not an isogenic experiment, an independent validation cohort, or a test of the published translation mechanism. These limits do not delay the separate MacroMap or liver branches.

## Sources and acquisition boundary

The [Poch/Bahn et al. study](https://doi.org/10.1016/j.xcrm.2024.101620), [E-MTAB-14013](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14013), and [E-MTAB-14103](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14103) were selected in the preceding [frontier review](psc-analysis-frontier-review.md). Its cached full text, BioStudies JSON and request receipts were reused and rehashed: **252,794 bytes**, not charged again. The qualification plan was saved before new metadata acquisition and any new expression effects.

New acquisition used **9,203,580 response-body bytes of the 15,000,000-byte authorization**: 20 logical requests, 34 HTTP hops, eight GETs and 12 HEADs. The body budget includes a 21,364-byte HTML challenge response, which is excluded as supplement evidence. A first ordinary PMC attachment HEAD returned 404. The exact published PMC attachment link returned the challenge page; it was not bypassed. The ordinary public publisher attachment succeeded. The Europe PMC supplementary-files endpoint was HEAD-only; its zero Content-Length is not evidence of an empty archive.

Successful bodies include both studies' IDF/SDRF files, `scrna_meta.tsv` (5,914,294 bytes), an article page used only to locate its supplement, and [Document S1](https://ars.els-cdn.com/content/image/1-s2.0-S2666379124003240-mmc1.pdf) (2,928,108 bytes; SHA-256 `5d37298b5a7efa1b0708f25fd2842c6eb1b0b129f4abfe3467cfb466de7d11ba`). Only Table S2, PDF page 8, was interpreted from the supplement. Its clinical values were not copied into public donor-level rows.

All new URLs, dates, HTTP statuses, redirects, headers, bytes and hashes are in ignored `data/raw/psc-bach2-qualification/request-ledger.jsonl`. Reused pins and the original scope are in `plan.json`. The [aggregate JSON](../data/derived/psc-bach2-input-qualification.json) carries source pins and qualification statuses. Raw bodies, full joins, clinical rows and independent checks remain under ignored `data/raw/` and `work/` paths. HEAD ETags and empty-response hashes are **not** remote archive content hashes.

## What the scRNA interface establishes

- The SDRF has three repeated `Derived Array Data File` columns. All were retained: one raw-feature archive, one filtered-feature archive and the shared cell metadata. These are explicit links, not inferred sample order.
- The metadata `sample_name` matches the linked filtered-archive basename after removal of `_filtered_feature_bc_matrix.tar.gz`. That supplies a unique source sample, individual and genotype. Cell `condition` and sex agree with the SDRF in every row.
- All required sample, genotype, sex, state and barcode fields are nonmissing. All 55,460 prefixed barcodes and all sample-plus-stripped-barcode pairs are unique. Removing only the exact `sample_name + '_'` prefix gives a candidate per-archive barcode. **1,769 bare barcode strings occur in multiple samples**; a global barcode-only join would be wrong.
- Matrix barcode membership, the pre-QC denominator, feature dimensions and count validity remain untested. Counts retained in annotation alone cannot measure gene detection or RNA coverage.
- The metadata has 16 published major-state labels. Paper `C15: MAIT` maps to Seurat `16`, not `15`. Metadata `C9: CD4+ TC` differs from Methods `C9: NKT`; neither was relabeled. Other source wording differences also remain in the original records.
- The four later naïve-CD4 activation subclusters are **not supplied** in `scrna_meta.tsv`. `C0: CD4+ TN RTE` and `C1: CD4+ TN mature` must not be renamed as those subclusters. UMAP coordinates were not interpreted and no clustering was repeated.

The following are retained-annotation eligibility counts, **not genotype composition effects or expression contrasts**. Both thresholds are inclusive. Each slash-separated entry is eligible `SNP`/`noSNP` donors; the maximum is 4/4.

| Source paper label | Metadata cells | ≥20 cells/donor | ≥50 cells/donor |
|---|---:|---:|---:|
| C0: CD4+ TN RTE | 9,320 | 4/4 | 4/4 |
| C1: CD4+ TN mature | 7,339 | 4/4 | 4/4 |
| C2: CD4+ TCM | 6,203 | 4/4 | 4/4 |
| C3: CD8+ TN | 5,774 | 4/4 | 4/4 |
| C4: CD8+ TC | 4,286 | 4/4 | 4/4 |
| C5: CD8+ TC term. diff. | 4,089 | 4/4 | 4/3 |
| C6: CD8+ TEM TC1 | 3,615 | 4/4 | 4/4 |
| C7: CD4+ TEM TH2 | 2,981 | 4/4 | 4/4 |
| C8: CD8+ TEM NK-like | 2,512 | 4/4 | 4/4 |
| C9: CD4+ TC | 2,146 | 3/2 | 3/2 |
| C10: CD4+ TEM TH1 | 1,807 | 4/4 | 4/4 |
| C11: CD4+ TEM TH17 | 1,802 | 4/4 | 4/4 |
| C12: CD4+ TREG naive | 1,358 | 4/4 | 4/4 |
| C13: γδ T cells | 1,068 | 4/4 | 4/3 |
| C14: CD4+ TREG effector | 814 | 4/4 | 4/2 |
| C15: MAIT | 346 | 4/3 | 3/1 |

Both naïve compartments and the source TH1, TH2 and TH17 compartments retain all eight donors at both thresholds. This permits a later fixed-compartment feasibility test, not a claim that BACH2 or any program is sufficiently measured. No endpoint or program was selected from expression here.

## Assay, reference, and retention limits

The source material is FACS-sorted peripheral **CD3+ T cells from cryopreserved PBMCs**, not whole blood, whole PBMCs, liver cells or a CD4-only stimulation experiment. The source specifies Chromium Next GEM **5′ v2**, gene-expression/ADT and V(D)J library methods, and NovaSeq6000. Methods state “TotalSeq A” while the key resources list TotalSeq-C products. Actual deposited RNA/ADT/TCR contents are unknown; the article says ADT processing occurred “where … available.” Do not assume every sample has every modality or that ADTs include BACH2 protein.

The scRNA alignment reference is literally **GRCh38 human reference genome (release 93)** with **Cellranger 3.0.2**. This does not identify the rs56258221 coordinate, genomic strand, REF/ALT or risk nucleotide. The source assay code `C__88670967_10` and VIC/FAM detection do not supply forward-strand allele coding. Those fields remain null. Neither rs56258221 nor these donor labels was equated with the project's other BACH2 anchors.

Original processing drops genes seen in fewer than 1% of a sample's cells, uses sample-dependent UMI/gene/mitochondrial gates and adjusted Scrublet thresholds, then SCTransform. Reported ranges are minimum UMI 800–2,000; minimum genes 450–1,000; maximum genes 2,500–3,000; mitochondrial proportion 5–7.5%. Exact per-sample settings and executed code are not deposited in the inspected sources. The SDRF attaches the transformation protocol to both raw/filtered filenames. **The filenames alone do not prove untouched integer UMIs or a complete pre-gene-filter feature universe.**

Methods report **56,209 cells** (31,014/25,195), while Figure 2 and the deposited annotation have **55,460** (30,609/24,851). This 749-cell difference is preserved, not silently repaired or attributed to a guessed exclusion. A matrix retention audit is still needed. PolyA 5′ RNA does not directly measure mature miR4464 or establish BACH2 protein translation.

## Clinical confounding: newly recovered Table S2

Table S2 explicitly refers to Figure 2A–J, the eight-person CITE-seq comparison. It supplies clinical records grouped by published SNP status. These aggregate facts replace the earlier *uninspected-supplement* missingness boundary; they do not establish a patient-to-archive clinical key.

| Source field | Homozygous, n=4 | Non-carrier, n=4 |
|---|---|---|
| IBD labels | `PSC-assoc. colitis` 3; `/` 1 | `Crohn's disease` 1; `/` 3 |
| `immunosuppr. therapy` | `Vedolizumab` 1; `/` 3 | `/` 4 |
| Birth-year range | 1969–1987 | 1966–1979 |
| Elastography range, kPa | `5,3`–`20` | `4,7`–`7,8` |
| MELD range | `6,0`–`11,8` | `6,0`–`6,0` |

Decimal-comma source strings remain unchanged; conversion to decimal points was used only to order range endpoints. Birth years are not ages at sampling. Gender cells literally say `cis`; male sex is supported separately by the main caption and SDRF. Table S2 does not define `/` as absent, untreated, unavailable or not applicable. It must not become a negative clinical value without evidence.

The table's patient numbering and genotype blocks agree with the SDRF's numbering pattern, but no explicit patient-number → Source Name/individual/archive key is printed. That is a **candidate correspondence**, not a verified row-level covariate join. No clinical rows were merged into expression inputs. Aggregate genotype-group imbalance is already clear without that merge. Exact collection ages, ancestry, run/batch, genotype-batch balance, complete medication timing/doses and clinical selection remain unresolved. Eight donors cannot supply a credible unrestricted adjustment for all these factors; more cells do not add independent clinical units.

## E-MTAB-14103: useful metadata, no genotype replication

The SDRF contains **12 unique accession-local individual labels**, **9 male/3 female**, and 12 nonmissing ages, range **18–59 years**. All are PSC liver samples. The shared `count_table.txt` has an HTTP-200 HEAD with the manifest's **2,681,853-byte** length. No body or count header was read, so its sample-column and feature interfaces remain untested. There are **zero genotype fields**, and numeric individual labels were not transferred from the blood accession.

The source describes 2011–2016 liver biopsies, TruSeq stranded mRNA, HiSeq4000, FastQC 0.11.5, Trimmomatic 0.36 and STAR 2.7.3a against the literal `ensemble 87` reference/annotation. This is not the blood alignment release. Preserve three source limits: age 18 versus the paper's >18 inclusion wording; SDRF library strand `not applicable` versus stranded-mRNA Methods; and a Methods sequencing-summary scope of 47 samples versus the 12 deposited liver samples. The paper says de-identified bulk RNA is available on request, whereas the current public manifest/HEAD exposes this processed file. Neither statement was erased.

Participant overlap with blood or other assays is unverified. This liver file is not an independent genetic replication or a direct CD4 translation assay. Reopen that question only after a public, explicit liver sample–genotype key and participant-overlap information are available.

## Validation and reopening rule

The native `.venv` parser has **30 passing synthetic tests**, covering duplicate SDRF fields, exact sample/barcode joins, collisions, conflicting genotypes, missing states, source-label crosswalks, clinical slash/decimal preservation, header-only acquisition, redirect constraints, cumulative byte caps and invalid supplement format. An independent implementation, which did not import the qualification module, passed **104 source-only checks**. All **69 aggregate comparisons** with its results agree, including the complete state/threshold/group table. Separate page-8 extraction and visual review agree with Table S2. The aggregate output replays byte for byte offline.

The full repository command was also run in `.venv`: **581 tests executed, seven errors** from missing older ETS2/GSE84161 caches or the missing historical Rscript runtime. No BACH2 test failed. The exact errors are preserved in `work/psc-bach2-qualification/full-tests-final.log`; no dependencies were installed or unrelated files repaired. Concurrent branches were adding tests, so this is the recorded run's count, not a frozen repository-wide total. Catalog import logic was not changed.

Reproduce this qualification from pinned caches:

```text
.venv/bin/python scripts/qualify_psc_bach2_inputs.py qualify
.venv/bin/python -m unittest discover -s tests -p test_psc_bach2_inputs.py -v
```

The PDF text derivative and its `pdftotext 26.08.0` command/hash are pinned in `supplement-text-provenance.json`. Full source joins and verification records are kept in ignored `work/psc-bach2-qualification/`.

**Next authorization request:** at most **325,315,169 additional bytes** for the eight filtered count archives, initially stopping after the smallest 26,893,378-byte archive if its full feature/assay/count interface is unsuitable. This would be payload qualification only: archive integrity, complete feature identifiers/types, count quantity, sparse structure and retained-barcode membership. No archive is approved by this report. An expression question, population definitions, feature mapping and all endpoints must be frozen separately before effects. The other eight listed raw-feature archives total 442,738,767 bytes and were neither requested nor acquired; raw reads are not publicly deposited for privacy.

No normalization, program score, fit, donor gene contrast, coordinate/allele inference, model call, private access, researcher contact, enrollment, installation or Git operation was performed. The work was not blinded to the already published narrative: a source-term search exposed neighboring published result wording, but it did not determine state selection or any new contrast. No deposited expression values were inspected. An association among the eight PSC donors would remain observational; RNA alone cannot validate the published miRNA–protein mechanism.

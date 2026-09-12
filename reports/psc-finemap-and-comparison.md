# PSC fine-mapping archive audit and prospective comparison

Completed September 12, 2026. Source retrieval began September 11 in America/New_York. **The complete released fine-mapping tables are imported, five candidate substitutions are reference-verified, and the comparison rules are fixed. Comparator matching and model inference have not been performed.**

The archive supports a limited comparison of AlphaGenome with statistical prioritization. It does **not** provide enough consistent run metadata to reproduce the published signal-specific credible sets. We retain the discrepancies instead of treating an inferred set as the authors' result.

## What the archive contains

We downloaded the 16,716,205-byte fine-mapping archive from [Goode et al.'s Zenodo release](https://zenodo.org/records/11143561), verified its publisher MD5 and pinned its SHA-256. The separate 7.64 GB expression-QTL archive was not downloaded.

| Content | Count |
|---|---:|
| GWAS regions | 18 |
| Molecular-QTL datasets | 8 |
| GWAS variant rows | 46,248 |
| Molecular-QTL variant rows | 24,835 |
| Nonempty configuration tables | 25 of 26 |
| Explicit credible-set files | 0 |
| Original association-input or LD-matrix files | 0 |

All 71,083 rows are preserved in the reproducible [compressed variant table](../data/derived/psc-finemap-variants.csv.gz). A row belongs to a dataset; repeated variants across datasets are not additional independent observations. Original identifiers, row/index values, four-decimal probabilities and log Bayes factors remain intact. A separate column contains the configuration-summed marginal probability. Missing configuration values stay blank, and source `-inf` log Bayes factors remain `-inf`.

The [82-member inventory](psc-finemap-archive-members.csv), [26-dataset summary](psc-finemap-datasets.csv) and [detailed audit](psc-finemap-audit.json) make the import inspectable. GPR35, HDAC7 and SIK2 are additional archive regions outside the 15 regions in Table 1. Their presence alone does not establish why they were omitted from that table.

## Source consistency and reconstruction limits

For **all 25 nonempty configuration files**, summing configuration probabilities for each included variant reproduces the matching SNP-table marginal probability within its four-decimal rounding precision. These two file types agree internally. This summation follows the definition in the [FINEMAP authors' documentation](http://www.christianbenner.com/); we do not replace the original published values with the reconstructed numbers.

The accompanying logs do not consistently describe those files:

- All 26 logs identify **FINEMAP v1.1**; the [paper's Methods](https://www.nature.com/articles/s41467-024-53602-w) identify v1.3.
- Every GWAS log specifies a maximum of one causal variant, but **16 of 18 GWAS configuration files contain configurations with two or more variants**.
- The posterior distribution over the number of causal variants differs between logs and configuration files in **23 of 25 comparable datasets**, including all seven molecular-QTL datasets with configurations.
- One PRKD2-region molecular-QTL configuration file is empty: `19_46705707_47705707_gene_nor_combat_mono.config`. The distinct `_2CV` file has data and is preserved separately.
- GWAS logs record 11,386 in the field labeled “Number of individuals in GWAS.” This cannot be silently equated with the figure's full case/control sample scope. The archived logs also append the unit `SNPs` to that individuals field.

These are source-version and provenance questions, not evidence that the underlying study is invalid. The archive cannot establish which run generated each published table or whether the saved configuration probabilities cover every model considered.

Marginal probabilities can sum above one when multiple variants may be causal. Dividing them by their sum and accumulating to 95% would not reconstruct a per-signal credible set. The archive has no explicit signal-membership files or original `.z`/`.ld` inputs to resolve this.

Only UBASH3A and SIK2 have exhaustive, singleton-only configurations. We reconstructed separate 95% and 99% sets, conditional on those archived one-causal-variant models, retaining all exact boundary ties:

| Region | Reconstructed 95% size | Included probability | Reconstructed 99% size |
|---|---:|---:|---:|
| UBASH3A | 7 | 0.9615628537 | 103 |
| SIK2 | 790 | 0.9500477350 | 1,803 |

The [membership table](../data/derived/psc-finemap-singleton-reconstructions.csv) is explicitly a reconstruction. Table 1 lists five UBASH3A variants; our calculation does not reproduce that number or establish the authors' exact set-construction rule. We do not relabel their five-member result as a 95% set without documentation.

The [Table 1 crosswalk](../data/derived/psc-finemap-table1-crosswalk.csv) also preserves unresolved aliases and probability differences. For example, CD28 signal 2 lists rs231799 with probability 0.17, while the archive's literal rs231799 record has 0.0004. UBASH3A has 0.62 in Table 1 and 0.6104 in the SNP file. Coordinate-only matches remain provisional. `BCL211`/`BCL2L11`, `IL2-IL21`/`IL21` and `ETS2`/`PSMG1` are recorded as region-label links, not assertions that the labels denote equivalent genes.

The full archive has no explicit REF/ALT or build columns. We have verified the shortlisted rsIDs below, **not all 71,083 rows**. Anonymous coordinates and array identifiers require additional reconciliation.

## Fixed candidate selection

The [protocol](../config/psc-controlled-comparison.json) selects previously unscored regions with one published signal and a published top probability from 0.1 to below 0.5. It orders them by published set size, then decreasing top probability. This yields **IL2RA, BACH2 and BCL2L11**; the original Table 1 label `BCL211` remains preserved. Set size is a selection field here, not proof of recovered membership. The [selection audit](../data/derived/psc-comparison-locus-selection.csv) includes every Table 1 region and its inclusion/exclusion reason.

Within each selected archive region, records with marginal probability at least 0.1 are assessed. The rule retains at most two unambiguous, noncoding, biallelic SNVs in decreasing probability order. Nine rsIDs were checked against Ensembl release 116 on both genome builds. Six SNVs passed identical-allele and reference-sequence checks against both Ensembl builds and AlphaGenome's own GRCh38.p13 FASTA; one remains a reserve.

| Region | Nominated variant | Archived probability | GRCh38 position, 1-based | REF → ALT |
|---|---|---:|---|---|
| IL2RA | rs7923054 | 0.1084 | chr10:6169124 | T → A |
| BACH2 | rs7750271 | 0.2009 | chr6:90326506 | A → G |
| BACH2 | rs72928086 | 0.1964 | chr6:90320043 | C → A |
| BCL2L11 | rs72837826 | 0.1771 | chr2:111175424 | G → T |
| BCL2L11 | rs13405741 | 0.1695 | chr2:111155479 | T → C |

These are region labels and candidate substitutions, not established causal gene assignments. Disease-risk direction remains unspecified.

The [complete 12-record candidate audit](../data/derived/psc-comparison-candidates.csv) shows the six exclusions and reserve. In particular, IL2RA's highest-probability record, rs4147359, is currently G/A/C/T in Ensembl; we do not guess which pair was analyzed in the fine-mapping study. rs68021656 is a variable-length sequence, and rs72837816 has three alleles. Two coordinate-only identifiers and one `kgp` identifier lack sufficient direct mapping. BACH2 rs7750559 passes reference checks but falls below the fixed two-anchor cap. These exclusions limit the region coverage of the planned analysis.

## What the controlled comparison will test

For each included anchor, select two **low-PIP comparison variants**, with archived probability no greater than 0.001. These are not experimentally established negative controls. There are 7,374 source records below that threshold across the three regions, including 6,117 literal rsIDs; none has yet passed the complete matching process.

The protocol fixes these rules before any new AlphaGenome scores:

- Use 1000 Genomes phase 3 EUR frequencies and genotypes. Match common-variant frequency within 0.05, mutation class, noncoding annotation class and nearest-TSS distance. Require at least 400 jointly called EUR donors for an LD calculation. Missing LD is not zero LD.
- Require comparator dosage-correlation r² ≤ 0.1 against all reference-resolved high-PIP SNVs, including the reserve, and other assigned comparators. The ambiguous/indel records are outside this screen; report that incomplete coverage. Do not reuse comparators or relax failed criteria.
- Use one identical 1,048,576-base window and one fixed gene universe per region. The [provisional windows](../data/derived/psc-comparison-windows.csv) are already defined. Pin the exact CD4 RNA-seq tracks and model gene annotation before inference. Missing outputs make a comparison unavailable.
- Use the fixed RNA-seq gene-mask log-fold-change scorer. For each variant, take the largest gene effect after calculating the median absolute score across the fixed tracks. Compare each anchor with the median of its two comparators, then report contrasts by region. Preserve all underlying signed scores and nearest-TSS comparisons.

The intended batch is **15 requests if all five anchors obtain two valid comparators**, within a protocol cap of 18. No prediction request has been made. The protocol, source manifest, complete variant table and candidate table are recorded by checksum in the [freeze record](../config/psc-controlled-comparison-lock.json).

This is a prospective comparison with **statistical prioritization**, not an independent accuracy benchmark. The archived probabilities are uncertain labels; region names are not causal-gene truth. The fixed `ALL_FOLDS` model has not held these genomic regions out of training, as described in the [official model-version documentation](https://www.alphagenomedocs.com/api/generated/alphagenome.models.dna_client.ModelVersion.html). Donor overlap remains unestablished. Three regions cannot support a broad generalization claim or a treatment conclusion.

## Reproduction and next action

The new importer and preparation workflow use Python's standard library. In the existing project environment:

```bash
python scripts/fetch_sources.py --manifest config/finemap-sources.json
python scripts/audit_finemap_archive.py
python scripts/prepare_controlled_comparison.py
python -m unittest discover -s tests -v
```

Add `--offline` to the download command to verify the cache. The 43 pinned sources include the full fine-mapping archive and the small reference/documentation files; raw downloads stay outside Git. The importer refuses unsafe archive members, invalid probabilities, duplicate identifiers and configuration references to absent SNPs. The new tests exercise multicause versus singleton reconstruction, missing data, log disagreements, source matching, variant selection and reference failures.

Verification completed: all **43 source hashes** checked, **36 tests passed**, and **11 generated outputs** reproduced with identical hashes on replay. Local document links resolve. The original Table 1 importer and previous prediction results were not changed.

**The next executable phase is comparator preparation:** obtain the required EUR frequency/LD and annotation data, assign the two comparison variants per anchor, and lock the final matching table and common gene/track universe before inference. The [preparation status](psc-comparison-preparation.json) explicitly remains `ready_for_inference: false`. Short [source-provenance questions](../docs/finemap-provenance-questions.md) are drafted and unsent; a reply would improve confidence but is not required for this limited exploratory comparison.

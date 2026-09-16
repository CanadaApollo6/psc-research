# PSC BACH2: single-archive count qualification

Date: 16 September 2026 UTC. **One authorized archive; no expression effects.** The four preceding [metadata qualification](psc-bach2-input-qualification.md) deliverables remain unchanged.

## Decision

The smallest filtered archive passes integrity, complete deposited feature-axis, sparse-value and fixed-retained-barcode checks. It contains a **mixed RNA/antibody matrix**, not an RNA-only matrix. All **4,278** matching retained metadata cells are present among **4,703** source barcodes; **425** source barcodes are outside that retained metadata.

The matrix explicitly declares **`cellranger-6.1.1`**, whereas the paper and IDF state **Cellranger 3.0.2**. This supports a Cell Ranger MEX-form count interface beyond filenames and integer values. It does **not** certify the historical export step, reference bundle, or untouched original RNA UMI versus a later SCT-derived export. The fixed metadata has no RNA/SCT count margins for a direct comparison. No count-likelihood or expression-analysis gate is declared passed.

**The one-archive phase is complete.** The other seven archives were not acquired. A larger payload-qualification phase could establish all-eight feature/assay and retained-cell interfaces. It could then support a separately frozen fixed-compartment question, conditional on its RNA-quantity definition. It would not by itself recover missing historical RNA/SCT margins or settle the version/export discrepancy if the same sparse layout is repeated.

## Authorized acquisition and source pins

The [public source](https://www.ebi.ac.uk/biostudies/files/E-MTAB-14013/sample05_filtered_feature_bc_matrix.tar.gz) is `sample05_filtered_feature_bc_matrix.tar.gz`:

- **26,893,378 bytes**, exactly the approved inventory size.
- SHA-256 **`5c57cc1c8eb48ce0208062b92fcd056b39e0f90e8eb6cb712ccd0d76f95dda29`**.
- Retrieved 16 September 2026, 02:00:40–02:03:20 UTC, through the ordinary BioStudies → EBI FTP HTTPS redirect.
- The separate allowance was at most 1,000,000 new metadata bytes. **Zero metadata body bytes** were read. Two bounded ordinary official Cell Ranger documentation requests returned 403 and 429; their headers/failure receipts are retained, but neither is treated as format-definition evidence. No challenge or restricted route was bypassed.

This archive allowance is additional to the earlier metadata budget. No other filtered archive, raw-feature archive, read archive or liver count body was requested. The initial plan, acquisition receipt and original compressed body are under ignored `data/raw/psc-bach2-qualification/counts/`. The [aggregate JSON](../data/derived/psc-bach2-count-qualification.json) provides source pins and qualification statuses. Full source features, matrices, barcodes and technical margins remain ignored.

## Integrity and actual schema

The outer gzip is consumed through full EOF, including CRC/ISIZE validation. It expands to a **26,900,480-byte TAR**, SHA-256 `96de2de1d7600902db0bb92c70c76a762335108dc9c44168b37d61f04f5f895b`. The TAR has one directory and three regular files, all inside `sample05_filtered_feature_bc_matrix/`:

| Source member | Compressed bytes | Decoded bytes |
|---|---:|---:|
| `barcodes.tsv.gz` | 24,126 | 89,357 |
| `features.tsv.gz` | 333,835 | 1,480,573 |
| `matrix.mtx.gz` | 26,530,656 | 85,682,856 |

All three nested gzip streams also pass full CRC/ISIZE/EOF checks. Every TAR header, path, member type and size is checked. There are no links, sparse/special files, path traversal or normalized path aliases. A complete zero TAR trailer is required. The independent parser also checks all 512-byte header checksums, payload padding and terminal zero blocks. Source member names are **never used as extraction destinations**; controlled work-file names and exclusive creation prevent path or source overwrite. Compressed and decoded member hashes are pinned. Byte, member-count, dimension and coordinate caps are enforced.

The matrix is exactly:

```text
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-6.1.1", "format_version": 2}
```

It has **36,639 feature rows × 4,703 barcode columns** and **6,835,238 declared and observed sparse coordinates**. Every coordinate is one-based and within its axis. Every stored value is finite, positive and exactly integral, with minimum 1 and maximum 21,236. There are no duplicate coordinates, explicit zero entries or fractional values. Coordinates are strictly column-major. Numeric tokens are checked as integers or exact decimal/rational values **before any float conversion**; nothing is rounded, summed to repair duplicates, or relabeled as an integer count.

## Complete deposited feature axis and assay separation

All 36,639 original feature IDs are unique. There are 36,629 unique source names, with ten duplicated name values; no name-based collapse or alias mapping was applied. The three source columns are feature ID, feature name and feature type. All RNA IDs have unversioned `ENSG` syntax. This is an identifier-format check, not independent gene, build or reference-release reconciliation.

| Source feature type | Feature rows | Positive coordinates | Exact total source units | All-zero feature rows |
|---|---:|---:|---:|---:|
| `Gene Expression` | 36,601 | 6,714,934 | 20,228,805 | 18,253 |
| `Antibody Capture` | 38 | 120,304 | 13,473,321 | 2 |
| Combined | 36,639 | 6,835,238 | 33,702,126 | 18,255 |

These are structural/technical aggregates, not gene effects. **33,702,126 is not an RNA UMI total.** Any later RNA aggregation must use the original `Gene Expression` rows and exclude `Antibody Capture` rows from both RNA totals and features. RNA and antibody column margins were calculated separately only for technical verification and kept ignored. Every barcode has positive entries in both feature types.

The full deposited axis was read, including its all-zero rows. No gene was imputed or dropped under the paper's historical 1% filter. An intact deposited axis does not establish the complete pre-filter experimental reference universe. There is no reference bundle manifest, HDF5/AnnData/Seurat layer object, SCT matrix or normalization record in this three-file payload. The feature count is not used to guess a different genome release.

## Fixed retained-cell interface

The published SDRF uniquely links the archive basename to the source sample. Exactly that basename supplies `sample05`; exactly `sample05_` is removed from matching retained metadata barcodes. This reproduces the preceding fixed join without a barcode-suffix repair or global barcode-only match.

- Source matrix: **4,703 unique barcodes**.
- Fixed retained metadata for this source sample: **4,278 unique barcodes**.
- Retained barcodes matched: **4,278**; missing from matrix: **0**.
- Matrix barcodes outside fixed retained metadata: **425**.

Those 425 barcodes are not silently discarded from the source or assigned new states. Their absence from the retained annotation is not a proven historical exclusion reason. This single sample cannot reconcile the paper's 56,209-cell versus deposited 55,460-cell discrepancy. No state selection or genotype/state gene contrast was performed.

## Quantity evidence and the remaining qualification limit

The exporter comment is direct source evidence for a declared **Cell Ranger MEX format**, not merely an inferred filename convention. Mixed `Gene Expression` and `Antibody Capture` axes, a full sparse integer matrix and complete barcode membership are consistent with that declaration. No normalized or SCT layer label is present.

However, the prior source Methods and IDF alignment protocol `P-MTAB-142881` state **Cellranger 3.0.2** and **GRCh38 human reference genome (release 93)**. The matrix instead says **cellranger-6.1.1**. Both statements remain unchanged. The source processing protocol also describes sample-level gene/QC filtering, Scrublet and SCTransform; the SDRF associates transformation protocols with both raw and filtered archive links. The exact export command is not supplied. Integer-valued transformed counts can exist, so integers alone do not exclude every SCT-derived layer.

The fixed `scrna_meta.tsv` has no `nCount_RNA`, `nFeature_RNA`, `nCount_SCT`, `nFeature_SCT` or comparable RNA/SCT margin columns. Thus historical RNA-versus-SCT margin reproduction is **unavailable**, not a failed match. Current reconstructed margins are not substituted for missing source margins. The data are characterized as **RNA-feature source count-format entries**, without certifying their untouched original-UMI lineage or applying a historical filter.

## Verification and reproducibility

The native `.venv` suite has **30 passing synthetic tests** for this new parser. These cover gzip corruption/truncation, TAR checksums/traversal/aliases/member types/EOF, size caps, exact values beyond float precision, fractions that would round in float, duplicate or malformed coordinates, complete axes, exact sample-prefix membership, cache-to-source lineage, exclusive output creation, redirect restrictions and cumulative metadata budgets.

An independent stdlib-only implementation did not import the main parser. Its **seven synthetic tests pass**, and two full source passes agree. The main and independent implementations agree on:

- 24 structural/assay/join aggregates;
- **146,556** full feature-axis fields, including source row positions;
- **14,109** barcode-position and retention fields;
- **9,406** exact RNA/antibody column margins.

This is independent implementation, not blinded verification: expected aggregates were shared before the independent scan. The aggregate JSON replays byte for byte offline. All four preceding metadata deliverable hashes remain unchanged.

The native full repository suite was also run: **637 tests executed, seven errors** from missing historical ETS2/GSE84161 caches or the historical Rscript runtime; no new BACH2 test failed. Its log is `work/psc-bach2-count-qualification/full-tests.log`. This count records that run while other branches were adding tests. No unrelated files were repaired and no packages were installed. Catalog import logic was unchanged.

Two implementation-only repairs are retained in logs: `Counter` dictionary-key construction in the member-role check, and quoting in a synthetic MatrixMarket header. Both were corrected before the successful final scan/tests; no source value, permitted archive, population or analysis rule changed.

```text
.venv/bin/python scripts/qualify_psc_bach2_counts.py qualify
.venv/bin/python -m unittest discover -s tests -p test_psc_bach2_counts.py -v
.venv/bin/python work/psc-bach2-count-qualification/compare_independent.py
```

The original compressed archive and receipts suffice to regenerate disposable work files. Existing work caches are rehashed and linked back to the original TAR payload; unverified existing paths are not overwritten.

## Next boundary

**Approval is required before the remaining seven filtered archives: 298,421,791 additional bytes.** If approved, their purpose should be full donor-interface, assay-type, exporter-comment and fixed-retained-barcode qualification only. Preserve source RNA/ADT separation and the unresolved version/reference/export record. Do not treat a repeated file schema as a missing historical layer-margin certificate. A later fixed-compartment expression question still needs its own explicit quantity definition and freeze before effects.

No other archive download, normalization, program score, model fit, gene contrast, expression-based exclusion, state selection, genotype inference, researcher contact, private access, installation or Git operation occurred. The earlier metadata and interpretation record is unchanged.

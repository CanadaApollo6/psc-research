# Source questions after the variant-coverage audit

Prepared September 12, 2026. **Unsent.** These requests concern public aggregate statistics and analysis metadata. The [completed report](../reports/variant-coverage-repair.md) and machine-readable tables make each question reproducible.

For the PSC study authors or GWAS archive maintainers:

1. Is an association record available for the exact four-base deletion **GRCh38 chr10:6,161,637 TAGAG>T**, equivalent to **GRCh37 chr10:6,203,600 TAGAG>T**? Published RNA annotations include rs149556258/rs71390193. Please retain REF/ALT and effect-allele definitions, genome build, imputation quality, allele frequency, per-variant case/control counts and original statistics. We found no equivalent deletion among the 7,639 original regional GWAS rows, including alternative repeat placements; it carries 47.55% of the strongest BLUEPRINT PFKFB3 RNA component.
2. Could the full PFKFB3-region association inputs and matching signed in-sample LD or score covariance be shared, including the two exact recovered deletions at GRCh37 6,203,538 TAC>T and 6,212,648 AGTATGGAGGTGGGGACAGG>A? Their original GWAS rows are rs138248884 and rs140499372. The historical file's allele_1 encodes the longer reference allele as risk in both records. Reference normalization and frequencies support that orientation, but an explicit original REF/ALT mapping would remove the remaining provenance ambiguity.
3. Which genotyping-platform samples entered each original fine-mapping run? The [earlier unsent questions](finemap-provenance-questions.md) explain why the Omni total of 11,386 and the archived run logs do not by themselves establish the exact input covariance.

For the eQTL Catalogue maintainers:

1. What variant-level QC or processing decision explains the absence of these four GRCh38 records across all genes in the nine selected Schmiedel_2018 CD4 datasets?

   | Variant | Exact REF>ALT | Position |
   |---|---|---:|
   | rs13405741 | T>C | chr2:111155479 |
   | rs72837816 | A>C | chr2:111159291 |
   | rs72836352 | C>T | chr2:111137297 |
   | rs144569746 | C>T | chr2:111150990 |

   The datasets are QTD000439, QTD000444, QTD000449, QTD000454, QTD000459, QTD000464, QTD000469, QTD000479 and QTD000484. BLUEPRINT QTD000031 contains all four exact BCL2L11 associations. The original DICE unfiltered downloads exist, but we have not extracted their target rows. Is there a compact aggregate extract or exclusion ledger with the actual imputation/filter metrics and run version? General profile defaults cannot establish the reason for a particular missing record.
2. For OneK1K_reannotated QTD000891, both the r8 dataset metadata and nominal AN/2 indicate **997**, while the study page lists **981 donors**. What is being counted, and what is the verified number of independent participants used for this association dataset? QTD000889 and QTD000892 similarly report 528 and 653. We have preserved all labels and used no exact-N-dependent new posterior.
3. The earlier Treg naive/memory mismatch for QTD000464/QTD000469 remains unresolved. Which cell-state labels should be used, and is there a corrected, versioned sample mapping?

Responses would be retained as new source evidence alongside the current files. None of these requests has been sent.

# Questions arising from the completed mechanism audit

Updated 2026-09-11 after the [direct-junction follow-up](../reports/ubash3a-junction-followup.md). These are draft questions for review, not sent correspondence. A [short, ready-to-review email](ubash3a-data-request-draft.md) and [six-row endpoint specification](../data/derived/ubash3a-data-request-endpoints.csv) now make the immediate request concrete. Public professional contacts remain in [contacts.md](contacts.md).

## eQTL Catalogue / Lepik cohort — most immediate data question

1. In QTD000377, phenotype `21:42434983:42437488:clu_31404_-` is assigned to UBASH3A with `strand=1`. The cluster suffix and annotated gene strand disagree for 115,561 of 116,015 unique single-gene-assigned phenotypes in the linked metadata. Is this an expected library-orientation convention, an old processing issue, or an annotation problem? Which archived workflow parameters or corrected release should we use? We have preserved both original labels.
2. Can we access the unfiltered rs1893592-C coefficients and uncertainty for the canonical `21:42434954:42437488` and alternative `21:42434983:42437488` junctions, with their cluster denominator definitions? The selected alternative row has beta 0.523483, SE 0.063608 and nominal P 1.98464e-15. We need both endpoints rather than a single selected tag trait.
3. Are de-identified per-genotype junction summaries, overhang/mapping-quality summaries or a public donor-level release available for a direct check? In DICE QTD000483/QTD000488 we found no UBASH3A-assigned phenotype in the published metadata; is fuller quantification available for these coordinates?

New evidence for question 1: both reference boundaries and Snaptron junctions are positive-strand GT–AG. In the published workflow, the reverse-stranded branch selects RF for HISAT2 but `-s 2` for regtools; the declared regtools 0.6.0 recipe defines RF as `-s 1`. This can explain reversal if those parameters were used, but actual historical run settings remain unverified. For questions 2–3, our GEUVADIS reanalysis has only 18 adequately covered European donors and no retained CC donors; it cannot supply a reliable effect estimate. A 2023 maintainer discussion confirms that full splicing files were too large to publish routinely and describes arranging a targeted transfer.

## PSC genetics / molecular-QTL investigators

1. For **rs313839**, can the conflicting original source labels be clarified? The GCST004030 GWAS row labels C as risk (OR 1.322, control frequency 0.84), whereas the thesis calls G risk and Supplementary Data 2 labels A. Reprocessed BLUEPRINT and separate DICE monocytes both show higher PRKD2 RNA per G allele (betas 0.908576 and 0.921304). Does that match the exact original molecular-QTL coding? Is the supplement's A label associated with the original lead rs60652743 rather than rs313839?
2. Does the **rs313829** label in the PRKD2 regional plot mean rs313839? Current database mappings put rs313829 on chromosome 7. Please distinguish the locus-level colocalisation probability from the individual variant's fine-mapping probability.
3. Are allele-coded regional summary statistics or full credible-set memberships publicly released for this locus? We can continue independently once the source fields are clear; individual-level genotypes are not required for this immediate check.

## UBASH3A expression / splicing investigators

1. Can the higher exon coverage reported by Newman et al. and lower total-transcript qPCR reported by Ge and Concannon be aligned to the same transcript segments, stimulation state and background haplotypes?
2. Which RNA quantity best distinguishes the **+29-nt donor** event from other retained-intron transcripts, including the transcript assayed in the 2018 experiment? We should not assume all “intron retention” measurements report the same isoform.
3. Which direct donor-level data could distinguish the +29-nt event, exon skipping and other retained-intron forms? The model predicts increased +29 use; Lepik blood gives provisional boundary/direction agreement with unresolved strand labels; the small 2020 lymphoblastoid experiment reports a different direction. DICE total gene RNA increases with C in both resting and activated CD4, contrary to the model's total-RNA prediction. These different endpoints should be tested separately.

No question assumes a therapeutic consequence. The purpose is to resolve reproducible molecular comparisons and identify an informative independent analysis.

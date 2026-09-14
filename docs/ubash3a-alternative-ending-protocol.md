# Conditional reconstruction of the TRAILS alternative ending

This is the second computational direction, performed while the raw-read audit runs. The fixed input is the previously identified TRAILS transcript `e09f04f5-c617-4bb3-b69c-519ccb16d0f1` (T29), not a new structure selected from the current read results.

Reconstruct all eleven exons from the pinned GRCh38 reference. Translate from the shared UBASH3A reference start at chr21:42403946 under the standard genetic code, ending at the first in-frame stop. This start is an explicit assumption for T29; it is not a measured initiation site. Retain both hypothetical rs1893592 A and C sequences without assigning the TRAILS donor's genotype. Do not invent a replacement start if the shared ATG is missing.

Report complete reference RNA, conditional ORF/protein, genomic stop coordinates, exon-junction position features using the preceding >55-nt primary and >50-nt sensitivity conventions, and codon-based correspondence to the pinned UniProt protein. Compare each candidate with all twenty previous normal/+29 scenarios. Sequence/domain correspondence does not establish folding or activity; a terminal-exon stop does not establish RNA survival.

Generate tryptic candidate peptides with cleavage after K/R except before P, zero or one missed cleavage, and lengths 7–35 residues. Preserve peptides shared with prior products as such, and label any sequence absent from the previous twenty products only as a **candidate for a broader specificity search**. Neither uniqueness across the human proteome nor detection in mass spectra is established here. Treat leucine/isoleucine as equivalent when checking this limited mass-spectrometry specificity.

An independent check reconstructs RNA through a base-by-base genomic coordinate map, translates with Biopython, checks the first stop and RNA/allele difference, and compares results with the primary literal-codon implementation. Preserve all prior outputs and pin the new inputs and calculations.

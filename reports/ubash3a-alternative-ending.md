# A distinct conditional product from the TRAILS alternative ending

The fixed TRAILS T29 model reconstructs to **4,392 RNA nucleotides** and, if translation starts at the shared UBASH3A reference ATG, a **482-amino-acid product**. Its first in-frame stop lies in its final exon. This gives a specific sequence to test for, but does not show that the protein is produced or that the RNA survives degradation.

This completes the conditional sequence reconstruction in computational direction 2. It uses the previously identified catalog structure. The subsequent [raw-read audit](ubash3a-raw-read-audit.md) did not confirm its exact +29 or distinguishing downstream junction in the two selected CD4 libraries; the sequence remains hypothetical.

## What was reconstructed

The source is [TRAILS](https://www.nature.com/articles/s41467-024-48615-4) model `e09f04f5-c617-4bb3-b69c-519ccb16d0f1`, previously isolated in the [annotation audit](ubash3a-hypothesis-test.md). Its eleven exons include the exact +29 boundary, chr21:42434983→42437488, followed by 42437580→42437855. The terminal exon ends at 42440816. Coordinates use GRCh38, positive-strand, one-based inclusive conventions.

The [fixed plan](../config/ubash3a-alternative-ending-plan.json) pins the source exon table, full-gene reference, UniProt record and previous products. Translation starts conditionally at chr21:42403946, RNA nucleotide 30. This is a shared reference start assumption, not an observed T29 initiation site. Both hypothetical rs1893592 alleles are retained; no donor genotype is assigned.

| Feature | Hypothetical A | Hypothetical C |
|---|---|---|
| Reconstructed RNA length | 4,392 nt | 4,392 nt |
| Conditional protein length | 482 aa | 482 aa |
| First stop | TGA | TGA |
| Stop, genomic coordinates | 42437900–42437902 | 42437900–42437902 |
| Stop, RNA coordinates | 1476–1478 | 1476–1478 |
| Amino acid 428 | M | L |
| Stop in final exon | Yes | Yes |

The A/C RNAs differ at the single specified retained base; their proteins differ only at amino acid 428. The terminal exon has 2,962 nt. The stop ends 48 nt after the final exon junction, leaving no downstream junction in this catalog model. Consequently, the earlier >55-nt and >50-nt downstream-junction flags are both false. Those limited positional rules do not establish NMD escape or measured RNA half-life.

## How it differs from the previous products

The appropriate comparison for separating the alternative ending is the upstream-matched U model, ENST00000291535.12. Relative to its normal 623-aa product, T29 shares the first 427 amino acids and has a 55-aa replacement tail. Relative to its previously reconstructed 498-aa +29 product, T29 shares the first 467 amino acids and has a distinct 15-aa ending, `LGIFKIHLESGIRKP`.

T29 also differs from the canonical model in upstream exon structure. It should not be described as simply truncating the canonical 661-aa sequence. All forty pairwise comparisons against the twenty previous scenarios are retained in [the comparison table](../data/derived/ubash3a-alternative-ending/previous-product-comparisons.csv).

Exact genomic-codon projection preserves all 46 reference UBA residues and all 66 reference SH3 residues. Only 70 of the 267 residues in the annotated phosphatase-like region retain an exact reference codon and amino acid. These are sequence correspondences, not measurements of folding, binding or activity. See [feature correspondence](../data/derived/ubash3a-alternative-ending/protein-features.csv).

## Candidate targets for the protein-data search

The declared tryptic digest permits zero or one missed cleavage and 7–35 residues. Four distinct sequences are absent from all twenty earlier modeled products, including an I/L-equivalent comparison:

| Candidate peptide | T29 residues | Missed cleavages |
|---|---|---|
| RPNSSWKLGIFK | 461–472 | 1 |
| LGIFKIHLESGIR | 468–480 | 1 |
| IHLESGIR | 473–480 | 0 |
| IHLESGIRKP | 473–482 | 1 |

These same four occur in both hypothetical allele products. They do not distinguish A from C. Their specificity across the human proteome and other possible isoforms has **not** been established, and no mass spectra have been searched. The [complete peptide dictionary](../data/derived/ubash3a-alternative-ending/tryptic-candidates.csv) also retains shared peptides and their matches instead of discarding them.

## Verification and reproduction

An independent implementation reconstructs RNA base by base from genomic coordinates and translates with Biopython. Both 4,392-base RNAs, both 482-aa proteins, the first stop and the single allele difference agree with the primary literal-codon calculation. Synthetic digest tests cover cleavage blocked by proline, missed cleavage and inclusive length limits. The full repository suite passes 212 tests at this stage.

Run `.venv/bin/python scripts/analyze_ubash3a_alternative_ending.py` from the research repository. The [audit](../data/derived/ubash3a-alternative-ending/audit.json) pins all nine primary outputs, including [RNA](../data/derived/ubash3a-alternative-ending/transcripts.fasta), [ORFs](../data/derived/ubash3a-alternative-ending/orfs-with-stop.fasta) and [proteins](../data/derived/ubash3a-alternative-ending/proteins.fasta). No new model prediction, experimental protein measurement, treatment inference or researcher contact was required.

# UBASH3A: conditional transcript and protein consequences

## Question and scope

If the previously identified alternative donor extends the upstream UBASH3A exon by 29 bases, what RNA and translated sequence would follow for each transcript carrying that junction? This is a deterministic sequence reconstruction conditional on that splice event. It does not estimate splice-form abundance, allele effects on splicing, protein abundance or clinical benefit.

The existing physical hypothesis is fixed: GRCh38 chromosome 21, positive strand; normal upstream exon ends at **42434954**, the next exon starts at **42437488**, and the alternative donor retains genomic bases **42434955-42434983**, inclusive. The rs1893592 A>C substitution at **42434957** lies inside the retained segment. Other alleles at this multiallelic site are outside this A/C analysis.

## Transcript selection and source freeze

Start from the preserved Ensembl release-116 expanded gene record for ENSG00000160185. Inventory every transcript, selecting every model with those exact consecutive exon boundaries, regardless of canonical or nonsense-mediated-decay annotation. The five previously recorded qualifying transcript versions are preserved in the plan. Other models receive an explicit exclusion reason, not an assumption of biological absence.

Retrieve current Ensembl release/assembly metadata and expanded annotation before sequence interpretation. Require matching versions, exon structures and translation boundaries for the selected models; report any drift and freeze an amendment before dependent effects. Pin the complete gene-region reference, each selected transcript's cDNA/CDS/protein responses and UniProt's reviewed human UBASH3A record. Record source URLs, retrieval times, response representation, bytes and SHA-256. Existing source files and earlier results are not overwritten.

## Sequence reconstruction and independent checks

Use the positive-strand reference gene sequence and ordered exon intervals to reconstruct mature cDNA. Require exact agreement with the Ensembl cDNA response. Map the annotated translation start/end through exon coordinates, independently validate the CDS response and reproduce the annotated protein. Identify explicitly whether the API's CDS includes the terminal stop; do not silently infer a stop or truncate a mismatch.

For each selected transcript construct four named scenarios: normal splicing with A, normal splicing with C, +29 retention with A, and +29 retention with C. The two normal RNAs must be identical because the variant is intronic and removed. The two retained segments must differ only at genomic position 42434957. Add precisely 29 bases at the selected exon junction, keeping all other splice junctions fixed. This synthetic construction is not asserted to be an observed full-length isoform.

Translate each scenario from its annotated start through the first in-frame stop, continuing into downstream transcript sequence if necessary. Preserve a no-stop result if one occurs. Report cDNA/CDS/peptide coordinates, the first changed amino acid, unchanged prefix, novel peptide tail, first stop codon and any removed reference residues. Distinguish a coding frameshift from UTR retention or a pre-existing stop. Use independent Biopython translation and a separately implemented sequence reconstruction to check every scenario; do not compare a function with itself.

## RNA-decay assessment

Report the first stop's distance to the final exon junction and the closest downstream junction on each altered mature RNA. Count the distance from the stop codon's final nucleotide to the exon boundary, and retain both stop start/end coordinates so the convention is auditable. Flag whether the stop ends **more than 55 nucleotides before the final junction**, and show a **50-nucleotide sensitivity flag**. Also report proximity to the annotated start (first 150 coding nucleotides), stop-bearing exon length and whether that exon exceeds 400 nucleotides. These are transparent sequence features relevant to nonsense-mediated decay, not a fitted probability or proof of degradation. Apply the same calculations to the normal transcript, including models already annotated for nonsense-mediated decay.

The primary literature establishes exceptions to simple positional rules; sequence context and expression can modify decay. No tumor-trained NMD efficiency is transferred numerically to UBASH3A in CD4 cells. [Lindeboom et al., 2016](https://doi.org/10.1038/ng.3664), [Lindeboom et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6858879/).

## Protein features and interpretation

Use the reviewed UniProt reference for domain, active-site and interaction-region annotations, preserving evidence types. Map features directly only where the Ensembl reference peptide exactly matches the UniProt sequence; do not assign canonical residue numbers to another isoform by position alone. Any additional projection must have an explicit sequence alignment and ambiguity accounting. Distinguish intact reference sequence, interrupted feature and absent sequence from demonstrated loss of function. A shortened product is a possible translation product, not evidence that it accumulates in cells.

Compare the reconstruction with the specific +29-base literature and with broader intron-retention studies without equating their RNA endpoints. State prior knowledge and disagreements. Consider structural modeling only after establishing whether a plausible surviving protein and a specific unresolved structural question justify it.

## Finished result

Deliver the source/parameter manifest, all-transcript inventory, scenario FASTA files, consequence and protein-feature tables, a readable sequence/domain figure, an evidence report and independent verification. Include negative and unavailable outcomes. Preserve previous analyses and unsent correspondence; no researcher contact, new model prediction or wet-lab intervention is part of this task. Run the repository test suite after adding source/analysis parsing and verify that existing outputs are unchanged except for dated navigation/evidence updates.

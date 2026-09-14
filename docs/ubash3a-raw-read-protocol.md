# Initial TRAILS raw-read audit

The first audit uses the complete naïve and memory CD4 libraries from the single TRAILS donor. Selection follows deposited cell labels and precedes target-read inspection. This is an exploratory examination of a previously observed catalog candidate. It is not independent validation of the catalog or a donor-level allele/decay comparison.

## Frozen question and inputs

Can uncorrected individual read alignments support the exact normal or +29 splice junction, and link it to the downstream boundary, canonical, or alternative-ending structures? All coordinates are GRCh38, chromosome 21, positive-strand exon endpoints, one-based inclusive. Normal is 42434954→42437488; +29 is 42434983→42437488. B, C, U and T29 retain the identities in the earlier assay-target table. Their full junction chains are pinned before execution.

Retrieve only DRR465268 (NCD4), DRR465269 (MemCD4), and Ensembl release 112's GRCh38 primary assembly. The approximately 7.12 GB compressed transfer fits existing free space. Raw inputs, genome, index and alignments stay in ignored storage. Record source URL, byte count, ETag, retrieval time and SHA-256; require complete decompression and valid FASTQ records. Source checksums are distinguished from newly calculated hashes. The primary assembly includes unplaced/unlocalized sequence but excludes alternative haplotypes and patches; this limitation remains explicit.

## Alignment and counting

Use approved Minimap2 2.31-r1302 with the long-read `splice` preset, both candidate transcript orientations, eight mapping workers, 50-million-base query batches, secondary alignments retained, and no transcript annotations or junction bonuses in the primary alignment. Index the complete primary assembly in one part. A gene-only mapping cannot establish uniqueness. Do not deduplicate PCR reads as though molecular barcodes existed. Retain target-read sequences and all of their output alignments, including secondary/supplementary alignments to other loci, in ignored storage; retain aggregate and read-level audit tables in the repository.

The primary local-junction category requires an exact skipped-reference (`N`) CIGAR boundary, primary non-QC-failed alignment, MAPQ ≥20, no supplementary alignment/SA tag, positive inferred transcript strand where supplied, and at least 20 contiguous aligned query/reference bases on each junction side. Adjacent mismatches count toward anchor length; insertions/deletions interrupt the strict anchor. Record all excluded candidates and reasons. Report the fixed ≥10-base anchor sensitivity separately. A deletion is not a splice, and low coverage is not an absence result.

For linkage, require the exact target and downstream chain on the same alignment, with qualifying anchors at each required junction. A partial read can establish local linkage but not a complete RNA ending. Distinguish a read containing the first discriminating downstream junction, one covering the entire downstream chain, and one covering the full annotated splice chain plus the shared coding-start codon. For full-chain matches, require at least 20 aligned bases in the model's terminal exon; the catalog's exact transcription start/end are not demanded or inferred from fragment ends. Exact sequence identity, aligned start/end and read names are technical diagnostics, not independent-molecule counts.

## Artifact and sensitivity checks

Before real target reads, verify synthetic normal/+29/T29 sequences, reverse complements, indel-versus-intron coordinate handling, short anchors, and chimeric/secondary exclusion. Compare recovered raw target junctions under a second Minimap2 run with `--splice-flank=no`; retain disagreements rather than using a favorable setting. For each event-bearing candidate, compare the observed local sequence against independently constructed normal and +29 cDNA junction templates with a separate dynamic-programming alignment, reporting score margins without treating an uncalibrated score as a probability. This is algorithmic sensitivity, not an independent sequencing experiment.

Record per-read alignment identity, clipping, supplementary/secondary hits, exact duplicate sequence hashes, poly-A/T end runs, and internal A-rich genomic windows near the terminal alignment. These are warning features, not proof of a biological terminus or automatic grounds for excluding all support. Manual inspection of exceptional reads can explain categories but cannot change the frozen primary count without a versioned amendment.

## Reporting boundary

Report both libraries separately, including zero and unqualified categories, normal evidence, other donor boundaries and incomplete endings. Counts represent amplified cDNA sequencing reads from one donor. Do not compute donor-replication P values, infer genotype from noisy retained bases, infer NMD response, or claim that a full-length translated protein exists. Any additional subset, changed threshold or new method is a separately declared extension.

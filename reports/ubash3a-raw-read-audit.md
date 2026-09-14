# Initial TRAILS RNA-read audit: +29 not confirmed

The complete naïve and memory CD4 libraries contain **8,401,078 reads**. Under the fixed primary analysis, **80 reads support the normal junction and zero support the exact +29 junction**. No exact +29 alignment appears even before the quality filters, including secondary/supplementary records. The second alignment setting agrees for the nominated target reads.

This completes the bounded two-library audit in [computational direction 1](../docs/computational-followup-roadmap.md). It does not confirm the TRAILS alternative-ending catalog candidate in these libraries. It does not establish biological absence, RNA degradation, or a lack of protein in other samples.

## What the reads show

The [DRA016285 archive](https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA016/DRA016285/) explicitly links DRR465268 to NCD4 and DRR465269 to MemCD4. Both are amplified cDNA libraries from one donor. Read counts are neither original-molecule counts nor independent donor replication.

| Measurement | Naïve CD4 | Memory CD4 |
|---|---:|---:|
| Complete input reads | 5,611,025 | 2,790,053 |
| Complete input bases | 5,469,014,753 | 2,581,343,875 |
| Read groups with any locus alignment | 397 | 248 |
| Exact normal-junction alignments, before quality filters | 86 | 49 |
| Normal-junction reads passing the primary criteria | 47 | 33 |
| Normal-junction reads at the fixed 10-base-anchor sensitivity | 57 | 36 |
| Exact +29 alignments, before quality filters | 0 | 0 |
| Exact T29-specific downstream-junction alignments | 0 | 0 |
| Full normal U0 splice chain plus aligned shared start, primary criteria | 0 | 1 |

Coordinates are GRCh38, chromosome 21, positive-strand, one-based inclusive exon endpoints. Normal is **42434954→42437488**; +29 is **42434983→42437488**. T29's distinguishing downstream junction is **42437580→42437855**. The naïve normal total includes one secondary/supplementary alignment; the other 85 and all 49 memory alignments are primary.

All 80 qualifying normal reads remain qualifying under the second alignment setting. Three naïve and two memory reads cover the complete normal downstream chain shared by U0, C0 and C20. These are the same reads across those model labels and must not be added together. At the 10-base-anchor sensitivity, the counts are 13 and nine; three memory reads cover the full U0 chain plus shared start. No full or downstream +29 reference chain, including T29, qualifies at either anchor threshold. The boundary model's downstream chain is also not supported with normal splicing.

See the [complete summaries](../data/derived/ubash3a-raw-read-audit/raw-read-summary.csv), [architecture counts](../data/derived/ubash3a-raw-read-audit/raw-architecture-counts.csv), and [all-setting event counts](../data/derived/ubash3a-raw-read-audit/raw-all-setting-event-counts.csv). The [junction-neighborhood table](../data/derived/ubash3a-raw-read-audit/raw-junction-neighborhood.csv) preserves other exact boundaries without rounding them to the target.

## Method and checks

The [plan](../config/ubash3a-raw-read-plan.json) was frozen before target-read inspection. Minimap2 2.31-r1302 aligned every input read against all 194 sequences of Ensembl release 112's GRCh38 primary assembly, without supplied transcript annotations. All alignments for nominated target queries were retained, including competing loci. The reference excludes alternative haplotypes and patches.

Primary eligibility requires MAPQ ≥20, a primary alignment without QC/duplicate flags or supplementary/chimeric evidence, no inferred negative transcript strand, and ≥20 contiguous aligned bases on each side of the exact junction. The fixed ≥10-base-anchor analysis and `--splice-flank=no` run are sensitivities. The second run examines nominated reads only; it cannot recover queries missed by primary-setting nomination. Full-chain counts require aligned shared-start coordinates and terminal-exon coverage, but do not certify physical RNA ends or an error-free coding sequence.

All 134 evaluable primary normal-junction reads favor the normal sequence in the separate local dynamic-programming comparison. The +29-minus-normal score ranges are −33 to −25 in naïve and −33 to −22 in memory. There are no exact +29 candidates to score. These uncalibrated scores are not probabilities, and this is not an alignment-free search of every unmapped read.

Complete FASTQ validation, bzip2 integrity and input/output query totals pass. Read/base totals agree exactly with independently retrieved ENA records. The [independent verifier](ubash3a-raw-read-verification.json) checks 29 sample/run joins, all 661 locus alignment rows using a separate textual-CIGAR implementation, and 20 deterministic local-score comparisons. Eight [synthetic forward/reverse examples](ubash3a-raw-read-synthetic-verification.json) recover their expected junction chains and genomic transcript strand. All 212 repository tests pass; the four primary tables and audit replay byte for byte. [Reproduction instructions](../docs/ubash3a-raw-read-reproduction.md) identify required source caches and the separately declared coding-coordinate descriptors.

## Consequence for the working hypothesis

The preceding [catalog audit](ubash3a-hypothesis-test.md) preserves T29-assigned matrix values of 14 for naïve and 16 for memory CD4. Those are source-assigned count units with unresolved export/normalization details, not counts of reads traversing T29's unique junctions. This raw audit does not reproduce direct junction support in either selected library. Assignment ambiguity, assembly behavior and coverage limitations are possible explanations; none has been established here.

The [conditional 482-aa T29 reconstruction](ubash3a-alternative-ending.md) remains a sequence hypothesis. Its final-exon stop is not evidence of RNA survival or protein production. Other TRAILS subsets and donors were outside this fixed initial scope. Protein-evidence work can proceed as a separately qualified test with its own specificity and error controls; no mass spectra have yet been searched.

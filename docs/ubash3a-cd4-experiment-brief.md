# UBASH3A RNA survival: laboratory feasibility brief

Prepared September 14, 2026. **Status: experimental design prepared; laboratory, samples and bench conditions not yet arranged. No new biological measurements have been collected.** This brief requests a staged feasibility assessment before a confirmatory study. It is not a validated laboratory SOP or an order for materials.

## The question

Does the exact +29-nucleotide splice form of UBASH3A occur in a complete RNA with the ENST00000398367.2 downstream structure, and does that RNA undergo less nonsense-mediated decay (NMD) than the corresponding canonical-derived +29 RNA in primary human CD4 cells?

The practical outcome is an identified RNA and an interpretable decay comparison. An increase in total UBASH3A, a short +29 amplicon, or a predicted protein structure cannot answer this question. A result in healthy-donor CD4 cells would establish a mechanism in that experimental context; it would not establish an effect on PSC or a treatment benefit.

The prior public-data audit found a +29 catalog candidate with a different terminal exon. Its assigned CD4 support comes from one donor. The available UPF1-depletion libraries lacked diagnostic reads, so the survival contrast was not measurable. Those results motivate the targets below; they are not a positive NMD result.

The sequence rationale is conditional: the reconstructed B29 stop occupies 42442545–42442547, placing its end 49 nucleotides before the final exon junction at 42442596 (51 nucleotides from the stop's first base). Other reference structures retain more downstream junctions. This is a positional boundary hypothesis, not demonstrated NMD escape. T29 ends before that stop region and needs separate interpretation.

## RNA identities that the assay must distinguish

All coordinates are GRCh38 chromosome 21, positive strand. Junction notation is **last exonic base → first base of the next exon**, both one-based. It is not an intron BED interval.

| Label | Reference or source | Identifying structure | Role |
|---|---|---|---|
| B29 / B0 | ENST00000398367.2, with +29 / normal target splicing | Common downstream junction 42437580→42442452, followed by 42442596→42447057 | Original boundary hypothesis and its normal-splicing counterpart |
| C29 / C0 | Canonical ENST00000319294.11, with +29 / normal splicing | Additional upstream exon 42413410–42413523; downstream junctions through 42443312 and 42444534 | Original canonical comparison; upstream coding structure must also be identified |
| U29 / U0 | ENST00000291535.12, with +29 / normal splicing | Same reference upstream exon structure as B; downstream junctions through 42443312 and 42444534 | Secondary comparison that reduces the upstream-structure difference in B versus C; this model is **not** the canonical transcript |
| T29 | TRAILS e09f04f5-c617-4bb3-b69c-519ccb16d0f1 | 42437580→42437855, then terminal exon to 42440816 | Separate exploratory candidate; its deposited NMD label is a prediction |
| Other / unresolved | Other published reference models, novel ends, partial molecules, ambiguous assignments | Preserve their actual exon chains | Never assign these to B, C or T merely because they contain +29 |

The normal target junction is **42434954→42437488**; +29 is **42434983→42437488**. The retained bases are **42434955–42434983**. rs1893592 at **42434957**, analyzed previously as A/C, lies inside the retained segment and is removed by normal splicing. Normal RNA cannot report that allele at this site. A/C reference sequences are hypothetical templates, not donor genotypes; other alleles remain separately labeled.

The accompanying target table contains every exon chain for ten reference structures and the one TRAILS candidate. The supplied 12 short planning windows establish local junction linkage only. The 278/133-nt windows omit rs1893592; the 297/152-nt versions include it. **None is a validated primer or proof of a full RNA.**

## Stage 1: establish feasibility and identify the RNA

The proposed starting pilot uses **four independent, coded healthy donors**, with purified CD4 cells split into a resting condition and one standardized activated condition. This is a feasibility allocation, not a powered sample size or a genotype-association study. Record cell purity, activation state and viable-cell recovery; paired aliquots and repeated sequencing runs remain one donor. The laboratory should state whether existing appropriately sourced donor material can support this allocation.

Qualify one acute NMD suppression method using pathway controls before interpreting UBASH3A. A selective SMG1 inhibitor is a candidate method, subject to validation in these cells. Match its vehicle control and establish a usable exposure window from target engagement, at least two independently validated NMD-substrate assays, non-NMD comparison transcripts and cell-state measurements. The responsible laboratory must specify the actual compound identity, dose, timing and acceptance thresholds. Cancer-cell experiments support the general perturbation approach, but do not supply validated CD4 conditions. See [the SMG1-inhibitor study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11832170/).

Each donor/state is split between control and the qualified perturbation: **4 donors × 2 states × 2 conditions = 16 biological aliquots**, before additional qualification controls. Separate reverse-transcription reactions from each RNA preparation provide a technical reproducibility check. This pilot can reveal a transcript that becomes detectable only when NMD is suppressed; that observation alone is not its measured decay rate. Do not transfer the earlier auxin treatment to ordinary CD4 cells: that study used engineered degron-tagged UPF1 cell lines. [Boehm et al.](https://pubmed.ncbi.nlm.nih.gov/40934927/)

Use targeted long-read RNA/cDNA sequencing with enrichment across shared UBASH3A regions and all three candidate endings. The library must retain the target junction, all discriminating upstream and downstream junctions, and evidence of the actual 5′ and 3′ ends on the same original molecule. Use molecular identifiers introduced before amplification if the selected workflow supports validated molecule counting. Preserve raw reads, consensus families and partial-read exclusions. A consensus or assembly made from overlapping molecules is not evidence that the complete chain occurred on a single molecule.

The reference coding start is 42403946. Reads that miss the upstream coding region cannot distinguish all of the reference RNAs. Independently confirm unexpected ends and junctions using a second RNA preparation/RT and an appropriate end-mapping or amplicon-sequencing assay. Separate short 5′ and 3′ assays cannot by themselves connect those ends to +29. Targeted long-read enrichment has been demonstrated for immune transcripts, but the UBASH3A assay still needs its own validation. [Penter et al.](https://www.nature.com/articles/s41467-023-44137-7)

Include no-template and no-RT controls, normal-splicing and +29 template controls, both A/C template versions where relevant, and mixtures of competing architectures to measure cross-assignment. Validate extraction/enrichment recovery, amplification bias, genomic-DNA contamination, internal priming and chimera formation. Controls should be distinguishable from biological molecules. Do not design an enrichment primer that requires the canonical terminal exon and thereby excludes T29. Whole-transcriptome read counts alone are unlikely to solve the previously observed coverage problem.

**Proposed structural progression rule:** both primary +29 RNAs must be independently resolved in at least two of the four donors in the same cell state, with each positive supported in two independent RT preparations by at least three distinct original-molecule families per preparation. These are proposed conservative pilot criteria, not a validated detection limit. Final molecule/error thresholds must first be justified by dilution and negative-control experiments; UMI families qualify only when they demonstrably represent distinct input RNA molecules. No threshold may be relaxed after seeing the desired architecture. If only local linkage is established, report that narrower result.

If B29 is not detected, record the number of eligible UBASH3A molecules, measured detection limits and state/condition. Do not infer its biological absence. If only T29 is verified, continue it as a new exploratory question; it cannot replace B29 in the original test. A non-detection interval applies only to the demonstrably sampled molecule population and its recovery assumptions.

## Stage 2: measure decay, with donors as the comparison unit

Proceed after structural and assay qualification. The primary comparison remains **B29 versus C29**. U29 is a predeclared secondary comparison because it shares B's reference upstream structure; this improves interpretation of the downstream difference without silently changing the original question. Measure the corresponding normal-splicing RNAs independently rather than assuming they are stable normalization controls. ENST00000635325.1 already carries an NMD annotation and is not a presumed negative control.

For each donor, split the same cell state between control and the qualified NMD perturbation. Use a metabolic RNA pulse/chase or another laboratory-validated direct turnover assay, with a baseline and at least four informative post-chase times. Determine the actual schedule in feasibility work so it covers the observable decline. Track the pre-existing labeled **mature RNA**, including recovery controls, residual label incorporation, maturation of precursors, changing cell numbers and cell division. Report an apparent turnover rate if degradation cannot be separated from these contributions.

The turnover readout must retain **both RNA age and full-structure identity**. Default gene-level or 3′-tag SLAM-seq is insufficient for these shared-end isoforms. A candidate implementation is labeled-RNA enrichment followed by targeted, molecule-resolved sequencing; another is an appropriately validated isoform-resolved metabolic sequencing workflow. A short calibrated assay may quantify a structure only after its specificity is established against every observed competing RNA. The laboratory must demonstrate recovery and specificity at the target's abundance before selecting either route. Metabolic labeling has been used to estimate RNA stability, including in human CD8 cells, but adapting it to this primary-CD4 target remains feasibility work. [Herzog et al.](https://www.nature.com/articles/nmeth.4435), [human CD8 study](https://www.nature.com/articles/s41467-025-67762-w), [nano-ID study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7545145/).

Freeze the final assay, cell state, times, donor number, exclusion rules and analysis before confirmatory donor measurements. Use pilot recovery and variance to calculate the required donor count for a specified meaningful effect; four donors cannot be presented as sufficient power. If there are N donors and T time points in one state, the two-arm decay study requires **N × 2 × T biological harvests**, plus controls. Include these harvests and cell input in the laboratory's feasibility response.

The [analysis plan](ubash3a-cd4-analysis-plan.md) specifies the contrasts and failure states. A second, independently qualified way of suppressing the NMD pathway is required to strengthen a single-perturbation result into an NMD mechanism claim. A lone UPF1 or SMG1 response remains dependence on that intervention, with possible effects outside NMD. A claim that the downstream junction itself causes the difference would additionally require a matched, splice-competent sequence intervention/reporter that preserves the relevant processing context; an intronless cDNA construct does not reproduce the proposed exon-junction-complex mechanism.

## Stage 3: protein, after an RNA is verified

Use an orthogonal protein assay with sequence-specific detection and appropriate normal-product controls, such as targeted mass spectrometry with reference peptides plus an antibody-based size assay where validated. Validate peptide uniqueness against the current human proteome and all candidate products before interpretation. Several of our conditional RNAs encode the same altered protein: a shared peptide cannot identify its source RNA. An altered peptide establishes translation evidence; persistence or a protein half-life requires a protein-turnover measurement. Protein detection is not proof of preserved UBASH3A function.

## What the laboratory needs to return before execution

1. Named scientific lead, specimen route, available donor/cell input, and the applicable institutional requirements for that material.
2. A proposed platform/enrichment approach and a demonstration that molecule identity and RNA age can be joined at the required sensitivity.
3. The qualification plan for NMD suppression, cell state, molecular counting, controls and detection limits; finalized primers/probes follow that validation.
4. Separate feasibility and confirmatory-study estimates, including donor material, culture/perturbation, sequencing, assay development, orthogonal confirmation and raw-data delivery. No cost or completion date is assumed here.
5. A completed method-freeze record and a secure route for coded raw measurements. Populate the supplied templates in the laboratory's study storage; they contain no patient or donor data now.

This packet is prepared for that feasibility discussion. No laboratory has been contacted, no specimen access has been claimed, and no material has been ordered.

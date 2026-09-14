# UBASH3A: public-data test of transcript-dependent RNA survival

## Question and current knowledge

Does the exact +29 splice form occur in complete human RNA isoforms, and does its downstream structure alter sensitivity to nonsense-mediated decay (NMD)? The earlier sequence reconstruction predicts a boundary case for ENST00000398367.2 and more clearly positive positional flags for four other models. Those calculations are prior information, not new experimental evidence.

This follow-up uses existing public observations. It may test only a necessary premise if a suitable decay experiment is unavailable. RNA annotations, read counts, an NMD perturbation response, ribosome occupancy and stable protein detection are separate evidence levels. No researcher contact, new enrollment, laboratory order or model prediction is included.

## Source discovery and selection, before target values

Two initial web searches found long-read CD4 references and long-read NMD perturbation resources. Only study descriptions, file inventories and general methods were inspected before this freeze; no new UBASH3A rows, transcript structures, abundances or responses were inspected. Search history and this boundary are recorded explicitly; the entire search is not described as preregistered.

The initial source families are:

1. Activated primary human CD4 T cells, GSE229971: four time points from one donor.
2. Resting primary human CD4 T cells within GSE202327/GSE202329.
3. TRAILS human immune-cell long-read reference, including its independent PBMC annotation and public isoform count table.
4. Rapid UPF1-depletion / consolidated human NMD transcriptome resource, including long-read annotation and accessible experimental data.
5. Leshem et al. primary-human-lung long-read SMG1-inhibition dataset, GSE329233 and its versioned public data archive.

Resolve current source versions, assembly, sample identities, filtering and count definitions before target interpretation. A bounded targeted PubMed/GEO check may add up to two directly relevant human CD4 plus NMD-perturbation datasets. Include every selected source, including unavailable, low-coverage and negative annotation results; do not select datasets according to UBASH3A results. Do not treat multiple algorithms, cell subsets or time points from one donor as independent donors.

Prefer deposited transcript coordinates and quantification tables. Download the inputs needed for this locus, not whole software environments or unrelated model outputs. If a raw-read audit is necessary, first establish public availability and a workable indexed or otherwise bounded extraction. A reference-guided or filtered transcript catalog cannot by itself establish that an unreported splice form is biologically absent.

## Fixed RNA endpoint and classification

Keep the previous GRCh38 positive-strand endpoints: normal upstream exon end 42434954, retained upstream exon end 42434983, next exon start 42437488. Compare exact consecutive exon pairs. Search the whole source locus by coordinates as well as gene IDs, retaining discordant labels separately. Nearby splice junctions, complete intron retention and exon skipping do not count as the exact +29 event.

Construct a reference junction-chain specification from all five prior transcript models. Keep the normal/+29 version of each and retain full exon coordinates. For an observed model report: target junction status; exact full internal-junction-chain agreement; agreement only downstream of the target; coverage of the reference coding start and new stop; terminal-coordinate differences; strand/assembly; and the original source ID and support fields. Matching downstream chains alone does not establish an entire named transcript. A single +29 junction does not link it to a particular downstream chain.

Count observed source transcript models separately from supporting reads, molecules, samples and donors. Preserve original support units and fractional estimated counts; never relabel them as direct molecule counts. For catalogs using junction correction or reference-based inference, describe that dependency and inspect the corresponding uncorrected or independent annotation if available. Absence in a filtered catalog is an annotation-level non-detection, not a zero molecular abundance or a disproved hypothesis.

## Conditions for a decay-response test

Before examining target response values, require explicit treated/control identities and independent replicate or donor correspondence, a documented NMD perturbation, linked +29/downstream transcript structure, and a usable count/abundance table with its normalization and filtering defined. CD4 observations are the primary biological context; other human cells provide a separate transfer test and cannot substitute for CD4 validation.

If a source passes these metadata gates, freeze its sample/contrast specification and statistical procedure in a separate dated amendment before looking at UBASH3A response values. Retain all planned isoforms and zero/low-count outcomes. If target isoforms are not identifiable or covered, report the intended contrast as not estimable. Do not use gene-level UPF1 responsiveness or an annotation-based NMD label as the requested isoform-specific test. UPF1 perturbation alone is not automatically a specific measurement of decay rate.

Protein-sequence distinguishability must also be checked. If several altered RNAs encode the same peptide, detecting that peptide cannot identify the originating RNA. Translation evidence and stable protein evidence remain distinct. New structure modeling does not fill these measurement gaps.

## Acceptance and verification

Deliver the pinned source inventory, reference discriminators, every eligible observed UBASH3A structure and its original support fields, an explicit dataset-by-question assessment, and any estimable frozen decay contrast. Otherwise state exactly which premise was tested and why the full hypothesis remains untested. Preserve all previous analyses. Validate new coordinate/parsing logic against synthetic truth and independently check real-data locus extraction and chain matches. Run the repository suite after parsing changes, record hashes and update the report/evidence navigation. No biological claim should exceed the evidence level reached.

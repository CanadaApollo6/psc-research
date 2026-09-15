# T29 peptide specificity and public-spectrum qualification

This is computational direction 3. The candidate is the previously fixed TRAILS alternative-ending T29 sequence, conditional on the shared translation start. Its exact junctions were not confirmed in the two preceding TRAILS read libraries. The present analysis asks a separate protein-evidence question and does not infer RNA stability.

## Fixed sequence question

Before inspecting new peptide matches, retain all four previously generated targets: RPNSSWKLGIFK, LGIFKIHLESGIR, IHLESGIR and IHLESGIRKP. They arise under the earlier tryptic rule (zero or one missed cleavage, 7–35 residues). Do not change the candidate list or length/cleavage criteria after viewing results.

Search the complete retrieved human UniProtKB sequence set, including reviewed/unreviewed entries and available curated isoforms, plus GENCODE release 50's comprehensive protein translations and a separately pinned common-contaminant collection. Retain the twenty earlier UBASH3A modeled products and both T29 allele products as explicit reference groups. Record source release, retrieval time, sequence counts, completeness evidence and hashes before searching targets. If a planned source cannot be obtained or validated, report the missing scope rather than treating it as an empty set.

Perform exact substring matching and a separate comparison with I and L treated as equivalent. Retain every matching record and position, sequence identity, source header, flanking residues and tryptic compatibility. Substring matches outside expected tryptic boundaries remain possible sequence competitors and are not discarded. Scan complete sequences before digest filtering. A nonmatch means absence from these retrieved sequences, not uniqueness across all possible human variants, unannotated proteins or organisms. Ambiguous source residues and uncertain translation models remain explicit limitations.

Independently verify target matches with a different search implementation and recheck all candidate sequences against their parent proteins. Preserve overlapping/nested peptides as such; four peptide strings are not four independent pieces of protein evidence. Shared peptides cannot assign a source RNA or distinguish the hypothetical A/C alleles.

## Public-file qualification

Refresh PXD000376 project metadata and use the already complete file inventory. Read the primary paper and supplemental methods, and inspect deposited pepXML search metadata to reconcile sample/fraction labels, instrument settings, enzyme/modification rules, database identities and file references. Metadata and prior database-search identifications are distinct from a new search of spectra. Do not infer donor count from file names or use old-database nonmatches as protein absence.

The initial qualification transfer is bounded to protein references, supporting methods, and deposited search metadata (at most 2 GB). Do not transfer all 92 GB of raw spectra as a default. After metadata and reference checks, fix a coherent public experiment/fraction set and its resource budget before any new target-spectrum scoring. State which sample identities are observed, inferred or unresolved. A database or sample-label uncertainty can limit attribution without preventing all useful technical checks.

## Spectrum-search gate

A new spectrum search requires readable spectra, a complete frozen target/background/contaminant database with decoys, source-supported digestion/modification and mass-tolerance settings, and declared error controls before inspecting candidate scores. Background and ordinary UBASH3A peptide evidence must be retained. A global 1% error estimate alone does not establish credibility for a small novel-peptide class. Candidate spectra need explicit competing-assignment and fragment evidence; one shared peptide is not proof of the proposed protein. Document any inability to measure this endpoint.

The outcome may be qualified sequence targets, evidence of non-specificity, a bounded new spectrum analysis, or a concrete unresolved measurement gate. A missing match does not prove that T29 protein is absent. No personal genotype, clinical effect, new laboratory experiment or researcher outreach is part of this work.

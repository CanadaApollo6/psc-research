# Working in this repository

## Scope

Help carry out the public-data PSC research protocol in `docs/research-plan.md`. Independent analysis does not require a researcher reply. Keep the completed work, proposed analyses, and biological conclusions clearly distinguished.

## Research practices

- Read `docs/evidence-log.md` before interpreting the imported tables.
- Preserve original source values and record corrections separately with supporting evidence.
- Record source URLs, retrieval dates, genome builds, allele conventions, model versions, and transformations.
- Do not equate lead SNPs, fine-mapped SNPs, credible-set members, or functional-QTL variants.
- Keep gene labels as published until identifiers are independently reconciled.
- Distinguish model predictions, observational associations, experimental results, and clinical outcomes.
- Evaluate at the locus or donor level when appropriate; correlated variants and cells from one donor are not independent replicates.
- Never invent missing coordinates, alleles, predictions, collaborators, or study-access claims.

## Workflow

- Use `data/raw/` for downloaded inputs, `data/derived/` for small reproducible tables, and `reports/` for findings.
- Keep personal health information, credentials, private correspondence, and large model outputs out of Git.
- Do not contact researchers or enroll in services without the user's authorization.
- Run `python -m unittest discover -s tests -v` after changing the import logic. Rebuild the catalog when source parsing changes and inspect the resulting diff.
- Do not use a model-generated score as a clinical recommendation.

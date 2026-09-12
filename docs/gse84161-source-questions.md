# GSE84161 source questions — prepared, unsent

Study: Senger et al., *Science Signaling* (2017), [DOI 10.1126/scisignal.aah4273](https://doi.org/10.1126/scisignal.aah4273); monocyte series [GSE84161](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84161), superseries GSE84162.

These are research-provenance questions for the corresponding study team or GEO submitter. No contact has been made and no recipient address has been guessed. The [study report](../reports/gse84161-study-investigation.md) and [sample map](../data/derived/gse84161-samples.csv) make the request reviewable.

1. **Anonymous donor grouping.** Do replicate numbers 1–5 in the monocyte records each correspond to one distinct donor across all six conditions? For example, are GSM2228386–GSM2228391 all donor group 1, with successive groups of six through GSM2228415? An anonymous sample-to-donor key or confirmation is sufficient; no personal identifiers are requested.
2. **Chemicals and array-specific doses.** What exact compounds, structures or catalog identifiers were used for the TPL2i and MEKi arms? GEO gives 2.8 µM and 0.31 µM, while Figure 2D gives 926 nM and 103 nM. Which concentrations generated the deposited monocyte CEL files, and is the selectivity supplement available?
3. **Historical probe mapping.** Which annotation release and Entrez mapping produced the 19,944 retained probe sets? What sample universe was used for highest-variance selection, and how were ties and multiple Entrez mappings handled? The current GPL570 snapshot maps 267 retained probes to no Entrez ID and 414 to multiple IDs; these differences may reflect annotation versions.
4. **Deposited quantity and processing code.** Are the per-sample VALUE entries the RMA summaries before donor adjustment, donor-model residuals, or another transformed quantity? Could the original R scripts, software versions and any array QC/exclusion records be shared? We understand that the paper separately describes donor adjustment and per-gene scaling for heatmaps.

Independent work can continue on original-array normalization, fixed annotation and program coverage while these questions remain open. A reply would support donor-paired interpretation and exact source reproduction; it is not a request for access to a clinical treatment.

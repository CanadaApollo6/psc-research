# Data provenance

Snapshot: 2026-09-11. The exact URLs, output filenames, byte limits and SHA-256 hashes are recorded in [config/sources.json](../config/sources.json).

The later [QTL follow-up manifest](../config/qtl-followup-sources.json) adds 39 pinned public sources; [its query manifest](../config/qtl-followup-queries.json) records eight indexed extractions and all 121 retained rows. Source metadata uses eQTL Catalogue release 7 from resources commit `65995fe96c8c6da1e750b51bac3b4c2f5ba9260f`; the queried archive files are pinned through their headers, indexes and exact subsets, without claiming a hash of the full remote files. The [source report](../reports/qtl-followup.md) distinguishes reprocessed BLUEPRINT, DICE, Lepik and TwinsUK cohorts, their RNA endpoints and missing measurements. The Catalogue API is retired; its official [archive access page](https://www.ebi.ac.uk/eqtl/Data_access/) documents indexed queries and ALT effect-allele semantics.

LeafCutter phenotype metadata comes from [Zenodo record 7850746](https://zenodo.org/records/7850746), Ensembl 105. The conversion to model junction intervals is explicit, and conflicting cluster/gene strands are preserved separately. The [AlphaGenome training-track extract](../data/derived/qtl-comparison-training-tracks.csv) comes from published Supplementary Table 2, not an inferred inventory of the live service. Raw papers, metadata workbooks and compressed indexes stay ignored; small derived tables retain source attribution.

The [direct-junction manifest](../config/ubash3a-junction-sources.json) pins 35 sources for the next analysis: two exact Snaptron junction responses, recount3/ENA study metadata, direct sample-ID responses, public 1000 Genomes metadata/index, reference DNA, workflow code and data-access evidence. The [run manifest](../config/ubash3a-junction-run.json) separately pins the one-variant VCF subset and the plan saved before target counts were inspected. The VCF's missing rsID is preserved; exact GRCh37 coordinates and A/C alleles match the previously verified mapping. All 667 RNA runs are matched explicitly and combined into 464 donor records, of which 447 have public genotypes. Public sample-level joins remain ignored; only aggregate results are versioned. [Full methods and limitations](../reports/ubash3a-junction-followup.md).

## Imported inputs

| Input | Source | Transformation |
|---|---|---|
| Article XML, PMC11541731 | [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11541731/fullTextXML) | Extract Table 1 and expand merged cells so secondary signals keep the correct locus. |
| Supplementary Data 1–3 | [Goode et al. article attachments](https://www.nature.com/articles/s41467-024-53602-w) | Cache all three small files; import Data 2's colocalisation evidence. Data 1 and 3 are supporting source material. |
| GSE243977 SOFT metadata | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243977) | Decompress; extract sample fields; derive explicitly provisional title labels and flag selected cross-field differences. |
| GSE247128 SOFT metadata | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247128) | Same metadata-only transformation. |

Hashes apply to the exact files stored in `data/raw/`, **after decompression** for SOFT files. Downloads are bounded in size and fail on a changed hash. Raw files and the run-specific retrieval log are ignored by Git. Derived outputs can be regenerated with the documented commands.

## Attribution and reuse

Goode EC, Fachal L, Panousis N, et al. *Fine-mapping and molecular characterisation of primary sclerosing cholangitis genetic risk loci*. Nature Communications 15, 9594 (2024). DOI: [10.1038/s41467-024-53602-w](https://doi.org/10.1038/s41467-024-53602-w). The article is licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Derived CSVs change the format and expand merged cells; source content and attribution are retained. Third-party material can have separate terms.

The GEO records originate from the [PSC liver-atlas study](https://pubmed.ncbi.nlm.nih.gov/38199298/) and its submitting investigators. Accession and sample links remain attached to the derived rows. Public access does not grant unrestricted rights to every linked dataset; check each source's terms before reusing additional material.

Individual-level genotype and expression files at EGA can require controlled access. They are not downloaded or required for this initial import. Large expression matrices, personal genomes and private medical records are not included.

## Variant-reference and prediction provenance

[config/reference-sources.json](../config/reference-sources.json) pins 27 additional responses. Ensembl release 116 provides variant mappings on both builds, individual reference bases, and genes overlapping fixed 1-Mb neighborhoods. Build-38 reference bases were also retrieved directly from AlphaGenome's GRCh38.p13 FASTA using its `.fai` index and exact HTTP byte ranges. Both the index and the one-base responses have recorded hashes. All alternative mappings remain in raw inputs; ambiguous primary mappings are rejected by the normalizer.

The variant table distinguishes published risk labels from reference/alternate alleles. Candidate-gene tables are derived from the pinned gene annotations. The [run provenance](../reports/pilot-run-provenance.json) records client commit, selected model version, protocol/input hashes, metadata hash, request dates and prediction-file hashes. The [pilot report](../reports/first-alphagenome-pilot.md) includes compact model-derived results; full output tables remain in ignored local storage.

AlphaGenome access and outputs have [separate terms documented by the official project](https://github.com/google-deepmind/alphagenome#terms). The current work is noncommercial research. Do not assume a future code license would override source-data or model-output conditions.

## Mechanism follow-up sources

[mechanism-sources.json](../config/mechanism-sources.json) pins article XML/HTML, the 2020 Goode thesis, GWAS Catalog study/association metadata, and Ensembl transcript/variant context. Where Europe PMC did not expose XML, the public PMC HTML was used. No browser challenge was treated as article content. The thesis's relevant allele passages on pages 60 and 80 were rendered and visually checked.

[association-sources.json](../config/association-sources.json) pins the original GWAS file header, directory metadata, 23 small byte ranges and rs313839 population frequencies. `fetch_gwas_audit_rows.py` verifies sorted complete interior lines, skips partial boundary fragments and requires an exact rsID/coordinate match. The two selected lines are copied without changing their original fields. `build_allele_audit.py` verifies those copies against pinned ranges and derives separate reviewed fields. The 662,022,204-byte file is not downloaded in full. Its header defines `allele_1` as risk; sample scope and conflicting source labels are described in the audit.

[splice-audit-sources.json](../config/splice-audit-sources.json) pins the official model's coordinate implementation to research-code commit `0db53bd4352c66d1e00a049a81da373a066e6670` and the 2020 Mucaki splicing paper. This code was read to confirm coordinates; it was not used to run a different model. The API client remains pinned separately.

The one additional model request has [separate provenance](../reports/splice-followup-provenance.json). Its original v0.1 protocol is archived unchanged; v0.2 records the post-request donor-index correction. Complete output arrays, native junction coordinates and metadata are saved locally with hashes. Summary joins use feature coordinates and assay identities rather than row order; unreturned features remain missing. The downstream donor is the largest positive change within the originally fixed 101-base window, selected for exploratory interpretation. It is not an independently chosen validation endpoint.

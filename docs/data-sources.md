# Data provenance

Snapshot: 2026-09-11. The exact URLs, output filenames, byte limits and SHA-256 hashes are recorded in [config/sources.json](../config/sources.json).

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

AlphaGenome access and outputs have [separate terms documented by the official project](https://github.com/google-deepmind/alphagenome#terms). No outputs are included in this snapshot. Do not assume a future code license would override source-data or model-output conditions.

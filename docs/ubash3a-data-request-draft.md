# UBASH3A data request — draft, not sent

Prepared 2026-09-11. Suggested public professional contact: eQTL Catalogue, eqtlcatalogue@ebi.ac.uk. Kaur Alasoo's public [2023 response](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/issues/35#issuecomment-1627249622) also offered targeted data transfers. These contacts have not agreed to collaborate.

**Subject: Small UBASH3A/rs1893592 splice-data request and QTD000377 strand question**

Hello Dr. Alasoo and eQTL Catalogue team,

I'm conducting an independent, reproducible public-data analysis of PSC-associated regulatory variants. Could you help with a small UBASH3A extraction from Lepik_2017 blood (QTD000377)?

For rs1893592 A>C (GRCh38 chr21:42434957), we need the unfiltered coefficients, standard errors and sample counts for the canonical `21:42434954:42437488:clu_31404_-` and +29-nt `21:42434983:42437488:clu_31404_-` junctions, including their cluster-denominator and normalization definitions. The published +29 row has beta 0.523483 and SE 0.063608. Per-genotype junction counts and donor counts would also help, if shareable; a small aggregate extraction would be sufficient initially.

The cluster suffix is minus, but the assigned UBASH3A strand and both GT–AG reference motifs are positive. In workflow v22.04.1, the reverse-stranded branch passes regtools `-s 2`, while the declared regtools 0.6.0 recipe defines RF as `-s 1`. Could you confirm the historical run settings and whether a corrected dataset exists? We have preserved the published labels and treat this explanation as provisional.

We also found no UBASH3A-assigned Leafcutter phenotype in the public DICE CD4 metadata (QTD000483/QTD000488). If these junctions were quantified before filtering, could the same small extraction be supplied? Our public GEUVADIS check was too sparsely covered to estimate the effect reliably.

I can provide the short audit, exact source rows and endpoint CSV. Thank you for considering this request.

[Your name]

---

Suggested attachments, reviewed before sending:

- [Endpoint specification](../data/derived/ubash3a-data-request-endpoints.csv): six rows, build and allele explicit; DICE cluster identifiers left unspecified.
- [Audit report](../reports/ubash3a-junction-followup.md): includes source links and limitations.
- [Motif check](../reports/ubash3a-junction-motifs.csv) and [conditional strand truth table](../reports/ubash3a-strand-flag-audit.csv).

The repository is private, so a repository link alone will not give a recipient access. Share only the selected report/tables after review. No API credentials, personal health information or sample-level genotype tables belong in this message.

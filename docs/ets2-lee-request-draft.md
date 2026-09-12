# Draft request to the ETS2 study team

Status: prepared for review; not sent. Correspondence address from the [published article](https://www.nature.com/articles/s41586-024-07501-1): James Lee, james.lee@crick.ac.uk. The sender's name is intentionally left for completion.

**Subject:** Reproducing the ETS2/MEK-inhibitor RNA-seq signatures: focused data/code question

Dear Dr. Lee and colleagues,

I am working on an independent public-data research project focused on primary sclerosing cholangitis. We have been checking the released gene signatures and pathway-analysis outputs from your 2024 Nature ETS2 paper before testing their transfer into other macrophage datasets.

Using your Zenodo v1.0 release, we recomputed all nine raw enrichment scores from the released 500-nM MEK-inhibitor rank vector within numerical precision. The archived NES values also match the published Figure 5 source data; we did not recreate the historical null distribution, NES normalization or P values. We also recovered the six released ETS2/enhancer gene lists from Supplementary Table S1 at adjusted P < 0.05.

We could not trace the exact export connecting the two MEK-inhibitor rank CSVs to the drug-minus-control contrasts in the released R script. Would you be able to share, or point us to:

1. The complete original 100-nM and 500-nM differential-expression tables before symbol collapse, retaining Ensembl IDs, symbols, logFC, t statistics, raw and adjusted P values, with confirmation of each contrast's orientation.
2. The code revision and annotation version used to filter/map genes and produce `PD_0.1_MEK.resSYMBOL.rank.csv` and `PD_0.5_MEK.resSYMBOL.rank.csv`, including missing/duplicate-symbol and tied-statistic handling, with the relevant `sessionInfo()` or package versions.
3. Whether the script's request for 12,801 topTable rows and the 12,095 symbols in each saved rank file reflect annotation filtering, duplicate handling or another step.

The requested complete per-gene tables and export details would be sufficient for this provenance check; we are not requesting participant-level data in this email. We have a short reproducibility brief and the numerical comparison table ready to share if useful.

Thank you for making the code and results available.

Best regards,
[Your name]

## Companion material

The [compact reproducibility brief](ets2-lee-reproducibility-brief.md) contains exact source IDs, the successful checks and the remaining questions. Its findings table is [ets2-inhibitor-reconstruction-pathways.csv](../data/derived/ets2-inhibitor-reconstruction-pathways.csv). The private repository link is not presented as if the recipient already has access. No files have been attached or shared externally.

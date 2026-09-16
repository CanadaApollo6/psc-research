# Independent review of the frozen PSC organoid results

**Read-only results review, 16 September 2026 UTC.** This is a new review of completed outputs, not an amendment to the frozen methods or pre-effect review. No primary/sensitivity model was changed or rerun. No new gene, pathway or program test was performed.

## Bounded conclusion

The aggregate counts are confirmed: **eight NB interaction candidates at BH q≤0.05; six also pass BY**. Four interactions are positive and four negative. All eight have matching directions in the paired log-CPM contrast, and those paired-change directions survive every complete-donor omission. However, **no gene passes the fixed-family Welch BH or BY sensitivity**. The result is model-dependent evidence for a small set of ex vivo response differences, not an independently replicated PSC-specific mechanism or clinical finding.

The numerical verification is strong within its stated scope. It does not establish small-sample P-value calibration, observational exchangeability or biological causality.

## Coverage, missingness and model diagnostics

Independent aggregation of the frozen public tables, through the project's `.venv`, confirms:

| Quantity | Result |
|---|---:|
| Complete raw features | 36,601 |
| Fixed eligible family | 17,897 |
| Excluded by the fixed count filter | 18,704 |
| All-zero features in the complete raw universe | 7,641 |
| All-positive eligible ratio-normalization basis | 17,463 |
| Finite main interaction tests | 17,894 |
| Count-boundary-unavailable interaction tests | 3 |
| Welch tests with finite P values | 17,897 |
| Primary NB interaction BH / BY candidates | 8 / 6 |
| Welch fixed-family BH / BY candidates | 0 / 0 |

The three masked features are `ENSG00000248994` / `AC124852.1`, `ENSG00000171711` / `DEFB4A`, and `ENSG00000234156` / `AL359636.2`. Each has a wholly zero **non-PSC CNT arm**, not both non-PSC arms zero. The interaction and non-PSC response remain unavailable; their correction-only entries give q=1. The supported PSC contrast remains separately available. Native finite package P values are preserved as diagnostics and must not replace the masked main fields. Excluded and unavailable rows are not negative biological findings.

The trend is **parametric**; the permitted mean fallback was not used. All 17,897 native LFC flags are true. Native and public diagnostics agree on **118 false genewise-dispersion flags and 67 false MAP-dispersion flags**, affecting **185 distinct eligible genes** with no overlap between those two flag sets. The fit-warning list is empty, but this does not erase the recorded flags. All eight BH candidates have both dispersion flags and the LFC flag true. Their own flags do not independently validate the shared dispersion model.

All Cook entries are finite. No eligible gene has a Cook value above the diagnostic cutoff **F(10,6)₀.₉₉≈7.8741**. Automatic Cook exclusions are zero; the design has no eligible replicated identical rows. “No cutoff exceedance” is not proof of no donor influence or good small-sample calibration.

Zero-total individual donor pairs occur for 97 eligible features in non-PSC and 81 in PSC; these feature counts can overlap. They remain flagged without a new filter or donor deletion. None of the eight BH candidates has such a pair.

## Complete primary candidate set and sensitivities

IDs and labels below are preserved from the source; no alias or functional reconciliation is implied. `γ` is the NB **log2 ratio of treatment ratios**. Its interval is a pointwise model-based 95% Wald interval, not a simultaneous discovery interval. `Δ` is the distinct equal-donor log2(CPM+1) change contrast. Welch P values below are unadjusted; **all eight Welch BH q values are approximately 0.9998, and their BY q values are 1**.

| Source ID / label | NB γ [95% CI] | NB BH q | NB BY q | Welch Δ | Welch P | Allocation tail |
|---|---:|---:|---:|---:|---:|---:|
| ENSG00000138061 / CYP1B1 | -1.550 [-2.110, -0.990] | 0.00106 | 0.01099 | -1.341 | 0.01725 | 4/70 |
| ENSG00000233954 / UQCRHL | -1.279 [-1.759, -0.798] | 0.001359 | 0.0141 | -1.058 | 0.003477 | 2/70 |
| ENSG00000140465 / CYP1A1 | -1.992 [-2.747, -1.238] | 0.001359 | 0.0141 | -0.657 | 0.07228 | 8/70 |
| ENSG00000231767 / AL136454.1 | -1.371 [-1.910, -0.831] | 0.002875 | 0.02981 | -1.132 | 0.01146 | 2/70 |
| ENSG00000173110 / HSPA6 | +5.576 [+3.344, +7.808] | 0.002925 | 0.03034 | +0.581 | 0.01885 | 2/70 |
| ENSG00000273132 / AL355312.4 | +0.998 [+0.598, +1.397] | 0.002925 | 0.03034 | +0.888 | 0.02406 | 2/70 |
| ENSG00000112303 / VNN2 | +1.383 [+0.799, +1.967] | 0.008896 | 0.09225 | +1.000 | 0.05246 | 4/70 |
| ENSG00000128016 / ZFP36 | +0.983 [+0.531, +1.434] | 0.04475 | 0.464 | +0.985 | 0.03743 | 4/70 |

The six BY candidates are the first six rows. **VNN2 and ZFP36 do not pass BY.** CYP1A1 and VNN2 have Welch 95% intervals that include zero. The absence of a corrected Welch discovery is a substantive qualification, not a reason to change the family or recorrect only the eight NB-selected genes.

All eight candidate paired-change effects retain the full-panel sign in **8/8** donor omissions, with no zero/opposite leaveout estimate. This is a useful influence check on the paired-change estimand. **There were no leave-one-donor-out NB refits**, so this does not show that their NB q values, dispersion estimates or NB significance survive omission. The overlapping leaveouts are not replications or confidence intervals.

Only four candidates have the minimum allocation fraction 2/70; three have 4/70 and CYP1A1 has 8/70. Across the full eligible family, 299 features attain 2/70. These are the already frozen, unstudentized, tie-aware **descriptive allocation tails**, not additional FDR discoveries or randomized-disease P values. Exhaustive enumeration does not establish exchangeability of clinically different donors.

Selected-candidate sign agreement is not universal agreement between methods. Among the 17,894 finite NB interactions, 16,436 have the same nonzero sign as the log-CPM contrast and **1,458 have the opposite sign**. The two estimands use different normalization, weighting and low-expression behavior. This descriptive difference is not, by itself, a computational error.

## Interpret the response scale correctly

A positive γ means that the modeled PSC treatment/control ratio exceeds the non-PSC treatment/control ratio. It does **not** mean a positive absolute RNA-count difference or necessarily induction in every PSC donor. Illustrations from the already fitted candidate contrasts:

- **VNN2:** PSC log2 response +3.293 and non-PSC +1.910; both fitted group responses are positive, with a larger ratio in PSC.
- **ZFP36:** PSC +0.391 and non-PSC −0.592; the positive interaction combines opposite fitted response directions. Its nonsignificant PSC within-group test (BH q≈0.224) does not invalidate or replace the direct interaction test.
- **CYP1A1:** PSC −2.163 and non-PSC −0.170; the negative interaction describes a more negative modeled response in PSC. A nonsignificant non-PSC response does not establish no response.
- **HSPA6:** γ≈+5.576 is a large ratio-of-ratios estimate, while the paired log-CPM contrast is only +0.581 and baseMean is about 18.5 normalized counts. Do not translate that model-scale interaction into a large absolute RNA gain or a uniform cell-level effect. The different scales and weighting must remain explicit.

The within-group tables contain **615 PSC versus 164 non-PSC BH responses**, and 298 versus 87 at BY. Their unequal list lengths do not test a disease-by-treatment difference and do not imply four times the biological responsiveness. The primary interaction family is the relevant comparison.

## What the independent verification does and does not establish

The completed [verification receipt](psc-organoid-independent-verification.json) reports **397 passed checks**. This results review independently checked the seven public payload hashes/byte lengths against their manifest and reconciled private native convergence flags. The root reports that the eight primary files, including the manifest, replay byte-identically; the replay was not repeated during this review.

All **16 prospectively selected conditional fixed-dispersion fits** meet the frozen numerical agreement criteria. Maximum contrast difference is approximately **8.18×10⁻⁵ native SE**, versus the 0.02 tolerance; maximum twice-objective improvement is **1.48×10⁻⁸**, versus 0.001. However:

- **Nine of 16 optimizer-success flags are false.** All 16 pass the separate prespecified stationarity criterion; maximum Newton decrement is **5.91×10⁻⁷**, below 1×10⁻⁵. Report this distinction rather than saying every independent optimizer returned success.
- The deterministic panel contains **none of the eight BH candidates**. Complete Wald arithmetic was checked for all eligible features, including those candidates; independent likelihood optimization was only the prespecified 16-feature panel.
- Dispersion estimation, its parametric trend, prior and MAP shrinkage were **not independently refitted**. Agreement is conditional on supplied native dispersions. Numerical replay is not a second biological cohort or proof of asymptotic P-value validity.

## Reporting boundary

The independent units remain **eight donors**, not 16 libraries, 46,343 cells, 17,897 genes or eight leaveouts. Controls are procedure controls, not healthy volunteers. Donor blocking does not remove disease-associated age, source, clinical-history, sampling or culture-selection effects on treatment response. Relative-RNA measurements and changes in the retained cell/state mixture are not absolute transcription rates or patient outcomes.

The raw full-feature universe and its **70,263-count difference** from historical RNA metadata remain the frozen measurement scope. No forced reconciliation occurred. Numerical verification cannot recover an undocumented historical feature filter or assay-export command.

A defensible summary is: **“The frozen donor-aware count model identified eight candidate disease-associated IL-17 response interactions, six under BY. Their paired-change directions were stable under donor omission, but no feature passed the complete-family Welch sensitivity. The findings require independent biological validation.”** Do not claim equivalence for the other genes, universal absence of an IL-17 response, a new established PSC mechanism, or benefit from IL-17 blockade. Published organoid IL-17 response differences were already known before this reanalysis.

### Provenance of this read-only review

- Execution freeze: `config/psc-organoid-il17-execution-freeze.json`, SHA-256 `89650c5f72875cd26913ada15f6874194b4f4694fd05b182418b1bf00d44a58b`.
- Public input tables: `data/derived/psc-organoid-il17/`; manifest SHA-256 `43c382dcf1ae6d11adac1e07a53a59c8a5000f4e617985fb222018328a0ce88a`.
- Verification receipt SHA-256: `30852357b68254ae3a1e93ab4f6bdfdc18463db87d9c923eb20395981536d34a`.
- Private reads were limited to native feature flags and the warning list in `work/psc-organoid-il17-primary/`. No donor-level expression values are reproduced here.
- Aggregation used `.venv/bin/python`, pandas 3.0.5 and NumPy 2.5.3 with literal-string source reads and explicit numeric conversion. Only this new report was written. No tests were rerun because this task changed no code, parser or frozen file.

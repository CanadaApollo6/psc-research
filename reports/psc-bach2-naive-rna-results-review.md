# BACH2 naive-CD4 RNA: independent aggregate results review

**16 September 2026. Accepted for complete aggregate reporting. Close this fixed branch.** No aggregate inconsistency or interpretation blocker was found. This is one same-cohort observational RNA check, not a nominal discovery, an equivalence result or independent biological replication.

## Scope and verification

This post-execution review read only the hash-pinned public aggregate, frozen metadata/protocol and aggregate completion, independent-audit and replay receipts. It did not open source matrices, archives, donor joins, private donor quantities or donor-indexed omissions. It computed no new statistic or endpoint and changed no frozen method, code or schema. The earlier pre-effect review and its access disclosure remain unchanged.

Root records checkpoint `a1b0a6e1d880f9b97c49bdbc0358aef7db408a47` before real endpoint access. The pinned independent-audit receipt reports a passing reconstruction from all eight complete MatrixMarket streams, all eight donor records and all eight whole-donor omissions. This review checked that receipt's bindings; it did not repeat that source audit. The replay receipt confirms byte-identical private scientific measurements. The public aggregate differs only in the predeclared authorization digest. Completion and destination-binding differences are operational provenance, not scientific changes. Neither audit nor replay adds biological observations.

The frozen verification protocol used 4,000-digit independent direct-ln/exact-rational reference calculations. Nonzero point, moment, SD, SE, df and P agreement required relative error within `1e-10`; zero, sign, constancy and availability required exact agreement, without an absolute epsilon. Standardized CI critical values required relative agreement within `1e-9`. The P reference used a 100-digit Decimal incomplete-beta series with float log-gamma normalization; CI inversion used `scipy.special.betainc` and `scipy.optimize.brentq`. Shared libraries and finite-precision components remain limits to software independence. Working digits and numerical agreement are not biological precision. No post-value relaxation is reported or used here.

The native aggregate checker passed 34 checks of approved input hashes, public schema, displayed numbers, complete diagnostics, provenance and the audit/replay receipts. Its receipt and log are under `work/psc-bach2-naive-rna-independent-review/post-execution-review/`.

## Complete fixed result

The endpoint remains the unique released `ENSG00000112182 / BACH2 / Gene Expression` feature in pooled source C0 `CD4+ TN RTE` and C1 `CD4+ TN mature` cells. Each donor contributes:

`Y = log2(1 + 1,000,000 × selected-cell BACH2 RNA units / selected-cell all-RNA units)`.

Donors have equal weight. The denominator includes all 36,601 released RNA rows in the same selected cells and excludes all 38 ADT rows. Values below are rounded displays of the frozen output, not replacement estimates.

| Source group | Donors | Mean donor Y | Donor sample SD |
| --- | ---: | ---: | ---: |
| SNP | 4 | 5.26602480 | 0.416728747 |
| noSNP | 4 | 5.71012493 | 0.313682619 |

| Sole SNP-minus-noSNP contrast | Result |
| --- | ---: |
| Point difference | −0.444100136 |
| Standard error | 0.260796680 |
| Actual Welch df | 5.57342779 |
| Nominal two-sided 95% CI | [−1.09427507, 0.206074803] |
| Nominal two-sided P | 0.143269126 |
| Availability | `available_nominal_Welch` |

**The point estimate is negative, but the nominal interval spans zero and both directions.** The result does not meet the nominal two-sided 0.05 threshold. This is not evidence of RNA equivalence, unchanged RNA or a biologically negligible effect. No equivalence margin was specified. The contrast is on the stated donor-average transformed-relative-RNA scale; it is not a raw RNA fold change, absolute RNA amount, transcription rate or protein difference.

### All fixed sensitivity points

All **eight** whole-donor omissions were included. Their complete point range is **[−0.631587772, −0.297536419]**. There are **8 negative, 0 positive and 0 zero** omission points. The full-point sign is negative, with **0 strict sign reversals**.

The negative point direction survives each single-donor omission. This is a point-direction diagnostic only. It does not establish population-level direction, precise magnitude, absence of confounding or freedom from multi-donor influence. No omission P value or CI was calculated. The overlapping omissions are not independent replicates, and no donor was removed from the primary result.

### Complete public QC

- All **8** fixed donors are represented; the source and membership gates passed.
- **16,659** selected cells come from the fixed **55,460** retained rows. Frozen metadata split these into **9,320 C0** and **7,339 C1** cells. No source-only cells were added.
- The denominator uses **36,601 RNA rows**; all **38 ADT rows** are excluded.
- **0** donors have zero pooled BACH2 units. The rule to retain zero-B donors was not changed.
- **2,408** selected cells have positive released BACH2 units, cohort-wide. This is a fixed assay diagnostic, not a genotype-group detection endpoint or a basis for new filtering.

Individual cell joins, donor B/R/Y, per-donor state and detection fractions, and indexed omissions remain private. The pinned full-output audit records agreement on the private metadata/units and all omission points. This aggregate review did not inspect those private diagnostics. Aggregate-only publication is not a formal differential-privacy guarantee.

## Interpretation limits

1. **Same small, selected cohort.** These are the original eight male PSC donor units, four per source genotype label. This is not a PSC-versus-control comparison or independent replication. Source independence and approximately Gaussian donor Y are working assumptions. Four donors per arm cannot establish calibration, exchangeability or generalizability; many cells do not increase the donor sample size. Published findings were known during design, and the historical figure's normalization and test units are not fully specified.
2. **Observed mixture and relative-RNA weighting.** Pooling preserves the observed C0/C1 mixture and RNA-content weighting. It cannot separate state redistribution, within-state regulation or other compositional effects. Equal donor weighting does not remove these within-donor features. No cell-intrinsic RNA mechanism is established.
3. **Observational labels, not causal alleles.** Genotype-group labels do not identify a causal nucleotide. LD, ancestry, batch, clinical selection, treatment/severity differences and undocumented retention history remain possible explanations. No verified participant-to-archive clinical row key was available. No clinical covariates or alleles were imputed or adjusted here.
4. **Released units, not certified original UMIs.** The source export/version/reference history remains unresolved. The original-UMI/count-likelihood gate is still unpassed. Passing source integrity and numerical checks does not resolve that history or validate a count likelihood.
5. **RNA does not settle protein or translation.** Neither this negative point nor the nonsignificant nominal result proves or refutes the separate BACH2 protein/miR4464 translation experiments. There is no mature-miR4464 or BACH2 ADT endpoint here. No causal PSC mechanism, clinical benefit or treatment recommendation follows.

## Decision and audit trail

**Report the full aggregate result and close this source-cohort branch.** Do not select a new gene/state, alter normalization, drop donors, add a test or search another source to rescue the outcome. Closure follows the frozen stop rule regardless of P, interval width, direction or sensitivity. No additional analysis is proposed or authorized by this review.

Authoritative approved inputs:

- Root reporting handoff: `work/psc-bach2-naive-rna-freeze/root-reporting-handoff.json`, SHA-256 `53f636477b1f8c33c44308b3b48ab3f3a9b448f53422dc17f1e07fbdd935fab6`.
- Public aggregate: `work/psc-bach2-naive-rna-plan/runs/primary-v1/public-aggregate.json`, SHA-256 `3f9383678afc4a576c12f89225e3148ca07ff65486c8187e8579f4dfd4ac9837`.
- Independent source/numerical audit: `work/psc-bach2-naive-rna-execution-review/primary-v1-audit/aggregate-verification.json`, SHA-256 `f8d70ba6bfa6abef270697c6f7e809b77a5360198dfeaec0c46767f763eacf65`.
- Root replay verification: `work/psc-bach2-naive-rna-freeze/root-mechanical-replay-verification.json`, SHA-256 `2d6a7c4718929ed29bf62ce00c2599e984491b01c33cbdc13c1b6bf4d9a2ffa4`.
- Frozen scientific and verification protocols: [`docs/psc-bach2-naive-rna-protocol.md`](../docs/psc-bach2-naive-rna-protocol.md) and [`config/psc-bach2-naive-rna-verification-protocol.json`](../config/psc-bach2-naive-rna-verification-protocol.json).

All approved read pins, the native aggregate-check receipt and this review's completion receipt are in `work/psc-bach2-naive-rna-independent-review/post-execution-review/`. They distinguish receipt-level interpretation review from the separately completed real-source reconstruction. No Git, network, installation or new data acquisition occurred in this review.

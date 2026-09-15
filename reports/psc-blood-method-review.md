# PSC blood method review

Reviewed 15 September 2026. The statistical/code concerns were sent at **22:48:20 UTC**, before the [22:56:48 execution freeze](../config/psc-blood-execution-freeze.json) and new effects. This note was written after permission to inspect aggregate results. It changes no frozen method, threshold or result. No new genes or sources were searched.

## Pre-effect concerns and dataset-scoped mitigation

1. **Error-control interpretation.** Direction-consistent maximum two-sided P is a conservative conjunction construction. It does not establish the cross-gene dependence assumptions needed by BH. The [protocol](../docs/psc-blood-replication-protocol.md) clarified this **before execution** and prespecified BY for all three families. Its current hash matches the execution freeze. BH findings remain assumption-dependent; BY permits arbitrary dependence given valid constituent P values. A fixed low-count filter does not itself prove conditional P-value calibration. Secondary disease-control families remain separate exploratory claims, not one omnibus error-control guarantee.
2. **Identity/audit edge case.** The integration function rejects duplicate RNA stable IDs before coverage accounting. A synthetic eligible ID plus two excluded blank IDs caused a failure; version-colliding IDs could likewise halt integration. The RNA producer explicitly represents such problematic identities. **This generic edge case was not repaired.** For this dataset, the [RNA QC](../data/derived/psc-blood-rnaseq/rnaseq-qc.json) records **63,677 unique versionless IDs**, with no invalid or colliding rows. The edge case is therefore absent here. Future datasets still need explicit source-feature identity accounting; they must not select among ambiguous rows by P value.
3. **Input-pin completeness.** The integration runner verifies supplied hash entries but does not require every consumed input to be listed. **This general guard was not changed.** The actual [mechanical integration plan](../config/psc-blood-integration-plan.json) explicitly includes both result files, the full crosswalk and the executing integration script. This review verified all four hashes. Its four contrast values are distinct; aliases select the raw primary Wald/Welch P values and PSC-minus-reference effects, not adjusted P values or sensitivity tests. The integration receipt was created after assay outputs existed; the conjunction methods/code were already frozen before effects.

The nine integration tests passed in the native project environment during the pre-effect review. These checks and the dataset-specific mitigations do not claim a universally hardened wrapper or independent biological validation.

## Aggregate-result peer check

The [complete summary](../data/derived/psc-blood-replication/summary.json) contains **11,074 shared eligible genes**. Three genes retain unavailable constituent tests and correction input 1; they were not removed from the family.

| Prespecified conjunction | BH q ≤ 0.05 | BY q ≤ 0.05 |
|---|---:|---:|
| Same-direction Norwegian and Polish PSC/control | 93 | 0 |
| Additionally same-direction Polish PSC/UC | 0 | 0 |
| Additionally same-direction Polish PSC/PBC and PSC/CD | 0 | 0 |

The primary list has 84 increases and nine decreases. None of the seven fixed benchmarks passes the primary criterion. ETS2 has BH q **0.0667181376321612**; it remains a failed criterion, not a near-positive or reason to change the threshold.

Among the 93 BH passers, one row has `rna_genewise_converged=False` and one has `rna_MAP_converged=False`; none has `rna_LFC_converged=False`. These counts need not represent different genes. Their names were not inspected or nominated in this review. Keep the frozen result and visible flags; do not describe every passer as fully numerically qualified or silently remove flagged rows after seeing effects.

**Supported wording:** “Ninety-three genes met the prespecified same-direction two-cohort BH conjunction criterion. None passed the BY sensitivity or either stronger disease-control conjunction.” This is a dependence-sensitive candidate list of observational whole-blood associations, not 93 causal PSC genes, a validated diagnostic or proof of intrinsic PSC specificity.

Zero stronger-conjunction discoveries means **no additional support under those criteria**. It does not establish that the primary associations are biologically nonspecific, that disease groups are equivalent, or that no distinguishing gene exists. Limited power, stringent same-direction requirements and the Polish age/sex/treatment confounding remain relevant. Separate study families are not identity-verified disjoint participants, and prior published use of the Polish dataset prevents calling this an untouched validation.

This review is not a negative-binomial refit or the separate whole-pipeline numerical verification. It does not certify completion of that verification; its status is maintained in the [independent verification record](psc-blood-independent-verification.json).

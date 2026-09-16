# GSE303271 fixed-program identifier preparation

## Status and boundary

Completed offline **identifier-only preparation**, not an analysis freeze. Root must choose the identity, normalization and matching rules before effects. No normalization, expression matching, scoring, effects or models were run. No network requests, aliases, ambiguity collapse, zero-filling or organoid-derived genes were used.

All **41,096 released RNA rows** remain in original order. Candidate identifier pools are **not normalization universes**. No abundance, detection or protein-coding filter was applied. Detailed labels, source memberships, ambiguous candidates and indexes remain under `work/psc-liver-program-inputs/`; this report and the public JSON contain aggregates.

## Source lineage

The complete original `source_value` memberships were recovered from the pinned `ets2-program-membership.csv`, including previously unmapped and ambiguous atlas members. SC/SN copies agree in order; they are duplicate exports, not independent biological sources. The four original Ensembl lists also agree with the complete source-gene-set export. Old mapped fields, MacroMap/atlas subsets and legacy source-statistic columns were not used to define programs.

There are **4,590 original memberships** across eight programs and **4,589 intended memberships** after the single fixed ETS2 exclusion. The original ETS2 ZIP remains absent; this is verified-export reuse, not a new archive read.

The full Ensembl 116 / June 2026 reference has **91,748 relationships and 86,411 Ensembl IDs**. All relationships are used before any feature/program restriction, including 54,919 Ensembl IDs without Entrez. The reference gzip SHA-256 is `9db59e5c3ccf1c7f70cbc98aa59e08c0b1168e15ce16ee91c8158e3816cf1aab`; the historical normalization verification independently pins it.

The absent historical BioMart response was reconstructed exactly from preserved original fields/order, its pinned header, tabs/LF and terminal `[success]`. The owned reconstructed copy is **2,236,970 bytes**, SHA-256 `a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167`. This is not a fresh server retrieval or a reread/restoration of the absent original cache. The count response also reconstructs to its historical pin. Success/count agreement supports the returned export but cannot rule out an unobserved remote backend omission.

Historical source URLs and retrieval times:

- RNA count file: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303271/suppl/GSE303271_raw_counts.txt.gz — 2026-09-16T00:24:24.155938+00:00.
- BioMart: https://jun2026.archive.ensembl.org/biomart/martservice — 2026-09-12T17:40:41Z; independent count response 2026-09-12T17:40:57Z.
- ETS2 archive: https://zenodo.org/api/records/10707942/files/JamesLeeLab/ETS2manuscript_Stankey_CT_et_al_2024-v1.0.zip/content — 2026-09-12T14:21:34.543387+00:00 (historical manifest only; ZIP absent).

Full source/member/export pins, query XML and reconstruction receipts are in the JSON and `source-audit/findings.json`. Source previews incidentally displayed existing legacy statistic strings; these were not converted, summarized, selected on or used for conclusions. Preparation/audit code retains only identifier/provenance fields. No new liver expression quantities were inspected.

## Two prospective identity interfaces

Both policies were recorded in `work/psc-liver-program-inputs/scope.json` before mapping. Neither was selected by these coverage results:

1. **Core** (`exact_symbol_ensembl_bijection`): exact case-sensitive submitted symbol ↔ one current-reference Ensembl ID, plus one submitted row for that ID. Entrez is reported, not required as an RNA bridge.
2. **Strict** (`plus_reciprocal_entrez`): additionally require globally reciprocal one-to-one Ensembl ↔ Entrez relations. This is a separate optional cross-namespace interface, not the chosen analysis policy.

All released labels are gene-name-style; none are literal Ensembl IDs. Ambiguity checks use the complete human reference, not just measured genes or programs. The reference has no coordinates, biotypes or assembly-region fields. Multiple Ensembl links do not by themselves prove separate biological loci; no primary-locus choice was made.

The table shows mutually exclusive eligibility statuses, not all annotation relationships. Entrez conditions do not exclude core rows; the core zeros below do not mean those relations are absent.

| Feature status | Core | Strict |
|---|---:|---:|
| `eligible` | 37,778 | 22,562 |
| `no_exact_reference_symbol` | 209 | 209 |
| `symbol_has_multiple_ensembl` | 3,108 | 3,108 |
| `ensembl_has_multiple_symbols` | 1 | 1 |
| `no_entrez_relation` | 0 | 14,781 |
| `multiple_entrez_for_ensembl` | 0 | 328 |
| `entrez_has_multiple_ensembl` | 0 | 107 |

Missing Entrez is not absence of an RNA row or biological absence. Missing current symbols are not automatically retired genes. Every unresolved relationship remains explicit; no alias, zero substitute or arbitrary ambiguous ID was used.

## Original denominators and coverage

The unchanged mapping gate is at least **80% of the original intended members and ten mapped genes**. All eight programs pass that annotation gate under both interfaces. This does not establish matching-pool sufficiency or authorize scoring.

| Fixed program | Original intended | Core mapped | Strict mapped |
|---|---:|---:|---:|
| `ets2_g1_dn` | 927 | 848 (91.48%) | 825 (89.00%) |
| `ets2_g1_up` | 668 | 612 (91.62%) | 536 (80.24%) |
| `ets2_g2_dn` | 875 | 809 (92.46%) | 789 (90.17%) |
| `chr21_dn` | 141 | 126 (89.36%) | 122 (86.52%) |
| `inflammation` | 815 | 750 (92.02%) | 744 (91.29%) |
| `interferon_gamma` | 139 | 121 (87.05%) | 118 (84.89%) |
| `oxidative_stress` | 436 | 403 (92.43%) | 400 (91.74%) |
| `apoptotic_signaling` | 588 | 544 (92.52%) | 542 (92.18%) |

The strict opposite-direction program has 536/668 members (80.24%); its original denominator is not reduced. Per-program missing-symbol, missing-current-Ensembl, absent-row, multiple-symbol and Entrez failure counts are separate in the JSON and complete ignored audit.

**Nonoverlap sensitivity:** subtract the mapped union of inflammation, interferon-gamma, oxidative-stress and apoptotic-signaling programs. Core retains **695** of 848 mapped primary genes after removing **153**; strict retains **673** of 825 after removing **152**. The gate remains the parent's 927-member coverage gate plus the separate ten-retained-gene minimum. The retained fractions 695/927 and 673/927 do not create a new denominator or a new 80% test. Both sensitivities pass only this annotation gate.

## Complete candidate matching-reference pools

| Interface | ID-qualified RNA rows | Eight-program union excluded | ETS2 excluded | Background candidates |
|---|---:|---:|---:|---:|
| Core | 37,778 | 2,944 | 1 | 34,833 |
| Strict | 22,562 | 2,831 | 1 | 19,730 |

Every source row appears in each pool audit, including identity-ineligible rows. All eight mapped original programs are excluded, regardless of their coverage gates. These are complete ID-only candidate lists, not expression-matched sets. No abundance bins/draws, pool-sufficiency result, normalization, scores or effects exist. Any later normalization retains the separately qualified full 41,096-row count universe; these lists must not replace it.

## Remaining limits

- The article declares GRCh38, STAR 2.7.2b and FeatureCounts 1.6.4. Exact historical GTF/provider/release/checksum and counting options remain unavailable. Current-reference symbol identity does not reconstruct that GTF or establish historical locus equivalence. Allele conventions and inference-model versions are not applicable here.
- Preserve **17 PSC, 17 ASC and 30 AIH** as released; ASC is not folded into PSC. This is whole-biopsy bulk liver, not sorted macrophages or cholangiocytes. The article asserts 64 patients, but an independent person/visit crosswalk and sample-linked clinical/stage/batch covariates remain unavailable.
- GSE159676 remains **HOLD** for expression: unreconciled mixed-scale literals and **six paired PSC patients**, not twelve independent patients. No repair or adult expression use occurred.
- NoPSC atlas remains **PARK**: no qualified complete count/feature/barcode/donor export. No atlas subset was used.
- External program coverage is not evidence for ETS2 specificity, causal mechanism, transfer, replication or clinical utility. Root owns the next method choice and prospective analysis freeze.

## Reproduce and verify

```text
.venv/bin/python scripts/prepare_psc_liver_program_inputs.py --prepare
.venv/bin/python -m unittest discover -s tests -p test_psc_liver_program_inputs.py -v
.venv/bin/python scripts/prepare_psc_liver_program_inputs.py --verify
```

`--prepare` writes only the four owned outputs/ignored work area; it does not acquire data. `--verify` checks pins read-only. The JSON schema is version 1, with source lineage, both interfaces, full feature-status counts, original program coverage/status denominators, nonoverlap gates, candidate pools, transformation/scope boundaries, preservation pins and validation receipts.

Independent source audit: **4,750 checks and 7 tests passed**, with byte-identical reconstruction replay. Closed liver/organoid versioned files and prior source/artifact pins are checked unchanged.

Final validation receipt:

```json
{
  "checked_artifacts": [
    {
      "bytes": 4771,
      "path": "work/psc-liver-program-inputs/identifier-tests-final.log",
      "sha256": "a037c6beb6bdc269e418d1b57829164a9b6663238b18e83ecfd3594c2c022e1d"
    },
    {
      "bytes": 106412,
      "path": "work/psc-liver-program-inputs/repository-tests.log",
      "sha256": "5cb0079f6c9c302fec86ce466a4cea8e21968cdde42cebae6f5995b160d6d239"
    },
    {
      "bytes": 157,
      "path": "work/psc-liver-program-inputs/closed-liver-readonly-verify.log",
      "sha256": "13475d53799d85b1e375a34547cdddb0b3b5cb62630a9b37e50c1ccb364b4145"
    },
    {
      "bytes": 3312,
      "path": "work/psc-liver-program-inputs/source-audit/verification.json",
      "sha256": "d7dc95c645f19216a1375aa0c9f70147e79daecc7d540446a818995271891ffd"
    },
    {
      "bytes": 49919,
      "path": "work/psc-liver-program-inputs/source-audit/findings.json",
      "sha256": "a11fae6df4a3fb42ce547d95605be4782ba5dab3d17317b8e4762d29f5eaf1fb"
    },
    {
      "bytes": 30358,
      "path": "work/psc-liver-program-inputs/source-audit/run_source_audit.py",
      "sha256": "abb9443103395043d1f16ec2ff2757ffd20300d707383b3811f2573a8ce69e1a"
    },
    {
      "bytes": 9463,
      "path": "work/psc-liver-program-inputs/source-audit/test_source_audit.py",
      "sha256": "8edbfae521985e44d1681dab56ae45046509832bc04488f1537f8e8b57097ac4"
    },
    {
      "bytes": 39639,
      "path": "work/psc-liver-program-inputs/independent-audit/independent-verification.json",
      "sha256": "63a402a4639d45f324044468e9add7269f89d213d80baa15085496fe86f2f093"
    },
    {
      "bytes": 37531,
      "path": "work/psc-liver-program-inputs/independent-audit/verify_identifiers.py",
      "sha256": "7fd9e758961a19e44bfbe37348e73f55dcccb3cc3cc2e33d330962bdf0b6a54b"
    },
    {
      "bytes": 716,
      "path": "work/psc-liver-program-inputs/independent-audit/audit-artifact-hashes.json",
      "sha256": "fc9b545751d024876176712893f02f0336c3946fd913dc08e391502dcd792dc7"
    },
    {
      "bytes": 934,
      "path": "work/psc-liver-program-inputs/context-snapshot-boundary.json",
      "sha256": "686d029f2f7053aeb18daa49bff571de334375212e9190e163e49177afde5c5c"
    },
    {
      "bytes": 3283,
      "path": "work/psc-liver-program-inputs/pre-validation-replay.json",
      "sha256": "4bbbdf0f0be0ff3be78443fb343a0cab601c37e9eb590e1e1ed4addc3d91f1c1"
    }
  ],
  "closed_liver_read_only_verifier": {
    "closed_liver_versioned_files": 4,
    "closed_organoid_versioned_files": 4,
    "command": ".venv/bin/python scripts/qualify_psc_liver_inputs.py --verify",
    "exit_code": 0,
    "liver_ignored_artifact_pins": 16,
    "liver_source_pins": 46,
    "no_closed_outputs_rebuilt": true,
    "prior_organoid_material_pins": 67
  },
  "context_snapshot_boundary": "The parent updated docs/evidence-log.md and reports/macromap-program-design.md after source-audit snapshot. These were contextual reading, not identifier inputs or frozen source tables. Original receipts are retained; this task did not edit either document.",
  "independent_identifier_audit": {
    "all_features_compared": 41096,
    "all_original_members_compared": 4590,
    "all_pool_rows_compared": 82192,
    "all_qualified_index_rows_compared": 9657,
    "checks_passed": 965516,
    "failed_checks": 0,
    "no_producer_import_or_execution": true,
    "shared_Entrez_ambiguities_hidden_by_subset_but_retained_by_whole_reference": 61
  },
  "native_python": ".venv/bin/python",
  "offline_preparation_replay": {
    "eight_ignored_preparation_artifacts_identical": true,
    "four_versioned_files_identical": true,
    "stage": "before embedding this validation receipt; the final published-byte replay is recorded in completion-receipt.json"
  },
  "original_source_audit": {
    "checks_passed": 4750,
    "direct_cached_historical_raw_reread": false,
    "exact_historical_BioMart_bytes_reconstructed": true,
    "independent_tests_passed": 7,
    "offline_result_replay_byte_identical": true
  },
  "repository_tests": {
    "command": ".venv/bin/python -m unittest discover -s tests -v",
    "error_test_labels": [
      "ERROR: setUpClass (test_ets2_inhibitor_reconstruction.TestETS2InhibitorReconstruction)",
      "ERROR: setUpClass (test_ets2_sources.TestETS2Sources)",
      "ERROR: setUpClass (test_gse84161_mapping.TestGSE84161Mapping)",
      "ERROR: setUpClass (test_gse84161_metadata.TestGSE84161Metadata)",
      "ERROR: test_check_only_preserves_all_thirty_metadata_ordered_inputs (test_gse84161_normalization.TestGSE84161NormalizationGuard.test_check_only_preserves_all_thirty_metadata_ordered_inputs)",
      "ERROR: test_execute_requires_frozen_input_gate (test_gse84161_normalization.TestGSE84161NormalizationGuard.test_execute_requires_frozen_input_gate)",
      "ERROR: test_summary_csv_and_entrez_medians_on_synthetic_inputs (test_gse84161_normalization.TestGSE84161NormalizationGuard.test_summary_csv_and_entrez_medians_on_synthetic_inputs)"
    ],
    "errors": 7,
    "exit_code": 1,
    "failures": 0,
    "missing_history_not_repaired": true,
    "note": "Repository-wide count reflects concurrent parent-owned work; the 30 new task tests all passed.",
    "tests_run": 606,
    "unrelated_missing_historical_paths": [
      "data/raw/ets2-benchmark/ets2-primary-article.xml",
      "data/raw/ets2-benchmark/ets2-public-release.zip",
      "work/ets2-independent-search/GSE84161-family-soft.txt",
      "work/gse84161-normalization/runtime/bin/Rscript"
    ]
  },
  "schema_version": 1,
  "status": "identifier_preparation_validated_with_separate_unrelated_repository_errors",
  "unit_tests": {
    "command": ".venv/bin/python -m unittest discover -s tests -p test_psc_liver_program_inputs.py -v",
    "exit_code": 0,
    "includes_synthetic_failed_program_gate_pool_exclusion": true,
    "tests_run": 30
  }
}
```

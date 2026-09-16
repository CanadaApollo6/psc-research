# BACH2 naive-CD4 RNA: independent pre-effect review

**Accepted for a separate root freeze, not for expression execution.** The final V7 method has no unresolved reviewer blocker. This review computed no real BACH2 selected-population quantity and opened no real matrix or technical-margin payload. Root must still authorize the exact execution contract separately.

This is one modest-value RNA-observation check in the original eight donors. It is not exact historical Figure S2B computation, independent biological validation, a composition-adjusted effect, a causal allele test, translation evidence or a clinical result.

## Access disclosure

The initial schema/source audit also loaded the closed cohort-qualification JSON. That document contains previously qualified donor RNA/ADT technical totals. Those totals were not extracted, compared or used to select or evaluate this endpoint. The required evidence log and qualification reports were also read. No real source matrix, technical-margin payload or BACH2 selected-population quantity was read.

The final verifier instead uses **32 individual integrity/acquisition receipts** and annotation files. Its guard explicitly blocks the quantity-bearing qualification JSON. This stricter final check does not erase or relabel the earlier access. The detailed disclosure is `work/psc-bach2-naive-rna-independent-review/closed-receipt-access-disclosure.json`.

## Independent source and selection audit: passed

The verifier imports no production helper for this audit. It independently reconstructs the SDRF sample/donor/group mapping, preserves repeated file columns and checks every retained sample-plus-barcode key. Only the exact matching sample-name prefix can be removed. Bare barcodes are never joined globally.

- Eight distinct accession-local donor labels: four `SNP`, four `noSNP`, all male PSC cases.
- **55,460 ordered retained keys** match the compiled plan. All **16,659 selected cells and six associated fields** match independently: C0 9,320 and C1 7,339. Both states occur in every donor. All 3,855 source-only cells remain excluded.
- Eight complete identical feature axes contain **36,601 RNA and 38 ADT rows**. The unique exact `ENSG00000112182` / `BACH2` / `Gene Expression` triple is at **one-based position 12,314**. All RNA denominator positions match; ADT stays separate.
- The full cached Ensembl 116 export confirms the reciprocal exact BACH2/Ensembl identity across 91,748 relationships and 86,411 IDs. Its SHA-256 is `a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167`. These are previously reconstructed exact historical BioMart bytes, retrieved 2026-09-12T17:40:41Z, not a new retrieval or source-build/variant reconciliation.
- All eight archives were hashed as **opaque bytes only**, as allowed by root. No archive member was decoded here. The final matrix pins are compared to individual integrity receipts; dimensions come from features/barcodes. Fifty-five annotation/receipt hashes in the compiled manifest pass.

No gene label, allele, dosage, genomic build, clinical-row link, activation subcluster or retention reason was invented. The historical exporter/reference uncertainty remains.

## Mathematical review: passed

For each donor, sum BACH2 source units B and all 36,601 RNA-feature units R over the **same fixed C0/C1 cells**, then compute Y=log2(1+10^6 B/R). The sole contrast is the equally weighted four-SNP-donor mean minus the four-noSNP-donor mean. No other gene/state endpoint, detection filter or alternate normalization is introduced.

With a=s_SNP²/4 and b=s_noSNP²/4, using donor sample variances, SE=sqrt(a+b) and df=(a+b)²/(a²/3+b²/3). The 95% interval uses the actual Student-t df, not a normal approximation. The independent oracle uses exact rational centered moments and the incomplete-beta two-sided tail, with a separately inverted critical value. It does not import production mathematics or `scipy.stats` Welch/t helpers. Decimal and SciPy remain shared numerical dependencies; this is not an independent library certification.

Production evaluates group means, the contrast and all eight omissions through the algebraically identical exact product of Q=1+10^6 B/R. A stable near-one logarithm prevents spurious cancellation signs. Variance remains the exact rational pairwise sample variance of represented donor Y values; no variance threshold is used. A separate **1,100-digit direct-log oracle** verifies exact product equality with positive variance and a genuine fourth-order difference near 10^-371. Neither true zeros nor genuine tiny differences are silently changed. Working precision does not imply biological precision or high-precision Student-t calibration.

The required boundary behavior passes:

- B=0 is retained as Y=0 when R>0. Missing donors, altered membership or invalid R block the whole endpoint.
- One constant arm is valid when SE>0; its actual df is three. Both constant arms give a descriptive point only, with CI/P unavailable. Nonfinite inference is not forced to P=0/P=1. Numeric tail underflow leaves P unavailable.
- All eight whole-donor omissions remain point contrasts only. There is no omission P/CI, donor deletion, cell bootstrap, permutation claim, count likelihood or normality-test switch.

The model still assumes independent donors and approximately Gaussian donor Y within groups. Four per group cannot establish calibration. The pooled C0/C1 mixture, RNA-content weighting, batch/clinical/selection/LD explanations and original-cohort reuse remain limitations. No equivalence margin is supplied.

## Gates, privacy and completion: passed

The actual default CLI passes under a guard with **zero study-matrix/quantity-JSON/margin open attempts and zero file writes**. Plan directories must have the preparation role before any JSON open. Authorization files must have the separate root-freeze role, a simple JSON filename, a bounded size and no symlink path before their body is opened. Explicit root authorization, an external authorization hash and exact plan/code/runtime/source/selection/schema bindings are required before the matrix-open site.

Fresh ignored destinations reject reuse, aliases, traversal, arbitrary nesting and symlink components. An unsuccessful run remains consumed. Public output is a recursively checked allowlist: group aggregates, the sole contrast, all-omission count/range/sign summary and aggregate QC/provenance. Donor quantities, joins, state fractions and donor-indexed omissions remain private. No individual-value plot is produced or approved.

A successful native run must contain exactly:

```text
measurements.private.json
public-aggregate.json
authorization.used.json
completion.json
```

`completion.json` is published last through an atomic rename. It binds the three payload byte counts/SHA-256 values, output directory and all provenance digests. Missing, extra, modified, linked or incorrectly bound files fail validation. Injected failures before or after completion remove the marker and leave a failure record. A completion inventory authenticates bytes and bindings, not mathematical or biological truth; the separate execution auditor remains required before interpreting real results.

## Native validation and replay

**126 bounded native tests pass: 80 author tests plus 46 independent tests.** They cover the endpoint mathematics, full RNA-only selected-cell extraction, sparse zeros, exact integer parsing, identity/membership failures, root-authority bindings, pre-open file roles, runtime tampering, destination safety, privacy and completion failures. No full-suite success or catalog rebuild is claimed by this review; import logic for existing research tables was not changed.

Two isolated native processes execute the production extraction, analysis, projection and output writer on invented sources only. They use the full 36,601-RNA/38-ADT structure, 16,659 invented selected cells and eight invented donors, including zero B and distracting ADT/unselected target units. Only source-plan injection is mocked. All **four native artifacts plus the reviewer receipt replay byte for byte**. The largest ordinary synthetic numerical difference from the independent oracle is 1.33e-15, in an interval bound. No synthetic authorization can be used with the real project plan.

Preparation repairs remain visible in earlier logs: a reviewer quantile-literal correction; strict privacy validation after canary failures; exact-product zero/sign handling; stronger installed-file hashing; pre-open authorization/plan roles; and the completion inventory. A filename-only reviewer guard initially matched package fixtures during opaque runtime hashing and was corrected to distinguish installed software from study inputs. Stale draft pins were rejected. None changed source values, membership or the scientific estimand.

## Exact acceptance pins

| Artifact/contract | SHA-256 |
|---|---|
| Production script | `ba7ca52b0ed357e9a1c245f99aefaad8b4d7ac107ccf04b5ef1c4664b399b1ec` |
| Production tests | `7263c3cd5c8d656647bf85ae87abe19a0e7e7cb9569445dae46240ba07145e1b` |
| `candidate-v4/plan.json` | `b6826aa5adb6a3ff26c6f81f5fb16017df5e9580487a96afc2ff491554257b06` |
| Compiled private selection | `63e01409919ae4bd7ab016417a6f285666395464a31a9c53d0129bac381a6257` |
| Source manifest | `10ddbcbbf9ab9b1a413361da403d934e63e5353f3b7dfb3635e7af06aabf7d26` |
| Runtime | `4e824ecf01cba865acc06363c49274e4a594e30fbf32d53606a77ca7d97b0d30` |
| Public schema | `b1a305315aeba1080b8c15ccf33b6d86641d5071e51d7e5d72bed6a983cc97af` |
| Independent source/math verifier | `15fcd8d2926a0bb83c076aa707c5fa50e175763f905039b90ed2e2802e73a18d` |
| Independent tests | `87f78b012ca6c5c4cd4bae0a42c951bc282a349c9feae60fe54869617e39d2d1` |

The plan directory is `work/psc-bach2-naive-rna-plan/candidate-v4/`. Native versions are Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. Independent runtime checking verifies the executable and 929 NumPy plus 1,433 SciPy installed-file hashes.

Machine receipt, source comparison, read guards, logs and synthetic replay are under `work/psc-bach2-naive-rna-independent-review/`. The final receipt is `completion-receipt.json`. Earlier candidates are superseded preparation snapshots, not execution authority. **No real count authorization is issued here.** Root must freeze separately; a successful real execution still supports only this fixed same-cohort RNA question, then this source-cohort branch closes regardless of outcome.

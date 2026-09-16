# PSC BACH2: fixed naïve-CD4 RNA method preparation

**16 September 2026 UTC. Method preparation only. No real matrix or margin payload was opened. No real BACH2 extraction, selected-population RNA quantity, normalization, score, contrast or detection result was computed.** The production preparation/validation also avoids prior quantity JSONs. The independent review's early qualification-receipt access is disclosed below. No network, installation, reference query, clinical join or Git operation was used. The twelve closed qualification deliverables and source receipts were not edited.

The root preparation decision does **not** authorize expression execution. This method implements the [one-candidate design review](psc-bach2-followup-design-review.md). It does not expand the scientific question. A separate root freeze and authorization are still required.

## Fixed scientific target

- **One feature:** `ENSG00000112182` / `BACH2` / `Gene Expression`. Its full-source row is **12,314, one-based**, in all eight feature files.
- **One population:** exact `paper_clusters` labels `C0: CD4+ TN RTE` and `C1: CD4+ TN mature`. Their required short/Seurat crosswalks are C0/0 and C1/1. The fixed annotation union is **9,320 + 7,339 = 16,659 cells** from the 55,460 retained cells. Source-only extras stay excluded.
- **One independent unit:** each of the eight accession-local source donors. Four source labels are `SNP`; four are `noSNP`. Samples are unpaired. Genotype nucleotides, dosage, strand and REF/ALT are not inferred.
- **One contrast:** equal-donor mean `log2(CPM+1)` in `SNP` minus that in `noSNP`. There are no extra genes, separate C0/C1 or CD8 contrasts, ADT proxies, program scores, clinical models or RNA-detection filters.

For donor *d*, pool the selected C0+C1 cells before transformation. Let **B_d** be the released BACH2 units and **R_d** be the released units from **all 36,601 Gene Expression rows in those same cells**. Every ADT row is excluded. Implicit zeros remain part of the full RNA axis.

```text
Y_d = log2(1 + 1,000,000 B_d / R_d)
Delta = mean(Y_d in 4 SNP donors) - mean(Y_d in 4 noSNP donors)
```

No donor is weighted by its cell count, library size or estimated precision. Pooling retains the released C0/C1 mixture and RNA-content weighting. Cell fractions are private context, not covariates or extra endpoints. This method cannot separate population redistribution from within-state RNA change.

## Annotation-only preparation

`scripts/analyze_psc_bach2_naive_rna.py` independently rebuilds the selection from the frozen cell annotation and SDRF. It strips only the exact matching sample prefix from each barcode. Sample identity stays in every key. It verifies every retained-cell join, all eight distinct donor/sample mappings, the four/four labels, the exact state crosswalk, complete ordered feature axes, barcode axes and target identity. Duplicate source symbols are not collapsed or renamed.

An independent cached-reference check confirms the reciprocal `ENSG00000112182` / `BACH2` identity in the **Ensembl 116** human crosswalk. That file has SHA-256 `a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167`. Its exact historical bytes were reconstructed previously from preserved original relationships. The original retrieval was 2026-09-12T17:40:41Z from `https://jun2026.archive.ensembl.org/biomart/martservice`. Provenance is in `data/derived/psc-liver-program-inputs.json`, under the reference reconstruction and historical-response fields. **No new reference was queried. This check does not identify the source's executed build, reconcile variant alleles, or replace the released axis.**

Expected matrix and archive hashes/bytes come from the closed member-integrity and acquisition receipts. The production preparation/validation code does not open, stat or rehash matrix bodies or archives. The production compiler/validator loads no prior donor RNA margins or cohort-qualification quantity JSON. Matrix dimensions come from annotation axes; the future declared nonzero-coordinate count is bound by the frozen matrix hash and will be checked against its streamed count.

The compiled plan includes the exact candidate/review/decision, annotation and receipt pins, all expected source-matrix pins, private selection digest, script/test pins, runtime fingerprint and public schema digest. `.gitignore` must explicitly exclude `work/`; its bytes and conservative no-reinclusion policy are pinned. The historical sample05 acquisition receipt is reused at its original path.

## Numerics and inference

Use sample variances of the four donor scores in each group:

```text
v_g = s_g^2 / 4
SE = sqrt(v_SNP + v_noSNP)
df = (v_SNP + v_noSNP)^2 / (v_SNP^2/3 + v_noSNP^2/3)
CI95 = Delta +/- t(0.975, actual df) * SE
P = 2 * t_survival(abs(Delta/SE), actual df)
```

This is one nominal, two-sided, model-based Welch test. It assumes independent donors and approximately Gaussian donor scores within groups. Four donors per group cannot establish calibration, representative sampling or clinical exchangeability. There is no normality-test switch, permutation claim, cell bootstrap, rank-test rescue or count-likelihood/NB model.

### Stable evaluation, not a tolerance rule

- Source values are parsed as exact positive integers. Plain integer and exactly integral decimal/scientific tokens are accepted within the fixed 100-character/100-digit value cap. Nonfinite, fractional, zero stored coordinates, negative, duplicate, out-of-bounds and out-of-order coordinates fail.
- The parser permits at most 50 million coordinates per matrix. These **predeclared parser bounds**, not inspected RNA values, imply each denominator is below `10^108`. Distinct individual B/R fractions differ by more than `10^-216`.
- Each Y uses a fixed **320-digit Decimal** context. Sample means/variances for the numerical core use exact `Fraction` arithmetic on the represented Y values. Pairwise squared differences avoid cancellation. There is no threshold that reclassifies small variance as zero. An exact source-ratio/score collision check fails closed.
- Production group means, Delta and all eight leaveouts use the algebraically identical product identity with exact `Q_d = 1 + 10^6 B_d/R_d` fractions. For group sizes n1/n0, `Delta = ln(prod(Q_SNP)^n0 / prod(Q_noSNP)^n1)/(n1*n0*ln(2))`. Equal products therefore give exact zero, not a tiny rounding sign. A stable atanh-series evaluation of the logarithm near one retains very small nonzero differences. This is a fixed numerical evaluation rule for the same mean-log estimand, not another contrast.
- Decimal square roots retain SE even when an ordinary floating-point variance would underflow. The actual df remains in [3,6]. Native SciPy evaluates the Student-t quantile and survival probability in double precision. Package versions, the Python executable and installed NumPy/SciPy RECORD hashes **and file contents** are checked.
- JSON numerical outputs are finite Decimal strings. The working precision preserves numerical distinctions; it does **not** imply biological precision or 320-digit accuracy of Student-t functions. A numeric tail underflow leaves P unavailable rather than reporting P=0. It can leave a finite CI available.

Zero B with positive R is retained as Y=0. An exactly constant arm with positive total SE uses ordinary Welch inference, with df=3. If both arms are constant, SE is zero and CI/P are unavailable, even when their descriptive means differ. Nonfinite inference also leaves CI/P unavailable as appropriate. A defined all-eight-donor point remains descriptive. A missing donor, invalid denominator or changed fixed membership/axis/source invalidates the **whole endpoint**. There is no n=7 fallback, donor replacement, imputation or rescue.

## Fixed diagnostics and public/private separation

Private run files contain the eight donor measurements, exact B/R units, selected-cell counts, positive-cell counts/fractions, Y values, C0/C1 metadata fractions and **all eight** donor-indexed leaveout points. Each leaveout recomputes only the 3/4 or 4/3 mean-log contrast. There are no leaveout P values or intervals. Leaveouts are overlapping diagnostics, not eight new replicates or a donor-deletion rule.

The small public projection contains only:

1. Each source group's n=4, mean Y and sample SD.
2. The sole Delta, SE, actual df, nominal CI/P and availability status.
3. Leaveout count=8, complete minimum/maximum, positive/negative/zero counts and strict sign reversals relative to the full point. A zero full point has sign zero and no strict reversal; both leaveout directions remain explicit.
4. Aggregate QC: fixed donor/cell/feature counts, number of zero-B donors retained, and cohort-wide B-positive selected cells. No group-detection endpoint is added.
5. Candidate, plan, selection, source-manifest, implementation, runtime, public-schema and root-authorization digests.

The projection constructs an explicit allowlist. A recursive validator checks every nested key, literal, integer, finite numeric string and exact two-bound CI. It rejects private containers smuggled into an allowed field. There are no donor-indexed expression outputs, joins, individual leaveout points or individual plots in the public schema. Even `public-aggregate.json` first stays inside the ignored run directory for root review; this script never publishes into `data/derived/` or `reports/`.

## Fail-closed CLI and future execution boundary

The default and `--validate-only` first require a simple named preparation directory under the owned work root, excluding run directories. Wrong source/run path roles are rejected before any file body is opened. They then rebuild and compare the source-only plan in memory. They do not open matrices, write result outputs or authorize execution. `--compile-plan` explicitly writes only source selection/pin/schema artifacts into a fresh preparation directory.

```text
.venv/bin/python scripts/analyze_psc_bach2_naive_rna.py
.venv/bin/python scripts/analyze_psc_bach2_naive_rna.py --validate-only
```

A future real run needs **all** of the following:

- Explicit `--execute`.
- A separate root-owned authorization file under `work/psc-bach2-naive-rna-freeze/`, plus its expected SHA-256 supplied by root. Its simple filename must end in `.json` and its size must be at most 8,192 bytes. Wrong-role paths, oversized files and symlinks are rejected **before body access**. This separate freeze namespace is not created by preparation.
- Exact bindings to the code, runtime, complete plan, source manifest, selection, public schema, candidate and named output directory. The preparation decision alone is rejected.
- A new directory under `work/psc-bach2-naive-rna-plan/runs/`. Absolute paths, traversal, symlink components, nested arbitrary names and existing destinations are rejected. Exclusive directory creation consumes the permitted destination. Failed runs retain a failure marker and cannot reuse it.

Only after these gates may the single matrix-open site run. It streams each frozen decoded matrix once, checking bytes/SHA, the MatrixMarket header, exporter comment, dimensions, exact coordinate count, bounds and strict column-major uniqueness. It accumulates **only** selected-cell RNA denominator, selected BACH2 numerator and selected BACH2-positive-cell count. It does not aggregate ADT units, other genes or unselected cells. All eight hashes and fixed gates must pass before the scientific result is returned. There is no repair, alternate source or endpoint selection on failure.

A successful run writes exactly these four files:

```text
authorization.used.json
measurements.private.json
public-aggregate.json
completion.json
```

`completion.json` is written **last**, by an atomic publish from a flushed/fsynced temporary file. It binds the actual bytes and SHA-256 of the three payloads, the exact allowed inventory, the output directory and all plan/selection/source-manifest/code/runtime/schema/authorization/candidate digests. A read-only check then verifies the inventory and payload bindings. Missing, extra, modified or symlinked payloads fail. Any execution failure removes completion markers and preserves `failure.json`; a partial run is not a completed result. CLI success text alone is not a completion certificate. This hash inventory is not a substitute for the separate prospective real-output numerical/source audit that root will prepare.

Future primary and replay runs require distinct root-authorized destinations. Their authorization files, authorization digests, related public provenance and completion hashes can therefore differ. Those are provenance-only differences, not permission to change scientific values or computations. No post-effect repair is allowed.

The authorization is an explicit root-controlled workflow capability, not a claim of cryptographic identity against an operator who can alter the interpreter or bypass this program. No executable authorization file was created in this preparation phase.

## Validation and handoff

- **80 native production synthetic tests pass**, including the fourth-order exact-product cancellation case, no-count default/failure paths, public canaries, same-RECORD runtime tampering, exact boolean authorization, RNA/ADT separation and all-eight gates. Log: `work/psc-bach2-naive-rna-plan/tests-final-v7.log`.
- The final compiled plan is SHA-256 `b6826aa5adb6a3ff26c6f81f5fb16017df5e9580487a96afc2ff491554257b06`. Selection is SHA-256 `63e01409919ae4bd7ab016417a6f285666395464a31a9c53d0129bac381a6257`. The schema is SHA-256 `b1a305315aeba1080b8c15ccf33b6d86641d5071e51d7e5d72bed6a983cc97af`.
- A guarded source-only compile and replay pass. The default native CLI also passes with `execution_authorized=false`, `measurements_executed=false` and zero matrix opens. Recomputed plan/selection bytes match exactly. Logs and opened-path audit receipts are `source-compile-final-v7.*`, `source-replay-final-v7.*` and `default-validate-v7.log` under the preparation work directory.
- **46 independent synthetic adapter tests pass** at the final code pin. An independent 1100-digit direct-log oracle verifies exact product equality and genuine negative contrast near `10^-371`. At the final V7 code pin, paired isolated synthetic-source native runs reproduce all four production artifacts plus the independent receipt byte for byte. The independent completion inventory and numerical comparisons pass. No real matrix is used in those runs. The [independent final review](psc-bach2-naive-rna-independent-review.md) accepts V7/candidate-v4 for a separate root freeze only. It compares all 55,460 retained keys, 16,659 selected rows/six fields, eight donor joins, 36,601 RNA positions, eight matrix pins and 55 annotation/receipt hashes. Its guarded default CLI has zero study-matrix/quantity-JSON/margin open attempts and zero writes. The combined bounded native suite passes **126 tests**.
- The latest guarded native repository suite ran **1145 tests**: seven errors and no failures. It attempted no forbidden BACH2 quantity opens and had no BACH2 failure. The seven historical ETS2/GSE84161 missing-cache/Rscript errors remain. The initial run also found two concurrent liver-verifier fixture issues. Those issues are absent from the final run. This author did not edit unrelated files. `full-tests.log` and `full-suite-receipt.json` contain exact results; the initial full-run log is retained separately.
- The twelve closed deliverables were not written. Their before/after size, inode and modification-time fingerprints match; previously recorded SHA pins are retained. This author deliberately did not rehash quantity-bearing qualification JSONs in this phase. No real run directory or executable root authorization was created.


**Independent receipt-access disclosure:** An early independent audit loaded `data/derived/psc-bach2-cohort-qualification.json` for source/member hashes, dimensions, barcode counts and strict-order fields. That receipt also contains prior donor RNA/ADT technical totals. The reviewer reports no extraction, comparison, selection or analysis of those totals and no matrix/margin-payload or BACH2 measurement access. Subsequent final comparisons use individual integrity/acquisition receipts instead. The reviewer also reports root-authorized opaque hashing of the eight compressed archives, without decoding members; the production preparation/validation code does neither. This historical receipt parse is recorded explicitly; it is not described as a new RNA analysis or as a file that was never read.

The stable source-only plan is `work/psc-bach2-naive-rna-plan/candidate-v4/`. Its selection and schemas are private. Exact final file and runtime pins, source-only audit logs and preservation checks are listed in `work/psc-bach2-naive-rna-plan/completion-receipt.json`. Intermediate drafts and failed test logs remain available; no source value was corrected.

Preparation defects were found before freeze. An early public projection copied nested CI/scalar values without recursive type checks; independent canary tests exposed this and now pass after strict schema enforcement. The first df assertion compared a 320-digit result with a default 28-digit expected Decimal; it was corrected to a rational-error check. Runtime checking was strengthened from RECORD digest alone to verifying installed bytes. Exact-product point evaluation removes cross-arm cancellation signs without changing the estimand. A final pre-authentication file-role check rejects a count or annotation path passed as `--authorization` before opening any body. The regression uses real-source path names under an open-denying mock, not real source bodies. All repairs used synthetic data or source-only annotations.

## Interpretation and stop rule

This is a precision/influence check of the paper's same-cohort broad naïve-CD4 RNA observation, not exact historical Figure S2B normalization/test replication or independent biological validation. Released source-count-format units are **not certified untouched original UMIs**. The count-likelihood gate remains unpassed. The Cellranger 6.1.1 header versus Methods/IDF 3.0.2 and reference-release uncertainty remain unchanged.

The endpoint cannot identify transcription rate, absolute RNA, mature miR4464, BACH2 protein, translation, causal allele effects, mechanism, clinical outcome or within-state regulation. The panel has no explicit BACH2 ADT feature. Treatment, ancestry, batch, clinical selection, LD and pooled-state explanations remain. A null or narrow interval does not establish equivalence without an independently justified margin.

If root later authorizes the fixed extraction: one primary contrast, these diagnostics, then close regardless of sign, P, width, sparsity or influence. A failed gate yields unavailable/PARK, not another gene/state, donor deletion, normalization change or new acquisition. **Real expression execution remains blocked until that separate root action.**

# BACH2 naive-CD4 aggregate reporting design

## Scope and current status

This is a presentation-only reporter for the closed, single BACH2 endpoint. It does not call the scientific implementation, real-source plan validation, native execution, or the independent matrix auditor. Root owns the biological narrative and any later publication. The prepared artifacts remain in ignored `work/psc-bach2-naive-rna-reporting/` until root gives a separate hash-bound publication approval.

Authority is `work/psc-bach2-naive-rna-freeze/root-reporting-handoff.json`, SHA-256 `53f636477b1f8c33c44308b3b48ab3f3a9b448f53422dc17f1e07fbdd935fab6`. The reporter reads exactly that file and its six pinned inputs: the outer freeze, primary native completion, approved public aggregate, public schema, aggregate audit receipt, and mechanical replay receipt. It follows no paths inside those metadata files. The original source URLs, acquisition history and scientific method remain in the [frozen protocol](../docs/psc-bach2-naive-rna-protocol.md) and its qualification references.

The approved public aggregate has 5,391 bytes and SHA-256 `3f9383678afc4a576c12f89225e3148ca07ff65486c8187e8579f4dfd4ac9837`. The outer freeze is `df5cc30e11fd79424d8b8c429f4c62e48c46cda657afb4601d35d7e316fc346e`. The native primary completion is `850ff604b819fb0aae660fabc9a42248d511a5a334b07764afdb46d3dd961c4f`. These are existing scientific artifacts, not products of the reporter.

No private measurements, matrices, archives, technical margins, donor joins or individual omissions are opened. No network access, dependency installation, Git operation, source change or new statistical computation is part of reporting. The default and `validate` commands write nothing.

## Small, complete presentation

Each fresh preview mirrors these intended repository paths. It contains ten payload files plus a manifest and a last-written completion record:

| Intended path | Content |
| --- | --- |
| `data/derived/psc-bach2-naive-rna/public-aggregate.json` | Exact approved source bytes; all decimal strings remain unchanged. |
| `data/derived/psc-bach2-naive-rna/public-schema.json` | Exact pinned public schema bytes. |
| `data/derived/psc-bach2-naive-rna/endpoint.csv` | Exactly one row: schema version, status, endpoint, quantity, Delta, SE, actual df, both 95% CI bounds, nominal two-sided P, inference status and all eight provenance hashes. |
| `data/derived/psc-bach2-naive-rna/groups.csv` | Exactly two rows in `SNP`, `noSNP` order. Each retains n, mean, **sample SD**, the common endpoint/status fields and all eight provenance hashes. |
| `data/derived/psc-bach2-naive-rna/influence-qc.json` | Exact common fields, every fixed aggregate omission and QC field, and complete provenance. No donor-indexed omission values. |
| `data/derived/psc-bach2-naive-rna/aggregate-verification.json` | Exact pinned independent-audit receipt, not a new audit. |
| `data/derived/psc-bach2-naive-rna/mechanical-replay-verification.json` | Exact pinned root mechanical-replay receipt, not a new scientific replay. Private-artifact names/hashes in this receipt are metadata, not private values. |
| `data/derived/psc-bach2-naive-rna/reporting-info.json` | Reporting/read scope, precision rule, renderer settings and disclosure limits. |
| `reports/figures/psc-bach2-naive-rna.png` and `.svg` | One contrast point and its existing nominal Welch 95% CI. No individual donor/cell points or distributions. |
| `data/derived/psc-bach2-naive-rna/manifest.json` | Exact payload inventory, byte counts/SHA-256, input and reporter pins, CSV column lists/row grid and JSON projection schemas. |
| `data/derived/psc-bach2-naive-rna/completion.json` | Complete expected inventory and manifest pin. Written atomically last. Its origin preview path is provenance and is preserved on later publication. |

CSV decimal fields are copied as strings without arithmetic or rounding. CSV uses UTF-8, LF and an explicit fixed column order. Generated JSON uses UTF-8, sorted keys, two-space indentation, a terminal LF and no nonfinite JSON constants. There is no gzip archive.

## Closed validation and integrity gates

The new validator recursively interprets only the pinned public schema's descriptor vocabulary. Every object has an exact key set. It rejects unknown nested keys, private canaries, duplicate JSON keys, wrong types, booleans used as integers, nonnumeric/nonfinite decimal strings and malformed CI lists. It checks the two-arm n=4 grid; fixed RNA/selection/donor denominators; CI order around the point; positive SE; actual-df and probability bounds; all eight omissions, sign/extrema/count consistency; aggregate QC bounds; and all eight frozen provenance bindings.

This is deliberately hash-specific. Only the actual `available_nominal_Welch` status is accepted. Any changed/unseen status fails closed rather than gaining fabricated intervals, P=1, zero values or an arbitrary `passed` subtree. The exact source digest also rejects otherwise numerically valid alterations. No mean, variance, effect, inference, omission or detection contrast is recomputed. Decimal comparisons are validation, not a new endpoint calculation.

A preview must have a new simple name directly under the ignored reporting namespace. The reporter requires the existing `work/` ignore rule. It rejects existing destinations, absolute/aliased/traversal paths, symlinks, hardlink aliases, special files and source-path collisions. It validates before creating the preview. Payload writes are exclusive. Failure removes the incomplete preview. A temporary completion file is atomically renamed only after all payloads and the manifest are written. An unauthenticated or incomplete directory is not a completed bundle.

`verify` checks the exact recursive inventory, rejects extra files/directories and symlinks, authenticates the supplied completion digest, and regenerates the presentation in memory for exact byte comparison. It does not execute the scientific method. Fresh preview replay preserves all payloads and the manifest byte for byte; only the completion's origin preview directory and consequent completion digest differ.

## Figure and precision

PNG and SVG come from the same Matplotlib figure. The renderer fixes Matplotlib `3.11.2`, normal DejaVuSans.ttf SHA-256 `3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22`, default rc settings plus explicit overrides, SVG hashsalt `psc-bach2-naive-rna-public-report-v1`, path-based SVG glyphs, 160 dpi, fixed axes/size and fixed creator/software metadata. SVG date metadata is absent; no date is added to PNG. A version/font or fixed display-range change fails closed.

Only plotting coordinates use float conversion. Labels round the existing contrast/CI/P to three decimal places. This is display formatting, not a new numerical result. The source's hundreds of computational digits are not biological precision. The plotted interval is the sole contrast's **nominal 95% CI**, not either group's sample SD and not an omission range. The PNG is visually inspected before handoff; root must also inspect it before approval.

## Interpretation boundary

The reporter does not supply a biological conclusion. The [protocol](../docs/psc-bach2-naive-rna-protocol.md) fixes one gene, the C0+C1 union, eight male PSC donor units, four per source group, equal donor weighting and the all-36,601-RNA-row denominator. Donors, not cells, are the replicates. Count-format released RNA units are not certified untouched UMIs. Pooling retains observed C0/C1 mixture and RNA-content weighting. Clinical imbalance, LD, batch and unverified participant-to-clinical linkage remain limits.

There is no equivalence, causal-nucleotide, protein/translation, treatment or independent-biological-replication claim. No omission CI/P, group detection contrast, fold-change conversion, power, alternative gene/state or model is added. Public n=4 aggregate release is not formal differential privacy or a guarantee against inferential disclosure. Close this branch after reporting.

## Commands and later root approval

Use the native project environment:

```text
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py validate
.venv/bin/python -B -m unittest discover -s tests -p test_report_psc_bach2_naive_rna.py -v
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py preview --output work/psc-bach2-naive-rna-reporting/preview-v1
.venv/bin/python -B scripts/report_psc_bach2_naive_rna.py preview --output work/psc-bach2-naive-rna-reporting/replay-v1
```

The ignored stable handoff records the actual completion pins, file inventory, targeted test receipt, exact verify commands, another fresh replay command, and exact root-only approval/publication commands. Do not reuse an existing preview directory. A clone alone is incomplete: the seven original pinned public/metadata inputs in their exact ignored paths must be restored for validation or reproduction.

After root reviews the complete preview and figure, root may use `approval-template` with the exact bundle/completion pin, then deliberately save that template as a new JSON file in `work/psc-bach2-naive-rna-freeze/`. Root approval must bind the reporter digest, public aggregate digest, preview origin, completion digest and all three fixed publication destinations. `publish` requires both that separate approval file and its SHA-256, re-verifies the exact presentation, rejects existing/aliased destinations and writes completion last. It copies only verified in-memory public bytes. The child reporter author does not create an approval file or run publication. Existing frozen files and source records remain unchanged.

## Targeted verification

The native test suite covers exact source and provenance pins; recursive privacy and finite-value failures; fixed groups/omission denominators and statuses; complete CSV/JSON projection without loss of decimal strings; bounded allowed input reads; exact full presentation/figure replay; directory, alias, symlink, hardlink and inventory rejection; wrong completion hashes; write-failure cleanup; and publication refusal without exact separate approval. Publication itself is not executed by these tests. Synthetic temporary test files remain inside the ignored reporting workspace and are removed after each test.

Current reporting checks: all 27 targeted native tests pass. Fresh `preview-v1` and `replay-v1` each have 12 files. All ten payloads and the manifest replay byte for byte; only the completion origin path and its digest differ. Both authenticated `verify` commands pass. The PNG was visually inspected: the single point, interval, zero reference, endpoint labels and rounded values are readable and unclipped. This is author inspection, not root publication approval. The original six input pins and frozen protocol digest remain intact.

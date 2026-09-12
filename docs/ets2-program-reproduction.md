# Reproduce the ETS2 program benchmark

This phase uses public summary statistics and the previously downloaded liver-atlas count matrices. It requires no AlphaGenome key, GPU, new sequencing or controlled-access application. The [methods](ets2-program-methods.md) describe the question, scoring rules and limitations; the [frozen plan](../config/ets2-program-plan.json) pins the original inputs and dependencies.

## Environment and source cache

Use Python 3.12 in an isolated environment. The recorded package versions are in [requirements-ets2-program.txt](../requirements-ets2-program.txt), which includes the preceding liver-atlas requirements and adds the workbook reader.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-ets2-program.txt
python scripts/fetch_liver_atlas_sources.py
python scripts/fetch_liver_atlas_sources.py --manifest config/ets2-benchmark-sources.json
python scripts/fetch_liver_atlas_sources.py --manifest config/ets2-program-sources.json
```

Run commands from the repository root. Existing valid caches are reused. Add `--offline` to each fetch command for verification without network access. The liver objects total approximately 2.66 GB; the original ETS2 public release is approximately 35 MB, and the two additional publisher workbooks total approximately 13.18 MB. Raw files stay under ignored `data/raw/`.

Exact bytes and hashes are required. A changed live metadata response may prevent a fresh download from reproducing the historical snapshot. Do not replace a pinned file with a newer response or rewrite the plan to make a checksum pass. The alternate publisher-host retrieval and its preceding DNS failures are recorded in the new source manifest; no failed HTML page was treated as a workbook.

## Rebuild the analysis

```bash
python scripts/audit_ets2_sources.py
python scripts/benchmark_ets2_experiments.py
python scripts/reconstruct_ets2_inhibitor.py
python scripts/analyze_ets2_program.py
python scripts/verify_ets2_program.py
python scripts/plot_ets2_program.py
python -m unittest discover -s tests -p 'test_ets2*.py' -v
```

The first step extracts and audits source gene sets, final differential-expression summaries, saved rank vectors and experiment-design metadata. It preserves disagreements instead of choosing the sign that best fits the expected biology. The experimental benchmark uses final publisher tables for the primary comparisons. These are analyses of published summary statistics; the original donor-level RNA fits are not recreated.

The separately requested inhibitor reconstruction audits the intended R workflow and calculates the nine deterministic weighted enrichment scores from the saved 500-nM rank vector. It compares those scores with the archived results and compares existing archived NES values with publisher Figure 5. It does not re-fit the missing counts, reproduce the stochastic fgsea null or revise the frozen treatment-direction exclusions. Original full differential rows and rank values are stored once in the source-audit gzip files; the experimental benchmark does not duplicate them.

The atlas step checks inputs, reads raw counts in chunks and aggregates them by donor, assay and original cell-population label. It applies the frozen membership and mapping rules, constructs the expression-matched references, and reports all nine programs at both cell-count thresholds. The gene universe, control exclusions, random seed, sampling order and comparison list are fixed. Matched reference sets are reused across donors and thresholds.

The figure reads aggregate outputs. Public files contain program membership/mapping, coverage and matching audits, context summaries, planned paired contrasts, descriptive correlations and experimental summaries. Donor-level counts, individual score rows and random-set selections remain in ignored `work/` and are not versioned. Their recorded hashes support local verification without publishing the underlying records.

The full repository test command is `python -m unittest discover -s tests -v`; earlier phases need their separately documented dependencies. The recorded research environment already contains them.

## Isolated replay

The numerical scripts support alternate output and audit paths so a replay can preserve the recorded result:

```bash
python scripts/analyze_ets2_program.py \
  --output-directory work/ets2-program-replay \
  --audit-path work/ets2-program-replay-audit.json \
  --chunk-size 8192
python scripts/benchmark_ets2_experiments.py \
  --output-directory work/ets2-experiment-replay \
  --audit-path work/ets2-experiment-replay-audit.json
```

Compare numerical output hashes by basename against the corresponding recorded audit. File paths, execution timestamps and processing chunk counts can differ without changing numerical results; those differences should be explained rather than removed from provenance. Verification of saved hashes is distinct from rerunning the calculations.

The recorded production run used chunks of 4,096 observations; its replay used 8,192. All six public atlas tables and all four experimental benchmark tables reproduced exactly. The [integration verification](../reports/ets2-program-verification.json) records those checks, independent count/statistic comparisons and the 120-test full-suite result. This is a record of the completed run; it does not substitute for executing the commands on a new environment.

Reproduction confirms that the calculations follow their recorded inputs and rules. It does not establish independent biological validation, prove ETS2 dependence in a liver cell, or turn a matched-gene percentile into a clinical or population probability.

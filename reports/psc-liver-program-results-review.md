# Independent results review: fixed pediatric liver programs

## Conclusion and finite close

**No group mean difference was detected in the prespecified Core panel.** All 30
Core tests have nominal P ≥ **0.382808**. None meets nominal P < 0.05, BH q < 0.05
or BY q < 0.05. Every Core BH q is **0.9593488442261611**; every BY q is **1**.
The lack of discoveries is not an adjustment-only finding. These results do not
establish equivalence, absence of relevant biology or absence of inflammation.

The verified bulk-tissue analysis provides no positive statistical support for a
PSC-associated increase in the focal ETS2 gRNA1-down composite relative to AIH.
It does not test or refute ETS2 dependence in a particular cell type. It is not
independent mechanistic confirmation of MacroMap.

**Close this fixed branch after reporting, regardless of its null result.** Do not
switch interfaces, signatures, normalizations, thresholds or subgroups to obtain
a signal. No rescue analysis is proposed. The epithelial IL17A-response transfer
question remains **HOLD and unfulfilled**; the unsigned ECM set did not replace it.

## What this review checked

This review read only the authenticated aggregate, frozen plan/root seal, root
reporting handoff and aggregate verification/replay receipts. It did not open raw
counts, private units, library totals, score vectors, bootstrap arrays or individual
omissions. It did not rerun effects, tests, multiplicity corrections or resampling.
Numbers below are direct projections of the accepted aggregate.

- All **10 programs × 3 contrasts × 2 interfaces = 60** unique planned rows remain.
  All **30 Core hypotheses** remain in the sole discovery family. All 30 Strict
  rows are diagnostic sensitivity results, with null discovery q fields.
- All 20 prespecified score-eligibility gates pass. All 60 effects and Welch tests
  are available, with finite positive arm variances. No Core P is unavailable or
  substituted. `discovery_eligible=true` means a Core test has a usable P value;
  it does **not** mean a discovery was made.
- Every Welch interval includes zero. All 60 bootstrap intervals also include
  zero, with status `available` and all **10,000** requested draws valid. These
  intervals are **pointwise**, not multiplicity-adjusted or simultaneous. The
  bootstrap supplies no additional P-value family.
- Every public omission summary uses **47** units for PSC−AIH, **34** for PSC−ASC
  and **47** for ASC−AIH. All relevant effects are defined. The 64-unit private
  inventory is not a public stability denominator.
- All **30 matched Core/Strict effect signs agree**. Both interfaces reuse the same
  patient-samples and largely overlapping genes. Agreement is annotation
  sensitivity, not independent replication or proof of zero effect.

The authenticated frozen-adapter receipt records a complete source/native audit:
2,630,144 count/cache cells, 64 totals, 1,280 scores, 60 rows, the 30-test correction
vector, 640,000 RNG indexes, 600,000 bootstrap differences and 3,840 omissions.
The authenticated mechanical replay receipt records **all 12 native artifacts
byte-identical**. This review relied on those receipts rather than repeating
private computations. Deterministic replay is not another biological experiment.
No native logCPM matrix comparison is claimed.

## Focal result and what its uncertainty means

The primary contrast is **PSC minus AIH**, with 17 and 30 source patient-samples.
Core exact symbol–Ensembl bijection is the primary interface.

| Focal ETS2 gRNA1-down result | Core | Strict, diagnostic only |
|---|---:|---:|
| Mean-score difference | −0.021094 | −0.021980 |
| Pointwise 95% Welch CI | [−0.154686, +0.112497] | [−0.158645, +0.114686] |
| Two-sided raw P | 0.751928 | 0.747493 |
| Pointwise 95% bootstrap CI | [−0.150820, +0.106392] | [−0.154530, +0.108432] |
| Relevant omissions with strict sign reversal | 1/47 | 1/47 |

Units are differences of positive equal-member mean **log2(CPM+1)** scores, not
absolute RNA abundance, a program-level RNA fold change, TF activity or a clinical
scale. The Core Welch interval is wide relative to the −0.021094 estimate and
contains both negative and positive differences. Under the working sampling/model
assumptions, its numerical bounds constrain this defined unadjusted bulk-score
contrast; values beyond them are less compatible with this pointwise procedure,
not impossible. The bounds cannot be transferred to cell-specific ETS2 activity,
fibrosis severity, treatment response or other clinical outcomes.

There was **no prespecified equivalence margin** or clinical calibration. A high
P value is not a probability that the null is true. BH q ≈ 0.959 is not 95.9%
certainty of no biology. This does not establish a precisely zero effect, clinical
equivalence or exclusion of every biologically relevant contrast.

Core single-unit omission estimates range from **−0.035960 to +0.003578**.
One omission crosses to the opposite sign; maximum absolute change is **0.024672**.
The estimate's direction is therefore not invariant to every single-unit omission.
This does not identify a bad sample or justify its removal. The omission range is
**not a confidence interval**, and 47 correlated omissions are not 47 independent
replications. No omitted identity was inspected.

## Complete fixed panel, without selecting favorable rows

Each cell is **reported effect; two-sided raw P**, rounded only for display. The
unrounded 60-row aggregate remains the reference. All 30 Core BH/BY values are
reported above; Strict has no discovery q values. All pointwise Welch and bootstrap
intervals include zero, as checked across the complete panel.

### Core: sole 30-test discovery family

| Program | PSC − AIH | PSC − ASC | ASC − AIH |
|---|---:|---:|---:|
| `ets2_g1_dn` | -0.021094; 0.751928 | +0.018065; 0.803419 | -0.039159; 0.631867 |
| `ets2_g1_up` | -0.016245; 0.692576 | -0.011236; 0.829505 | -0.005009; 0.924823 |
| `ets2_g2_dn` | -0.021150; 0.753835 | +0.024189; 0.736975 | -0.045339; 0.576813 |
| `chr21_dn` | -0.054015; 0.587119 | +0.046859; 0.665321 | -0.100874; 0.382808 |
| `inflammation` | -0.027186; 0.590666 | +0.008642; 0.873119 | -0.035828; 0.539689 |
| `interferon_gamma` | -0.050813; 0.449158 | +0.004815; 0.942167 | -0.055628; 0.435494 |
| `oxidative_stress` | -0.013897; 0.723955 | -0.017500; 0.713870 | +0.003603; 0.944505 |
| `apoptotic_signaling` | -0.019991; 0.663167 | -0.017075; 0.754230 | -0.002915; 0.959349 |
| `ets2_g1_dn_without_comparators` | -0.009448; 0.879057 | +0.014823; 0.827232 | -0.024271; 0.751942 |
| `reactome_ecm_organization` | +0.029016; 0.680822 | -0.026978; 0.758929 | +0.055995; 0.513375 |

### Strict: identity sensitivity, not a replacement discovery analysis

| Program | PSC − AIH | PSC − ASC | ASC − AIH |
|---|---:|---:|---:|
| `ets2_g1_dn` | -0.021980; 0.747493 | +0.019144; 0.795230 | -0.041124; 0.620559 |
| `ets2_g1_up` | -0.012326; 0.792054 | -0.007696; 0.893416 | -0.004630; 0.936374 |
| `ets2_g2_dn` | -0.021778; 0.752247 | +0.025473; 0.728205 | -0.047251; 0.567235 |
| `chr21_dn` | -0.055011; 0.589962 | +0.051632; 0.641804 | -0.106644; 0.368305 |
| `inflammation` | -0.027548; 0.587582 | +0.008428; 0.876858 | -0.035976; 0.539482 |
| `interferon_gamma` | -0.050597; 0.445135 | +0.005154; 0.937070 | -0.055751; 0.427883 |
| `oxidative_stress` | -0.014325; 0.716799 | -0.018265; 0.703668 | +0.003940; 0.939616 |
| `apoptotic_signaling` | -0.019924; 0.664011 | -0.017202; 0.752029 | -0.002722; 0.962043 |
| `ets2_g1_dn_without_comparators` | -0.010079; 0.874287 | +0.016097; 0.815823 | -0.026176; 0.737774 |
| `reactome_ecm_organization` | +0.028650; 0.685902 | -0.027218; 0.758299 | +0.055869; 0.516672 |

The minimum Core raw P, **0.382808**, belongs to `chr21_dn`, ASC−AIH. This is an
inventory summary, not a selected lead or rescued finding. Differences between
program P values do not test differences between program effects. Even a
significant program beside a nonsignificant comparator would not by itself
establish comparative or ETS2-specific activity. Here no program is nominally
significant. Correlated programs and contrasts are not independent replications.

## Coverage and source labels

| Program | Original denominator | Core mapped | Strict mapped |
|---|---:|---:|---:|
| `ets2_g1_dn` | 927 | 848 | 825 |
| `ets2_g1_up` | 668 | 612 | 536 |
| `ets2_g2_dn` | 875 | 809 | 789 |
| `chr21_dn` | 141 | 126 | 122 |
| `inflammation` | 815 | 750 | 744 |
| `interferon_gamma` | 139 | 121 | 118 |
| `oxidative_stress` | 436 | 403 | 400 |
| `apoptotic_signaling` | 588 | 544 | 542 |
| `ets2_g1_dn_without_comparators` | 927 | 695 | 673 |
| `reactome_ecm_organization` | 321 | 301 | 299 |

The nonoverlap row uses its **original parent coverage gate plus at least ten
retained genes**, as frozen before effects. Its component `coverage_gate=false`
for 695/927 and 673/927 is expected; its `parent_gate=true` and
`score_status=available` are not contradictory. Do not mislabel it as a failed
endpoint or claim every component-level 80% flag is true. All approved score-level
gates pass. This is annotation eligibility, not proof that every member is
expressed or that a signature is biologically specific.

ECM remains an **unsigned organization pathway**, not a directional fibrosis-
activity score. Its denominator is 321 official symbols, with 301 Core and 299
Strict mapped members. Literal ETS2 is absent from the original set; no new target
exclusion was applied. The nonoverlap ETS2 set still removes only the original
four comparator sets, not ECM. Source companion directions were retained, not
inverted after viewing these results.

## Interpretation limits that still apply

1. **No healthy reference group.** PSC, ASC and AIH are disease groups. Similar
   scores between them do not show normal expression, lack of inflammation or
   absence of a program in all groups. The comparator lists do not span all
   possible generic inflammatory biology.
2. **No clinical/stage/batch adjustment.** These unadjusted observational
   comparisons do not identify a PSC-specific causal effect, disease progression
   or a treatment effect.
3. **Independence is assumed, not independently verified.** The paper asserts
   patients, but an independent person/visit crosswalk is unavailable. Welch and
   patient-sample bootstrap validity remain conditional on that assumption and
   their other working assumptions. Numerical verification does not repair it.
4. **Whole-biopsy and gene averaging can conceal heterogeneity.** Cell mixtures,
   opposing cell-specific changes or opposing member-gene changes could yield a
   small bulk composite difference. This analysis did not estimate those patterns;
   it does **not** show that masking or dilution occurred. That possibility is
   not a rescue explanation demonstrated by these data.
5. **Relative, annotation-dependent composites.** Full 41,096-row library
   normalization, the pseudocount, fixed memberships and mapping choices define
   the estimand. It is not all transcription, cell-intrinsic pathway activity,
   protein abundance or TF dependence. Current identifier reconciliation does
   not reconstruct historical read assignment or counting-GTF choices.
6. **Computation is not mechanism or clinical evidence.** The audit and exact
   replay establish agreement within the frozen protocol, not biological
   replication, causal identification or a clinical recommendation.

## Freeze provenance and disposition

The actual external root freeze is present and authenticated. The plan is
`state=frozen`, has explicit root authorization and points to the authenticated
external seal. Native `root_freeze_written=false` only states that the **native
writer** did not create that external freeze; it does not mean no freeze exists.
The handoff records checkpoint `c1ce7f25934c3797eed4660cfed414749c8b46bf`.
This review did not run Git or independently inspect checkpoint history.

| Authenticated artifact | SHA-256 |
|---|---|
| Aggregate public summary | `350ee5d1e447b1cab4621595611de99f7a787af40102187ca0377257c9a598d0` |
| Frozen plan | `5e110cfabca0951846d5101063cada8761f021122796649d195174dc6f75315d` |
| External root seal | `59eb70cb6a539d84799300332c1ece54bf5384e9655d26884dacf06b62df7c53` |
| Full source/native verification receipt | `9a5827803aa44a9124e59bb491d9b9f6d021b1c5fdd79d82cfbed73c1b856004` |
| Complete mechanical replay receipt | `00c4e01b10c01f0345933b08bf0a178ba6d94aee8581ac4242f8406c6e6b79e8` |

**No genuine execution, inventory or interpretation inconsistency was found** in
this permitted review scope. The nonoverlap coverage flag and native freeze flag
need explanations, not repairs. The finite branch is closed without rescue. The
held epithelial IL17A transfer question remains unanswered. GSE159676 remains
HOLD; the NoPSC atlas remains PARK.

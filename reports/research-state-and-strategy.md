# PSC research: current state and next decision

Assessment started September 15, 2026, at the user's request to pursue useful research with public data and compute. This report records the initial assessment and the completed first analysis separately. A completed computational analysis is not itself a biological breakthrough.

## Completed first analysis and next decision

The [PSC blood reanalysis](psc-blood-replication.md) is complete: 93 of 11,074 shared genes meet the prespecified same-direction BH conjunction, but none passes BY or either stronger disease-control conjunction. No established PSC-specific signature or new causal target follows. All 236 independent numerical checks pass, and the 12 checked primary artifacts replay byte-for-byte.

The next proposed compute priority is a donor-level PSC-by-IL-17A response interaction in GSE239283 cholangiocyte organoids. It has not been executed. MacroMap is a qualified reserve context dataset; its program effects also remain unmeasured. See the [updated roadmap](../ROADMAP.md).

## What the existing research establishes

- The variant catalogue, allele/build reconciliation, public-source recovery and numerical verification are unusually thorough. The [evidence log](../docs/evidence-log.md) preserves both favorable and unfavorable results.
- AlphaGenome has **not** demonstrated consistent prioritization performance in the small matched comparison. More model predictions alone are not the next priority.
- PRKD2 has aligned monocyte RNA associations and plausible regulatory context. Missing executed-study model/covariance or complete-component inputs still prevent a qualified shared-cause conclusion. More reference-LD fits cannot recover missing study facts.
- UBASH3A has a concrete splice hypothesis and conditional sequence consequences. The later direct RNA/protein searches did not establish the proposed complete RNA or protein. Coverage-limited negative searches do not establish absence.
- The published ETS2 program is measurable in liver myeloid cells, but transfer varies by cell state and overlaps broad responses. The program is not yet a validated ETS2-specific activity score or PSC treatment-response measure.
- The existing PSC liver atlas supports cell context. Chemistry and processing confounding and small eligible control groups block the intended immune-cell case-control comparison.

## Selected path: PSC blood cross-cohort replication

The next result should add **biological independence and a discriminating comparison**, not merely more source files or verification counts. The selected path is the Norwegian GSE177044 PSC case-control stratum, followed by Polish GSE119600 source-labeled adult arrays. The [prospective protocol](../docs/psc-blood-replication-protocol.md) fixes gene-level replication and disease-control comparisons before new effects. It does not claim untouched validation: the original RNA paper already tested a classifier in the Polish cohort.

The completed parallel reviews considered:

1. The prepared GSE255234 four-donor macrophage context study. This is an execution pilot, not a strong PSC-specific discovery endpoint.
2. MacroMap, reported to contain many iPSC donor lines and conditions. Its accessible matrix, exact measurement and donor design still need qualification. It could test program specificity across stimuli at a useful donor scale, but it is not a PSC cohort.
3. Independent public PSC liver datasets with disease controls and a feasible replication design. Observational tissue signals must not be labeled causal or cell-intrinsic without further evidence.

Selection preceded target-effect inspection. All primary sources have been downloaded and structurally checked. Norway supplies 220 PSC cases and 77 controls with complete age, sex and within-plate metadata; Poland supplies 275 source-labeled adult arrays including healthy, UC, PBC and CD comparison groups. Exact age, sex and treatment adjustment is unavailable in the Polish sample records. Existing frozen protocols and results are not rewritten to make a new dataset pass.

## Local execution status

This checkout initially lacked the historical ignored `data/raw/`, `work/` and `.venv` caches. A new project Python 3.12.14 environment now contains the pinned ETS2-program dependencies and supporting analysis packages. The initial full test attempt ran 255 tests and reported seven errors: four class setup errors need missing historical caches, and three tests need the absent historical GSE84161 R runtime. No new biological code had been changed before this run. These failures are environmental/reproduction limits, not a passing verification claim.

No new effects had been calculated at that initial baseline stage. The subsequent blood analysis and independent verification are now complete. The final full suite reports 322 tests with the same seven historical environment errors; the 67 new analysis tests pass. No researchers have been contacted, no services enrolled in, and no paid compute provisioned.

"""Synthetic-contract tests; no real expression data or target effects are read."""
from __future__ import annotations

import copy
import gzip
import importlib.metadata
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.analyze_psc_blood_rnaseq import (
    AnalysisError, FILTER_SETTINGS, MODEL_SETTINGS, PRIMARY_COUNTS, PRIMARY_PLATES,
    SENSITIVITY_SETTINGS, adjust_pvalues, align_primary_samples, build_design,
    check_design, check_pin, expression_filter, fit_negative_binomial, fit_ols_hc3,
    fit_subgroup_sensitivity, leave_one_plate_out, log2_cpm, parse_metadata,
    ratio_size_factors, read_counts, read_metadata, run_analysis, sha256, validate_counts,
    validate_plan, versionless_ensembl,
)


def metadata_fixture(n: int = 96) -> pd.DataFrame:
    i = np.arange(n)
    source = np.asarray(["Control", "PSC", "PSCUC"])[i % 3]
    sex = np.where(i % 2, "male", "female")
    frame = pd.DataFrame({
        "sample_id": [f"sample{i:03d}" for i in range(n)],
        "sample_accession": [f"GSM{i:03d}" for i in range(n)],
        "source_disease": source, "disease": np.where(source == "Control", "Control", "PSC"),
        "age_source": (20 + (i * 7) % 57).astype(str), "age": (20 + (i * 7) % 57).astype(float),
        "sex_source": sex, "sex": sex, "plate": np.asarray(["P1", "P2", "P3"])[i // (n // 3)],
    }).set_index("sample_id", drop=False)
    return frame


def json_metadata(frame: pd.DataFrame) -> dict:
    return {"GSE177044": [{"sample_accession": row.sample_accession,
        "title": [row.sample_id], "characteristics_ch1": [
            f"disease: {row.source_disease}", f"age: {row.age_source}",
            f"Sex: {row.sex_source}", f"sequencing plate: {row.plate}"]}
        for row in frame.itertuples(index=False)]}


def plan_fixture() -> dict:
    source = {"path": "not-read", "sha256": "a" * 64, "url": "https://example.org/source",
              "retrieved_at_utc": "2026-01-01T00:00:00Z"}
    return {"status": "frozen", "frozen_at_utc": "2026-01-01T00:00:00Z", "rnaseq": {
        "counts": [dict(source)], "metadata": dict(source), "expected_sample_count": 297,
        "expected_feature_count": 63677, "annotation_columns": ["Geneid", "gene_name"],
        "metadata_series": "GSE177044", "selection": {"plates": list(PRIMARY_PLATES),
            "source_disease_counts": dict(PRIMARY_COUNTS), "metadata_samples": 1035},
        "filter": dict(FILTER_SETTINGS), "model": copy.deepcopy(MODEL_SETTINGS),
        "sensitivity": dict(SENSITIVITY_SETTINGS), "n_cpus": 1,
        "software": {"python": ".".join(map(str, sys.version_info[:3])), **{
            name: importlib.metadata.version(name) for name in ["numpy", "pandas", "scipy", "pydeseq2"]}},
        "known_target_source_labels": ["IFIT1", "G0S2", "UBASH3A", "ETS2", "PRKD2", "PFKFB3", "BCL2L11"]}}


class TestCountsAndIds(unittest.TestCase):
    def write(self, directory: str, text: str, compressed: bool = False) -> Path:
        path = Path(directory) / ("counts.csv.gz" if compressed else "counts.csv")
        if compressed:
            with gzip.open(path, "wt") as stream:
                stream.write(text)
        else:
            path.write_text(text)
        return path

    def test_original_labels_integer_values_and_full_zero_feature_retained(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write(temp, "Geneid,gene_name,B,A\nENSG00000000001.2,SOURCE,5,7\nother,,0,0\n", True)
            data = read_counts(path)
        self.assertEqual(list(data.counts.index), ["B", "A"])
        self.assertEqual(list(data.counts.columns), ["ENSG00000000001.2", "other"])
        self.assertEqual(data.counts.to_numpy().dtype, np.dtype("int64"))
        self.assertEqual(data.features.iloc[0].ensembl_gene_id_versionless, "ENSG00000000001")
        self.assertEqual(data.features.iloc[1].source_gene_name, "")
        self.assertEqual(data.features.iloc[1].ensembl_gene_id_versionless, "")
        filt = expression_filter(data.counts, min_count=5, min_samples=2)
        self.assertEqual(filt.eligible.tolist(), [True, False])
        self.assertEqual(filt.all_zero_selected.tolist(), [False, True])

    def test_duplicate_sample_header_never_silently_mangled(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write(temp, "Geneid,gene_name,A,A\ng,S,4,2\n")
            with self.assertRaisesRegex(AnalysisError, "duplicate"):
                read_counts(path)

    def test_extra_csv_field_cannot_become_an_implicit_row_index(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write(temp, "Geneid,gene_name,A,B\nextra,g,S,4,2\n")
            with self.assertRaises(AnalysisError):
                read_counts(path)

    def test_duplicate_or_missing_gene_id_rejected(self):
        for rows in ["g,S,4\ng,T,2\n", ",S,4\n"]:
            with self.subTest(rows=rows), tempfile.TemporaryDirectory() as temp:
                path = self.write(temp, "Geneid,gene_name,A\n" + rows)
                with self.assertRaises(AnalysisError):
                    read_counts(path)

    def test_version_collisions_and_duplicate_symbols_not_collapsed(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write(temp, "Geneid,gene_name,A\nENSG00000000001.1,S,4\nENSG00000000001.2,S,2\n")
            data = read_counts(path)
        self.assertEqual(data.counts.shape, (1, 2))
        self.assertEqual(data.features.ensembl_id_status.tolist(), ["ambiguous_versionless_id"] * 2)
        self.assertEqual(versionless_ensembl("ABC.3"), "")
        self.assertEqual(versionless_ensembl("ENSG00000000001.2x"), "")

    def test_invalid_counts_never_rounded_or_imputed(self):
        for values in [[[1, -1]], [[1, 1.2]], [[1, np.nan]], [[1, np.inf]],
                       [[1, ""]], [[1, "NA"]], [[True, False]], [[0, 0]],
                       [[2**53, 1]]]:
            with self.subTest(values=values), self.assertRaises(AnalysisError):
                validate_counts(values)
        np.testing.assert_array_equal(validate_counts([[1.0, 2.0]]), [[1, 2]])

    def test_missing_count_csv_is_not_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            path = self.write(temp, "Geneid,gene_name,A,B\ng,S,4,\n")
            with self.assertRaises(AnalysisError):
                read_counts(path)

    def test_cpm_denominator_contains_ineligible_source_feature(self):
        counts = pd.DataFrame([[10, 90], [20, 180]], columns=["target", "other"])
        np.testing.assert_allclose(log2_cpm(counts)[:, 0], np.log2(1 + 100_000))
        # Subsetting first would incorrectly give 1,000,000 CPM.
        self.assertNotAlmostEqual(log2_cpm(counts)[0, 0], log2_cpm(counts[["target"]])[0, 0])

    def test_fixed_filter_uses_all_selected_samples_without_outcome_labels(self):
        counts = pd.DataFrame(np.ones((31, 3), dtype=int), columns=["pass", "fail", "zero"])
        counts.loc[:29, "pass"] = 10
        counts.loc[:28, "fail"] = 10
        counts.loc[:, "zero"] = 0
        result = expression_filter(counts, **FILTER_SETTINGS)
        self.assertEqual(result.eligible.tolist(), [True, False, False])
        self.assertEqual(result.samples_count_ge_threshold.tolist(), [30, 29, 0])

    def test_ratio_normalization_matches_exp_median_log_ratio(self):
        counts = pd.DataFrame([[2, 8, 0], [8, 2, 7], [4, 4, 0]])
        factors, genes = ratio_size_factors(counts)
        self.assertEqual(genes, 2)
        np.testing.assert_allclose(factors, [1, 1, 1])
        with self.assertRaisesRegex(AnalysisError, "every eligible gene"):
            ratio_size_factors(pd.DataFrame([[3, 0], [0, 5]]))


class TestSampleJoinsAndDesign(unittest.TestCase):
    def setUp(self):
        self.frame = metadata_fixture()
        self.metadata = parse_metadata(json_metadata(self.frame))
        self.counts = pd.DataFrame(10, index=self.frame.index[::-1], columns=["g1", "g2"])

    def align(self, counts=None, metadata=None):
        return align_primary_samples(self.counts if counts is None else counts,
            self.metadata if metadata is None else metadata, ("P1", "P2", "P3"))

    def test_sample_join_reorders_metadata_by_exact_title(self):
        selected, audit = self.align()
        self.assertEqual(selected.index.tolist(), self.counts.index.tolist())
        expected = self.frame.loc[self.counts.index]
        np.testing.assert_array_equal(selected.age, expected.age)
        np.testing.assert_array_equal(selected.disease, expected.disease)
        self.assertTrue(audit.selected.all())

    def test_missing_extra_and_duplicate_sample_ids_fail(self):
        missing = self.counts.iloc[1:]
        extra = pd.concat([self.counts, pd.DataFrame(3, index=["unknown"], columns=self.counts.columns)])
        duplicate = pd.concat([self.counts, self.counts.iloc[:1]])
        for counts in [missing, extra, duplicate]:
            with self.subTest(n=len(counts)), self.assertRaises(AnalysisError):
                self.align(counts)

    def test_excluded_geo_rows_retained_but_not_added_to_primary(self):
        german = self.metadata.iloc[:1].copy()
        german.index = ["GermanUC"]
        german["sample_id"] = "GermanUC"
        german["sample_accession"] = "GSM_EXTRA"
        german["plate"] = "13"
        german["source_disease"] = "UC"
        selected, audit = self.align(metadata=pd.concat([self.metadata, german]))
        self.assertEqual(len(selected), len(self.counts))
        self.assertFalse(audit.loc["GermanUC", "selected"])
        self.assertFalse(audit.loc["GermanUC", "in_count_file"])

    def test_missing_covariates_not_dropped_or_imputed(self):
        for field, value in [("age_source", ""), ("age_source", "nan"), ("sex_source", ""),
                             ("source_disease", "UC"), ("plate", "")]:
            metadata = self.metadata.copy()
            metadata.loc[metadata.index[0], field] = value
            with self.subTest(field=field, value=value), self.assertRaises(AnalysisError):
                self.align(metadata=metadata)

    def test_raw_soft_metadata_parser_matches_json(self):
        rows = []
        for record in json_metadata(self.frame)["GSE177044"]:
            rows += [f"^SAMPLE = {record['sample_accession']}",
                     f"!Sample_geo_accession = {record['sample_accession']}",
                     f"!Sample_title = {record['title'][0]}", "!Sample_series_id = GSE177044"]
            rows += [f"!Sample_characteristics_ch1 = {c}" for c in record["characteristics_ch1"]]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "samples.txt"
            path.write_text("\n".join(rows))
            parsed = read_metadata(path)
        pd.testing.assert_frame_equal(parsed, self.metadata)

    def test_duplicate_metadata_titles_accessions_and_keys_fail(self):
        for change in ["title", "sample_accession", "characteristics_ch1"]:
            document = json_metadata(self.frame)
            if change == "characteristics_ch1":
                document["GSE177044"][0][change].append("Age: 55")
            else:
                document["GSE177044"][1][change] = document["GSE177044"][0][change]
            with self.subTest(change=change), self.assertRaises(AnalysisError):
                parse_metadata(document)

    def test_rank_deficiency_perfect_disease_plate_confounding(self):
        frame = self.frame.copy()
        frame["plate"] = frame.disease
        with self.assertRaisesRegex(AnalysisError, "Rank-deficient"):
            build_design(frame)

    def test_unit_leverage_is_an_explicit_hc3_failure(self):
        x = np.column_stack([np.ones(5), [1, 0, 0, 0, 0]])
        with self.assertRaisesRegex(AnalysisError, "unit-leverage"):
            check_design(x)

    def test_explicit_reference_direction_and_age_center(self):
        design, qc = build_design(self.frame)
        np.testing.assert_array_equal(design.disease_PSC, self.frame.disease.eq("PSC").astype(float))
        self.assertAlmostEqual(float(design.age_centered.mean()), 0)
        self.assertEqual(qc["rank"], design.shape[1])
        self.assertEqual(qc["disease_reference"], "Control")


class TestStatistics(unittest.TestCase):
    def test_vectorized_hc3_matches_statsmodels_including_se_p_and_ci(self):
        rng = np.random.default_rng(191)
        metadata = metadata_fixture()
        design, _ = build_design(metadata)
        beta = rng.normal(size=(design.shape[1], 11))
        beta[1] = np.linspace(-1, 1, 11)
        y = design.to_numpy() @ beta + rng.normal(size=(len(metadata), 11)) * np.linspace(0.3, 2, len(metadata))[:, None]
        actual = fit_ols_hc3(design, y, chunk_size=3)
        for gene in range(y.shape[1]):
            reference = sm.OLS(y[:, gene], design.to_numpy()).fit(cov_type="HC3", use_t=False)
            np.testing.assert_allclose(actual.iloc[gene][["effect_log2_cpm", "hc3_se", "normal_z", "pvalue"]].to_numpy(dtype=float),
                [reference.params[1], reference.bse[1], reference.tvalues[1], reference.pvalues[1]], rtol=2e-11, atol=1e-12)
            np.testing.assert_allclose(actual.iloc[gene][["ci95_lower", "ci95_upper"]].to_numpy(dtype=float), reference.conf_int()[1], rtol=1e-11, atol=1e-12)
        self.assertLess(actual.iloc[0].effect_log2_cpm, 0)
        self.assertGreater(actual.iloc[-1].effect_log2_cpm, 0)

    def test_bh_bonferroni_unsorted_ties_missing_and_fixed_denominators(self):
        p = np.asarray([0.05, np.nan, 0.001, 0.7, 0.05, 0, 1])
        valid = np.isfinite(p)
        for method, reference_name in [("bh", "fdr_bh"), ("bonferroni", "bonferroni")]:
            q = adjust_pvalues(p, method)
            self.assertTrue(np.isnan(q[1]))
            np.testing.assert_allclose(q[valid], multipletests(p[valid], method=reference_name)[1])
        corrected = adjust_pvalues(np.where(valid, p, 1))
        reference = multipletests(np.where(valid, p, 1), method="fdr_bh")[1]
        np.testing.assert_allclose(corrected, reference)
        self.assertEqual(corrected[1], 1)
        fixed = adjust_pvalues(p, "bh", family_size=len(p))
        np.testing.assert_allclose(fixed[valid], corrected[valid])
        self.assertTrue(np.isnan(fixed[1]))
        np.testing.assert_allclose(adjust_pvalues([0.02, 0.2], "bonferroni", family_size=10), [0.2, 1])

    def test_invalid_and_all_missing_pvalues(self):
        for p in [[-0.1], [1.1], [np.inf], [[0.1]]]:
            with self.subTest(p=p), self.assertRaises(AnalysisError):
                adjust_pvalues(p)
        self.assertTrue(np.isnan(adjust_pvalues([np.nan, np.nan])).all())
        with self.assertRaises(AnalysisError):
            adjust_pvalues([0.1, 0.2], family_size=1)

    def test_perfect_fit_is_not_an_invented_significant_measurement(self):
        design, _ = build_design(metadata_fixture())
        result = fit_ols_hc3(design, np.ones((len(design), 2)))
        self.assertTrue(result.pvalue.isna().all())
        self.assertTrue(result.ci95_lower.isna().all())
        self.assertTrue(result.padj_bh.eq(1).all())

    def test_leave_plate_out_sign_summary_keeps_fixed_gene_universe(self):
        rng = np.random.default_rng(84)
        metadata = metadata_fixture()
        design, _ = build_design(metadata)
        beta = np.zeros((design.shape[1], 2))
        beta[1] = [2, -3]
        y = design.to_numpy() @ beta + rng.normal(0, 0.01, (len(metadata), 2))
        full = fit_ols_hc3(design, y)
        results, qc = leave_one_plate_out(metadata, y, ["g1", "g2"], full.effect_log2_cpm)
        self.assertEqual(results.source_gene_id.tolist(), ["g1", "g2"])
        self.assertEqual(results.n_same_nonzero_sign.tolist(), [3, 3])
        self.assertEqual(results.n_opposite_sign.tolist(), [0, 0])
        self.assertEqual(set(qc), {"P1", "P2", "P3"})

    def test_joint_subgroup_model_keeps_distinct_effects_and_shared_controls(self):
        rng = np.random.default_rng(212)
        metadata = metadata_fixture()
        y = (2 * metadata.source_disease.eq("PSCUC").to_numpy() - metadata.source_disease.eq("PSC").to_numpy()
             + 0.03 * metadata.age.to_numpy() + 0.1 * metadata.sex.eq("male").to_numpy()
             + rng.normal(0, 0.01, len(metadata)))[:, None]
        results, qc = fit_subgroup_sensitivity(metadata, y, ["gene"])
        self.assertAlmostEqual(results.loc["gene", "ols_psc_alone_effect_log2_cpm"], -1, delta=0.01)
        self.assertAlmostEqual(results.loc["gene", "ols_pscuc_effect_log2_cpm"], 2, delta=0.01)
        self.assertTrue(qc["shared_control_and_nuisance_coefficients"])
        self.assertEqual(qc["n_samples"], len(metadata))

    def test_nb_simulation_positive_and_negative_contrast_and_full_q_family(self):
        rng = np.random.default_rng(13)
        metadata = metadata_fixture()
        design, _ = build_design(metadata)
        base = np.linspace(40, 300, 100)
        mean = np.broadcast_to(base, (len(metadata), len(base))).copy()
        is_psc = metadata.disease.eq("PSC").to_numpy()
        mean[is_psc, 0] *= 6
        mean[is_psc, 1] /= 6
        counts = rng.negative_binomial(12, 12 / (12 + mean)).astype(np.int64)
        # All other genes provide an adequate normalization basis.
        data = pd.DataFrame(counts, index=metadata.index, columns=[f"gene{i:03d}" for i in range(len(base))])
        results, factors, qc = fit_negative_binomial(data, metadata, design, n_cpus=1)
        self.assertEqual(results.index.tolist(), data.columns.tolist())
        self.assertGreater(results.iloc[0].log2FoldChange, 1)
        self.assertLess(results.iloc[1].log2FoldChange, -1)
        np.testing.assert_allclose(results.padj, adjust_pvalues(results.pvalue.fillna(1)))
        self.assertEqual(qc["bh_family_size"], len(data.columns))
        self.assertTrue(factors.size_factor.gt(0).all())
        self.assertIn(qc["dispersion_trend"], ["parametric", "mean"])
        with self.assertRaisesRegex(AnalysisError, "sample order"):
            fit_negative_binomial(data.iloc[::-1], metadata, design)


class TestFrozenContract(unittest.TestCase):
    def test_synthetic_pipeline_keeps_full_feature_universe_and_no_public_sample_rows(self):
        rng = np.random.default_rng(823)
        metadata = metadata_fixture()
        means = np.linspace(80, 600, 100)
        counts = rng.negative_binomial(10, 10 / (10 + means), size=(len(metadata), len(means)))
        counts[:, -1] = 0  # A deposited all-zero feature must remain in outputs.
        ids = [f"ENSG{i:011d}" for i in range(100)]
        with tempfile.TemporaryDirectory(dir=REPO / "work") as temp:
            base = Path(temp)
            count_path = base / "counts.csv"
            source = pd.DataFrame(counts.T, columns=metadata.index)
            source.insert(0, "gene_name", [f"name{i}" for i in range(100)])
            source.insert(0, "Geneid", ids)
            source.to_csv(count_path, index=False)
            metadata_path = base / "metadata.json"
            metadata_path.write_text(json.dumps(json_metadata(metadata)))
            plan = plan_fixture()
            arm = plan["rnaseq"]
            arm["counts"][0].update(path=str(count_path), sha256=sha256(count_path))
            arm["metadata"].update(path=str(metadata_path), sha256=sha256(metadata_path))
            arm["selection"] = {"plates": ["P1", "P2", "P3"], "metadata_samples": len(metadata),
                                "source_disease_counts": metadata.source_disease.value_counts().to_dict()}
            arm["expected_feature_count"] = 100
            plan_path = base / "plan.json"
            plan_path.write_text(json.dumps(plan))
            # Only this synthetic test bypasses the real cohort-size gate.
            with patch("scripts.analyze_psc_blood_rnaseq.validate_plan", return_value=arm):
                qc = run_analysis(plan_path, base / "output", base / "private", sha256(plan_path))
            primary = pd.read_csv(base / "output/rnaseq-primary.csv")
            self.assertEqual(len(primary), 100)
            self.assertEqual(primary.source_gene_id.tolist(), ids)
            self.assertEqual(primary.iloc[-1].model_status, "excluded_prespecified_count_filter")
            self.assertTrue(pd.isna(primary.iloc[-1].pvalue))
            self.assertEqual(qc["source_features"], 100)
            self.assertFalse(any(c in primary for c in metadata.index))
            self.assertTrue((base / "private/design-matrix.csv").is_file())
            receipt = json.loads((base / "output/rnaseq-output-manifest.json").read_text())
            for name, pin in receipt["outputs"].items():
                self.assertEqual(sha256(base / "output" / name), pin["sha256"])

    def test_pins_reject_missing_sha_changes_and_prefix_files(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "source.txt"
            path.write_text("source")
            spec = {"path": str(path), "sha256": sha256(path), "bytes": 6}
            self.assertEqual(check_pin(spec), path)
            for bad in [{"path": str(path)}, {**spec, "sha256": "0" * 64}, {**spec, "bytes": 7}]:
                with self.assertRaises(AnalysisError):
                    check_pin(bad)
            prefix = path.with_suffix(".prefix")
            prefix.write_text("source")
            with self.assertRaises(AnalysisError):
                check_pin({**spec, "path": str(prefix)})

    def test_exact_frozen_contract_accepts_final_methods(self):
        plan = plan_fixture()
        self.assertIs(validate_plan(plan), plan["rnaseq"])

    def test_not_frozen_or_method_drift_fails_before_count_access(self):
        plan = plan_fixture()
        changes = [lambda p: p.update(status="draft"),
            lambda p: p.update(frozen_at_utc="2999-01-01T00:00:00Z"),
            lambda p: p["rnaseq"]["filter"].update(min_samples=10),
            lambda p: p["rnaseq"]["model"].update(refit_cooks=True),
            lambda p: p["rnaseq"]["software"].update(pydeseq2="0.0.0"),
            lambda p: p["rnaseq"].update(expected_sample_count=296),
            lambda p: p["rnaseq"]["selection"].update(plates=["German"])]
        for change in changes:
            copy_plan = copy.deepcopy(plan)
            change(copy_plan)
            with self.subTest(change=change), self.assertRaises(AnalysisError):
                validate_plan(copy_plan)


if __name__ == "__main__":
    unittest.main()

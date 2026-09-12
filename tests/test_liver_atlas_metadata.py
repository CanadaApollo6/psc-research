import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_liver_atlas_metadata import (  # noqa: E402
    annotation_catalog,
    build_geo_index,
    build_library_rows,
    collection_inventory,
    donor_overlap,
    exclusion_summary,
    parse_title_fields,
    profile_dataset,
    read_geo_table,
    verify_source_file,
)


class LiverAtlasMetadataTests(unittest.TestCase):
    def test_pinned_source_verification_rejects_byte_or_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.bin"
            path.write_bytes(b"source")
            declaration = {"id": "source", "bytes": 6, "sha256": "bad"}
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                verify_source_file(path, declaration)
            declaration["sha256"] = ""
            declaration["bytes"] = 5
            with self.assertRaisesRegex(ValueError, "byte count mismatch"):
                verify_source_file(path, declaration)

    def test_title_parser_is_explicit_and_does_not_replace_raw_fields(self):
        self.assertEqual(
            parse_title_fields("PBC001_SN_3pr"),
            {"donor": "PBC001", "suspension": "single-nucleus", "chemistry": "3-prime"},
        )
        self.assertEqual(
            parse_title_fields("C51_flush_SC_3pr"),
            {"donor": "C51", "suspension": "single-cell", "chemistry": "3-prime"},
        )

    def test_profile_reports_assay_by_disease_and_fixed_exclusions(self):
        obs = pd.DataFrame(
            {
                "donor_id": ["C1", "C1", "P1", "P1"],
                "disease": ["normal", "normal", "primary sclerosing cholangitis", "primary sclerosing cholangitis"],
                "assay": ["10x 5' v1", "10x 5' v1", "10x 5' v2", "10x 5' v2"],
                "author_cell_type": ["Kupffer", "Kupffer--LSEC-Doublet", "Chol", "Chol"],
                "cell_type": ["Kupffer cell", "unknown", "intrahepatic cholangiocyte", "intrahepatic cholangiocyte"],
                "cell_type_ontology_term_id": ["CL:1", "unknown", "CL:2", "CL:2"],
            },
            index=["cell-a", "cell-b", "cell-c", "cell-d"],
        )
        profile = profile_dataset(obs, "sc", shape=(4, 10), expected_observations=4)
        self.assertTrue(profile["observation_count_matches_manifest"])
        self.assertEqual(profile["field_value_counts"]["disease"]["normal"], 2)
        self.assertEqual(profile["disease_by_assay"]["10x 5' v1"]["normal"], 2)
        self.assertEqual(profile["exclusion_summary"]["doublet_or_hybrid_cells"], 1)
        self.assertEqual(profile["exclusion_summary"]["unknown_or_missing_ontology_cells"], 1)
        self.assertEqual(profile["exclusion_summary"]["retained_by_fixed_metadata_rules"], 3)

    def test_geo_crosswalk_retains_source_values_and_does_not_claim_library_identity(self):
        fields = [
            "series", "sample", "title", "donor_label_from_title", "cohort_from_title",
            "assay_from_title", "chemistry_from_title", "source_name", "description", "flags", "source_url",
        ]
        record = [
            "GSE247128", "GSM1", "PBC001_SC_5pr", "PBC001", "PBC", "single-cell", "5-prime",
            "single-cell RNA-seq", "3pr scRNA-seq, 10x Genomics", "Title says 5pr; description says 3pr",
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "geo.csv"
            pd.DataFrame([record], columns=fields).to_csv(path, index=False)
            geo = read_geo_table(path)
        index = build_geo_index(geo)
        obs = pd.DataFrame(
            {
                "donor_id": ["PBC001"],
                "disease": ["primary biliary cholangitis"],
                "assay": ["10x 5' v2"],
                "assay_ontology_term_id": ["EFO:1"],
                "suspension_type": ["cell"],
                "library_uuid": ["lib-1"],
                "sample_uuid": ["sample-1"],
                "is_primary_data": [True],
                "author_cell_type": ["Chol"],
                "cell_type": ["intrahepatic cholangiocyte"],
                "cell_type_ontology_term_id": ["CL:2"],
                "mapped_reference_annotation": ["GENCODE 24"],
                "alignment_software": ["Cell Ranger count v3.1.0"],
                "suspension_dissociation_reagent": ["1% TWEEN 20"],
            }
        )
        row = build_library_rows(obs, "sc", "dataset", "version", "atlas.h5ad", index)[0]
        records = json.loads(row["geo_source_records_json"])
        self.assertEqual(records[0]["description"], "3pr scRNA-seq, 10x Genomics")
        self.assertEqual(records[0]["flags"], "Title says 5pr; description says 3pr")
        self.assertTrue(row["geo_donor_label_match"])
        self.assertFalse(row["geo_library_identity_established"])

    def test_annotation_catalog_keeps_author_and_ontology_labels_as_distinct_fields(self):
        obs = pd.DataFrame(
            {
                "author_cell_type": ["P-Hepato", "P-Hepato", "P-Hepato--Mac"],
                "cell_type": ["periportal region hepatocyte", "periportal region hepatocyte", "unknown"],
                "cell_type_ontology_term_id": ["CL:1", "CL:1", "unknown"],
            }
        )
        catalog = annotation_catalog(obs)
        self.assertEqual(catalog[0]["author_cell_type"], "P-Hepato")
        self.assertEqual(catalog[0]["cell_type"], "periportal region hepatocyte")
        self.assertEqual(catalog[0]["cells"], 2)
        self.assertEqual(catalog[1]["author_cell_type"], "P-Hepato--Mac")
        self.assertEqual(catalog[1]["cell_type_ontology_term_id"], "unknown")

    def test_overlap_is_donor_level_and_observation_ids_are_reported_separately(self):
        sc = pd.DataFrame({"donor_id": ["C58", "P1"], "disease": ["normal", "PSC"]}, index=["sc-1", "shared-label"])
        sn = pd.DataFrame({"donor_id": ["C58", "P1"], "disease": ["normal", "PSC"]}, index=["sn-1", "shared-label"])
        overlap = donor_overlap({"sc": sc, "sn": sn})
        self.assertEqual(overlap["shared_donors"], ["C58", "P1"])
        self.assertEqual(overlap["shared_donor_count"], 2)
        self.assertEqual(overlap["shared_observation_ids"], ["shared-label"])

    def test_collection_inventory_exposes_public_dataset_metadata_without_contact_fields(self):
        collection = {
            "collection_id": "collection",
            "collection_version_id": "version",
            "name": "PSC atlas",
            "published_at": "2024-01-01",
            "revised_at": "2026-01-01",
            "contact_email": "should-not-be-copied@example.org",
            "datasets": [
                {
                    "dataset_id": "dataset",
                    "dataset_version_id": "dataset-version",
                    "title": "All nuclei",
                    "cell_count": 3,
                    "primary_cell_count": 2,
                    "feature_count": 10,
                    "raw_data_location": "raw.X",
                    "schema_version": "7.1.0",
                    "suspension_type": ["nucleus"],
                    "assay": [{"label": "10x 3' v3"}],
                    "disease": [{"label": "normal"}],
                    "donor_id": ["C1"],
                    "cell_type": [{"label": "unknown"}],
                    "assets": [{"url": "https://example.org/data.h5ad"}],
                }
            ],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "collection.json"
            path.write_text(json.dumps(collection))
            inventory = collection_inventory(path, ["dataset"])
        self.assertEqual(inventory["dataset_count"], 1)
        self.assertEqual(inventory["datasets"][0]["asset_urls"], ["https://example.org/data.h5ad"])
        self.assertNotIn("contact_email", inventory)

    def test_missing_required_labels_are_counted_without_filling_them(self):
        obs = pd.DataFrame(
            {
                "donor_id": ["C1", ""],
                "disease": ["normal", "normal"],
                "assay": ["10x", "10x"],
                "author_cell_type": ["Kupffer", None],
                "cell_type": ["Kupffer cell", "Kupffer cell"],
            }
        )
        summary = exclusion_summary(obs)
        self.assertEqual(summary["required_metadata_missing_cells"], 1)
        self.assertTrue(pd.isna(obs.loc[1, "author_cell_type"]))


if __name__ == "__main__":
    unittest.main()

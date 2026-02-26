import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest.mock import patch

from starterator import phams, starterate


def make_snapshot(pham_entries, overall_hash="overall-hash", database_version=7):
    return {
        "database_version": database_version,
        "generated_timestamp": "Sun Mar  8 00:00:00 2026",
        "overall_hash": overall_hash,
        "total_phams": len(pham_entries),
        "phams": pham_entries,
    }


def make_entry(combined_hash, gene_count=2, membership_hash=None, annotations_hash=None, metadata_hash=None):
    return {
        "membership_hash": membership_hash or ("%s-membership" % combined_hash),
        "sequences_hash": "%s-sequences" % combined_hash,
        "annotations_hash": annotations_hash or ("%s-annotations" % combined_hash),
        "metadata_hash": metadata_hash or ("%s-metadata" % combined_hash),
        "combined_hash": combined_hash,
        "gene_count": gene_count,
    }


class TestCompareHashSnapshots(unittest.TestCase):
    def test_detects_added_phams(self):
        previous = make_snapshot({"1": make_entry("a")}, overall_hash="before")
        current = make_snapshot({"1": make_entry("a"), "2": make_entry("b")}, overall_hash="after")

        result = phams.compare_hash_snapshots(previous, current)

        self.assertFalse(result["overall_hash_matches"])
        self.assertCountEqual(result["phams_added"], ["2"])
        self.assertEqual(result["phams_modified"], [])
        self.assertEqual(result["phams_removed"], [])

    def test_detects_modified_phams_using_full_entry_dict(self):
        previous = make_snapshot({1: make_entry("same-hash", gene_count=2)}, overall_hash="before")
        current = make_snapshot({1: make_entry("same-hash", gene_count=3)}, overall_hash="after")

        result = phams.compare_hash_snapshots(previous, current)

        self.assertEqual(result["phams_added"], [])
        self.assertEqual(result["phams_removed"], [])
        self.assertEqual(result["phams_modified"], ["1"])

    def test_detects_removed_phams(self):
        previous = make_snapshot({"1": make_entry("a"), "2": make_entry("b")}, overall_hash="before")
        current = make_snapshot({"2": make_entry("b")}, overall_hash="after")

        result = phams.compare_hash_snapshots(previous, current)

        self.assertEqual(result["phams_added"], [])
        self.assertEqual(result["phams_modified"], [])
        self.assertCountEqual(result["phams_removed"], ["1"])

    def test_identical_snapshots_yield_no_changes(self):
        snapshot = make_snapshot({"1": make_entry("a"), "2": make_entry("b")}, overall_hash="same")

        result = phams.compare_hash_snapshots(snapshot, snapshot)

        self.assertTrue(result["overall_hash_matches"])
        self.assertEqual(result["phams_added"], [])
        self.assertEqual(result["phams_modified"], [])
        self.assertEqual(result["phams_removed"], [])


class _HashFakeDB:
    def __init__(self, pham_ids):
        self._pham_ids = pham_ids

    def query(self, query, params=None):
        normalized_query = " ".join(query.lower().split())
        if normalized_query == "select version from version;":
            return [(7,)]
        if normalized_query == "select phamid from pham order by phamid;":
            return [(pham_id,) for pham_id in self._pham_ids]
        if "from pham join gene on pham.phamid = gene.phamid" in normalized_query:
            requested_ids = tuple(params[0])
            rows = []
            for pham_id in requested_ids:
                rows.append(
                    (
                        "Phage-%s" % pham_id,
                        "A",
                        "A1",
                        "final",
                        pham_id,
                        "MSEQ-%s" % pham_id,
                        "gene-%s" % pham_id,
                        "LT-%s" % pham_id,
                        "note-%s" % pham_id,
                        "Gene-%s" % pham_id,
                    )
                )
            return rows
        raise AssertionError("Unexpected query: %s" % query)


class TestGeneratePhamHashesDeterminism(unittest.TestCase):
    def test_overall_hash_is_stable_for_same_ordered_db_content(self):
        with patch.object(phams, "get_db", return_value=_HashFakeDB([2, 1])):
            first = phams.generate_pham_hashes(output_file=False)
        with patch.object(phams, "get_db", return_value=_HashFakeDB([1, 2])):
            second = phams.generate_pham_hashes(output_file=False)

        self.assertEqual(first["overall_hash"], second["overall_hash"])
        self.assertEqual(first["phams"], second["phams"])


class TestDeltaPhamCliValidation(unittest.TestCase):
    def test_delta_phams_conflicts_with_all_phams(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as context:
            starterate.get_arguments(["--delta-phams", "previous.json", "--all-phams"])

        self.assertEqual(context.exception.code, 2)
        self.assertIn("mutually exclusive", stderr.getvalue())

    def test_delta_phams_conflicts_with_phage_mode(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as context:
            starterate.get_arguments(["--delta-phams", "previous.json", "--phage", "Demo"])

        self.assertEqual(context.exception.code, 2)
        self.assertIn("mutually exclusive", stderr.getvalue())

    def test_delta_phams_conflicts_with_single_pham_mode(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as context:
            starterate.get_arguments(["--delta-phams", "previous.json", "--pham_no", "123"])

        self.assertEqual(context.exception.code, 2)
        self.assertIn("mutually exclusive", stderr.getvalue())


class TestDeltaPhamMode(unittest.TestCase):
    def _generate_stub(self, output_dir, current_snapshot):
        def _generate(output_file=True):
            if output_file:
                output_path = os.path.join(
                    output_dir,
                    "pham_hashes_v%s.json" % current_snapshot["database_version"],
                )
                with open(output_path, "w") as handle:
                    json.dump(current_snapshot, handle)
            return current_snapshot

        return _generate

    def test_invalid_previous_hash_file_fails_cleanly(self):
        with patch.object(starterate.utils, "get_config", return_value={"count": 0}), \
                patch.object(starterate.phamgene, "check_protein_db"), \
                patch.object(starterate.logging, "error") as mock_error:
            with self.assertRaises(SystemExit) as context:
                starterate.main(["--delta-phams", "/tmp/does-not-exist.json"])

        self.assertEqual(context.exception.code, 1)
        self.assertIn("Could not read previous pham hash file", mock_error.call_args[0][0])

    def test_invalid_previous_hash_json_fails_cleanly(self):
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write("{not json")
            previous_path = handle.name
        self.addCleanup(lambda: os.path.exists(previous_path) and os.remove(previous_path))

        with patch.object(starterate.utils, "get_config", return_value={"count": 0}), \
                patch.object(starterate.phamgene, "check_protein_db"), \
                patch.object(starterate.logging, "error") as mock_error:
            with self.assertRaises(SystemExit) as context:
                starterate.main(["--delta-phams", previous_path])

        self.assertEqual(context.exception.code, 1)
        self.assertIn("Invalid JSON in previous pham hash file", mock_error.call_args[0][0])

    def test_noop_when_no_eligible_changed_phams(self):
        previous = make_snapshot({"1": make_entry("a")}, overall_hash="same")
        current = make_snapshot({"1": make_entry("a")}, overall_hash="same")

        with tempfile.TemporaryDirectory() as temp_dir:
            previous_path = os.path.join(temp_dir, "previous.json")
            with open(previous_path, "w") as handle:
                json.dump(previous, handle)

            stdout = io.StringIO()
            with patch.object(starterate.utils, "get_config", return_value={"count": 0}), \
                    patch.object(starterate.utils, "FINAL_DIR", temp_dir), \
                    patch.object(starterate.phamgene, "check_protein_db"), \
                    patch.object(starterate, "generate_pham_hashes", side_effect=self._generate_stub(temp_dir, current)), \
                    patch.object(starterate, "get_all_phams", return_value=[3, 2, 1]), \
                    patch.object(starterate, "process_pham_list") as mock_process, \
                    redirect_stdout(stdout):
                starterate.main(["--delta-phams", previous_path])

            self.assertFalse(mock_process.called)
            self.assertIn("No eligible changed phams found; nothing to process.", stdout.getvalue())
            self.assertTrue(os.path.exists(os.path.join(temp_dir, "pham_hashes_v7.json")))

    def test_only_changed_and_eligible_phams_are_processed_in_get_all_phams_order(self):
        previous = make_snapshot(
            {
                "1": make_entry("unchanged"),
                "2": make_entry("same-combined", gene_count=2),
                "3": make_entry("removed"),
            },
            overall_hash="before",
        )
        current = make_snapshot(
            {
                "1": make_entry("unchanged"),
                "2": make_entry("same-combined", gene_count=3),
                "4": make_entry("added"),
                "5": make_entry("ineligible-added"),
            },
            overall_hash="after",
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            previous_path = os.path.join(temp_dir, "previous.json")
            with open(previous_path, "w") as handle:
                json.dump(previous, handle)

            with patch.object(starterate.utils, "get_config", return_value={"count": 0}), \
                    patch.object(starterate.utils, "FINAL_DIR", temp_dir), \
                    patch.object(starterate.phamgene, "check_protein_db"), \
                    patch.object(starterate, "generate_pham_hashes", side_effect=self._generate_stub(temp_dir, current)), \
                    patch.object(starterate, "get_all_phams", return_value=[4, 1, 2, 6]), \
                    patch.object(starterate, "process_pham_list") as mock_process:
                starterate.main([
                    "--delta-phams",
                    previous_path,
                    "--no-pdfs",
                    "--compress-json",
                ])

            mock_process.assert_called_once_with([4, 2], no_pdfs=True, compress_json=True)
            self.assertTrue(os.path.exists(os.path.join(temp_dir, "pham_hashes_v7.json")))


if __name__ == "__main__":
    unittest.main()

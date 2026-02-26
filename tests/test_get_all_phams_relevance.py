import unittest
from unittest.mock import patch

from starterator import phams


class _FakeGene:
    def __init__(self, gene_id, is_valid):
        self.gene_id = gene_id
        self._is_valid = is_valid

    def has_valid_start(self):
        return self._is_valid


class _FakePhage:
    def __init__(self, sequence):
        self._sequence = sequence

    def get_name(self):
        return "fake_phage"

    def get_sequence(self):
        return self._sequence


class _FakeDB:
    def __init__(self, candidate_rows, gene_rows_by_pham):
        self._candidate_rows = candidate_rows
        self._gene_rows_by_pham = gene_rows_by_pham

    def query(self, query, params=None):
        normalized_query = " ".join(query.lower().split())
        if "group by gene.phamid" in normalized_query:
            return self._candidate_rows
        if isinstance(params, tuple):
            pham_no = params[0]
        else:
            pham_no = params
        return self._gene_rows_by_pham.get(pham_no, [])


class TestGetAllPhamsRelevance(unittest.TestCase):
    def _run_get_all_phams(self, candidate_rows, gene_rows_by_pham, valid_start_by_gene, sequence_by_phage):
        fake_db = _FakeDB(candidate_rows, gene_rows_by_pham)

        def fake_new_pham_gene(gene_id, start, stop, orientation, phage_id, name):
            return _FakeGene(gene_id, valid_start_by_gene[gene_id])

        def fake_new_phage(phage_id=None, **_kwargs):
            return _FakePhage(sequence_by_phage[phage_id])

        with patch.object(phams, "get_db", return_value=fake_db), \
                patch.object(phams, "new_PhamGene", side_effect=fake_new_pham_gene), \
                patch.object(phams, "new_phage", side_effect=fake_new_phage):
            return phams.get_all_phams()

    def test_includes_pham_with_two_valid_genes(self):
        relevant = self._run_get_all_phams(
            candidate_rows=[(100, 2)],
            gene_rows_by_pham={
                100: [
                    ("g100a", "pA", 0, 1, 30, "F", "1"),
                    ("g100b", "pB", 0, 31, 60, "F", "2"),
                ],
            },
            valid_start_by_gene={"g100a": True, "g100b": True},
            sequence_by_phage={"pA": "ATGATGATG", "pB": "ATGATGATG"},
        )
        self.assertEqual(relevant, [100])

    def test_excludes_pham_when_all_genes_are_filtered_after_unknown_status(self):
        # Simulate unknown-only pham: status filter in SQL removes all rows.
        relevant = self._run_get_all_phams(
            candidate_rows=[(400, 2)],
            gene_rows_by_pham={400: []},
            valid_start_by_gene={},
            sequence_by_phage={},
        )
        self.assertEqual(relevant, [])

    def test_excludes_pham_when_genome_contains_n(self):
        relevant = self._run_get_all_phams(
            candidate_rows=[(200, 2)],
            gene_rows_by_pham={
                200: [
                    ("g200a", "pN", 0, 1, 30, "F", "1"),
                    ("g200b", "pN", 0, 31, 60, "F", "2"),
                ],
            },
            valid_start_by_gene={"g200a": True, "g200b": True},
            sequence_by_phage={"pN": "ATGNATGATG"},
        )
        self.assertEqual(relevant, [])

    def test_excludes_pham_with_only_one_valid_gene(self):
        relevant = self._run_get_all_phams(
            candidate_rows=[(300, 2)],
            gene_rows_by_pham={
                300: [
                    ("g300a", "pC", 0, 1, 30, "F", "1"),
                    ("g300b", "pD", 0, 31, 60, "F", "2"),
                ],
            },
            valid_start_by_gene={"g300a": True, "g300b": False},
            sequence_by_phage={"pC": "ATGATGATG", "pD": "ATGATGATG"},
        )
        self.assertEqual(relevant, [])

    def test_includes_mixed_pham_when_two_genes_still_valid(self):
        relevant = self._run_get_all_phams(
            candidate_rows=[(500, 3)],
            gene_rows_by_pham={
                500: [
                    ("g500a", "pE", 0, 1, 30, "F", "1"),
                    ("g500b", "pE", 0, 31, 60, "F", "2"),
                    ("g500c", "pF", 0, 61, 90, "F", "3"),
                ],
            },
            valid_start_by_gene={"g500a": True, "g500b": False, "g500c": True},
            sequence_by_phage={"pE": "ATGATGATG", "pF": "ATGATGATG"},
        )
        self.assertEqual(relevant, [500])

    def test_get_all_phams_filters_with_stable_order(self):
        # Candidate order from DB query should be preserved after filtering.
        relevant = self._run_get_all_phams(
            candidate_rows=[
                (500, 8),
                (100, 5),
                (200, 4),
                (300, 3),
                (400, 2),
            ],
            gene_rows_by_pham={
                500: [
                    ("g500a", "pE", 0, 1, 30, "F", "1"),
                    ("g500b", "pE", 0, 31, 60, "F", "2"),
                    ("g500c", "pF", 0, 61, 90, "F", "3"),
                ],
                100: [
                    ("g100a", "pA", 0, 1, 30, "F", "1"),
                    ("g100b", "pB", 0, 31, 60, "F", "2"),
                ],
                200: [
                    ("g200a", "pN", 0, 1, 30, "F", "1"),
                    ("g200b", "pN", 0, 31, 60, "F", "2"),
                ],
                300: [
                    ("g300a", "pC", 0, 1, 30, "F", "1"),
                    ("g300b", "pD", 0, 31, 60, "F", "2"),
                ],
                400: [],
            },
            valid_start_by_gene={
                "g500a": True,
                "g500b": False,
                "g500c": True,
                "g100a": True,
                "g100b": True,
                "g200a": True,
                "g200b": True,
                "g300a": True,
                "g300b": False,
            },
            sequence_by_phage={
                "pE": "ATGATGATG",
                "pF": "ATGATGATG",
                "pA": "ATGATGATG",
                "pB": "ATGATGATG",
                "pN": "ATGNATGATG",
                "pC": "ATGATGATG",
                "pD": "ATGATGATG",
            },
        )

        self.assertEqual(relevant, [500, 100])
        self.assertNotIn(200, relevant)
        self.assertNotIn(300, relevant)
        self.assertNotIn(400, relevant)

    def test_pham_get_genes_keeps_existing_failed_validation_error(self):
        fake_validation = {
            "genes": {},
            "total_count": 2,
            "skipped_n": ["g1"],
            "skipped_start": ["g2"],
        }

        pham = object.__new__(phams.Pham)
        pham.pham_no = 230
        pham.count = 0

        with patch.object(phams, "_validated_genes_for_pham", return_value=fake_validation):
            with self.assertRaises(phams.StarteratorError) as context:
                phams.Pham.get_genes(pham)

        message = str(context.exception)
        self.assertIn("Pham 230 has 2 gene(s) but all failed validation.", message)
        self.assertIn("skipped due to N's in genome", message)
        self.assertIn("skipped due to invalid start codon", message)


if __name__ == "__main__":
    unittest.main()

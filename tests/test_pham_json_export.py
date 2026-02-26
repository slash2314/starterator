import unittest
from types import SimpleNamespace
from unittest.mock import patch

from starterator import phams


def build_stub_gene():
    return SimpleNamespace(
        gene_id="DemoGene_1",
        db_id=None,
        start=100,
        stop=160,
        orientation="F",
        draftStatus=False,
        locustag="DEMO_LOCUS_1",
        annot_author=1,
        alignment_candidate_start_nums=[1, 2],
        alignment_candidate_starts=[4, 7],
        alignment_feature_runs=((0, 2, "seq"), (2, 4, "gap"), (4, 7, "seq")),
        alignment_start_num_called=2,
        alignment_candidate_start_counts=[1, 1],
        alignment_start_conservation=[1.0, 1.0],
        alignment_annot_start_nums=[2],
        alignment_annot_start_counts=[1],
        calls_most_annotated=True,
        has_most_annotated=True,
        called_start_is_bad=False,
        bad_adjacent_start_nums=[],
        cluster="A",
        subcluster="A1",
        alignment="ABCDEFG",
        alignment_indices_to_coords_optimized=lambda indices: [105, 108],
        get_locustag=lambda: None,
    )


class TestPhamJsonExport(unittest.TestCase):
    def build_stub_pham(self, gene):
        pham = object.__new__(phams.Pham)
        pham.pham_no = "123"
        pham.count = 1
        pham.aligner = "MAFFT"
        pham.total_possible_starts = [4, 7]
        pham.genes = {gene.gene_id: gene}
        pham.stats = {
            "most_common": {
                "annot_list": [gene],
                "annot_counts": {2: 1},
                "possible": {1: [gene.gene_id], 2: [gene.gene_id]},
                "most_called_start": 2,
                "most_annotated_start": 2,
            }
        }
        return pham

    def test_annot_summary_exports_seq_runs_and_drops_legacy_fields(self):
        gene = build_stub_gene()
        pham = self.build_stub_pham(gene)

        with patch.object(phams, "get_version", return_value=99), \
                patch.object(phams.Pham, "group_similar_genes", return_value=[[gene]]):
            summary = pham.annot_summary()

        exported_gene = summary["Genes"][0]

        self.assertEqual(exported_gene["GeneID"], "DemoGene_1")
        self.assertEqual(exported_gene["CalledStartAlignmentIndex"], 7)
        self.assertIn("SeqRuns", exported_gene)
        self.assertEqual(exported_gene["SeqRuns"], [[0, 2], [4, 7]])

        for run in exported_gene["SeqRuns"]:
            self.assertIsInstance(run, list)
            self.assertEqual(len(run), 2)

        self.assertNotIn("FeatureRuns", exported_gene)
        self.assertNotIn("StartConservation", exported_gene)
        self.assertNotIn("StartConservationCounts", exported_gene)


if __name__ == "__main__":
    unittest.main()

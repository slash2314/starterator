import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from starterator import phamgene


def build_gene(alignment_sequence, candidate_starts, ahead_of_start, orientation="F", start=100, genome_length=1000):
    gene = object.__new__(phamgene.PhamGene)
    gene.alignment = SeqRecord(Seq(alignment_sequence), id="test_gene")
    gene.alignment.features = []
    gene.candidate_starts = list(candidate_starts)
    gene.ahead_of_start = ahead_of_start
    gene.orientation = orientation
    gene.start = start
    gene.genome_length = genome_length
    gene.alignment_start_site = ahead_of_start
    return gene


def legacy_compute_alignment_start_and_candidates(alignment_sequence, ahead_of_start, candidate_starts):
    start_count = -1
    start_site = 0
    start_found = False
    candidate_count = -1
    aligned_starts = []
    candidate_start_set = set(candidate_starts)

    for index, char in enumerate(alignment_sequence):
        if not start_found and char in "AGTC":
            start_count += 1
            start_site = index
            if start_count >= ahead_of_start:
                start_found = True
        if char != "-":
            candidate_count += 1
            if candidate_count in candidate_start_set:
                aligned_starts.append(index)

    if start_count > ahead_of_start:
        start_site -= 1
    return start_site, aligned_starts


def legacy_non_gap_prefix(alignment_sequence):
    prefix = [0] * (len(alignment_sequence) + 1)
    count = 0
    for index, nt in enumerate(alignment_sequence, 1):
        if nt != "-":
            count += 1
        prefix[index] = count
    return prefix


def legacy_alignment_coord(alignment_sequence, index, start, ahead_of_start, orientation, genome_length):
    index = int(index)
    if index < 0:
        index = 0
    if index > len(alignment_sequence):
        index = len(alignment_sequence)

    non_gap_count = legacy_non_gap_prefix(alignment_sequence)[index]
    if orientation == "R":
        return (start + ahead_of_start - non_gap_count - 1) % genome_length + 1
    return start - ahead_of_start + non_gap_count + 1


def legacy_feature_runs(alignment_sequence, start_boundary):
    sequence_len = len(alignment_sequence)
    if sequence_len == 0:
        return ()

    segment_start = 0
    segment_type = "seq" if alignment_sequence[0] in "ACGT" else "gap"
    run_list = []

    for index in range(1, sequence_len + 1):
        at_end = index == sequence_len
        at_boundary = index == start_boundary and index != segment_start
        if not at_end:
            current_type = "seq" if alignment_sequence[index] in "ACGT" else "gap"
        else:
            current_type = None
        type_change = (not at_end) and (current_type != segment_type)

        if at_end or at_boundary or type_change:
            run_list.append((segment_start, index, segment_type))
            if at_end:
                break
            segment_start = index
            segment_type = current_type
    return tuple(run_list)


class TestPhamGeneAlignmentCache(unittest.TestCase):
    def test_compute_alignment_start_and_candidates_matches_legacy(self):
        cases = [
            ("A--CGTTA", [0, 2, 4, 6], 0),
            ("--A-CN-GT-", [0, 1, 3, 4], 1),
            ("NN--AAGT", [0, 2, 3, 6], 3),
            ("----", [0, 1], 2),
            ("A--CGTTA", [0, 3, 6, 99], -1),
        ]

        for alignment_sequence, candidate_starts, ahead in cases:
            gene = build_gene(alignment_sequence, candidate_starts, ahead)
            expected_start, expected_candidates = legacy_compute_alignment_start_and_candidates(
                alignment_sequence, ahead, candidate_starts
            )
            actual_start, actual_candidates = gene._compute_alignment_start_and_candidates()
            self.assertEqual(expected_start, actual_start)
            self.assertEqual(expected_candidates, actual_candidates)

            # Repeat call to exercise cached path and ensure identical output.
            cached_start, cached_candidates = gene._compute_alignment_start_and_candidates()
            self.assertEqual(expected_start, cached_start)
            self.assertEqual(expected_candidates, cached_candidates)

    def test_get_alignment_coord_lookup_matches_legacy_prefix(self):
        gene = build_gene("A--CG-TN", [0, 2], 1)
        seq_len, prefix, _coord_cache = gene._get_alignment_coord_lookup()
        self.assertEqual(seq_len, 8)
        self.assertEqual(prefix, legacy_non_gap_prefix("A--CG-TN"))

        # Reusing the same alignment should return the same cached prefix object.
        _seq_len2, prefix2, _coord_cache2 = gene._get_alignment_coord_lookup()
        self.assertIs(prefix, prefix2)

        # Updating the alignment sequence should rebuild the lookup.
        gene.alignment.seq = Seq("ACGT")
        seq_len3, prefix3, _coord_cache3 = gene._get_alignment_coord_lookup()
        self.assertEqual(seq_len3, 4)
        self.assertEqual(prefix3, legacy_non_gap_prefix("ACGT"))

    def test_alignment_coord_optimized_matches_legacy(self):
        alignment_sequence = "A-CG--TTA"
        indices = [-5, 0, 1, 3, 7, 9, 20]

        forward_gene = build_gene(
            alignment_sequence,
            [0, 2],
            ahead_of_start=2,
            orientation="F",
            start=100,
            genome_length=1000,
        )
        reverse_gene = build_gene(
            alignment_sequence,
            [0, 2],
            ahead_of_start=2,
            orientation="R",
            start=100,
            genome_length=1000,
        )

        for gene in (forward_gene, reverse_gene):
            expected = [
                legacy_alignment_coord(
                    alignment_sequence,
                    index,
                    gene.start,
                    gene.ahead_of_start,
                    gene.orientation,
                    gene.genome_length,
                )
                for index in indices
            ]
            actual = [gene.alignment_index_to_coord_optimized(index) for index in indices]
            actual_batch = gene.alignment_indices_to_coords_optimized(indices)
            self.assertEqual(expected, actual)
            self.assertEqual(expected, actual_batch)

    def test_add_gaps_as_features_matches_legacy_runs(self):
        alignment_sequence = "A--CGTT--A"
        boundaries = [0, 3, 5, len(alignment_sequence)]
        template_cache = {}

        for boundary in boundaries:
            gene = build_gene(alignment_sequence, [0, 2], ahead_of_start=boundary)
            gene.alignment_start_site = boundary
            gene.add_gaps_as_features(feature_template_cache=template_cache)
            self.assertEqual(gene.alignment_feature_runs, legacy_feature_runs(alignment_sequence, boundary))
            self.assertEqual(gene.alignment.features, [])

        # Reuse template cache with a new object and ensure parity remains.
        second_gene = build_gene(alignment_sequence, [0, 2], ahead_of_start=3)
        second_gene.alignment_start_site = 3
        second_gene.add_gaps_as_features(feature_template_cache=template_cache)
        self.assertEqual(second_gene.alignment_feature_runs, legacy_feature_runs(alignment_sequence, 3))

    def test_add_gaps_as_features_lowercase_handling_matches_legacy(self):
        alignment_sequence = "Aa--cG"
        boundary = 2
        gene = build_gene(alignment_sequence, [0], ahead_of_start=boundary)
        gene.alignment_start_site = boundary
        gene.add_gaps_as_features()
        self.assertEqual(gene.alignment_feature_runs, legacy_feature_runs(alignment_sequence, boundary))


if __name__ == "__main__":
    unittest.main()

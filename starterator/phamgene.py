# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Class and functions for pham genes

from .phage import new_phage
from .database import DB, get_db
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
from . import utils
from .utils import StarteratorError, clean_up_files
import subprocess
import math
import os
import numpy as np


def get_protein_sequences():
    proteins = []
    results = get_db().query('SELECT GeneID, Translation from gene')
    for row in results:
        gene_id = row[0].replace("-", "_")
        translation = utils.decode_if_bytes(row[1])
        protein = SeqRecord(Seq(translation.replace('-', '')),
                            id=gene_id+"_", name=row[0], description=gene_id)
        proteins.append(protein)
    return proteins


def update_protein_db():
    clean_up_files(utils.INTERMEDIATE_DIR)
    proteins = get_protein_sequences()
    try:
        fasta_file = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")
        SeqIO.write(proteins, fasta_file, 'fasta')
    except:
        print("creating proteins folder in correct place")
        utils.create_folders()
        fasta_file = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")
        SeqIO.write(proteins, fasta_file, 'fasta')

    blast_db_command =  [
           'makeblastdb',
           '-in', fasta_file,
           '-dbtype', 'prot',
           '-title', 'Proteins',
           '-out', fasta_file
    ]
    # print blast_db_command
    # else:
    #     blast_db_command = [BLAST_DIR + 'formatdb',
    #                 '-i', "\""+ fasta_file+ "\"",
    #                 '-o', 'T',
    #                 "-t", "Proteins"]
    #     print blast_db_command
    subprocess.check_call(blast_db_command)


def check_protein_db(count):
    results = get_db().query('SELECT count(*) from gene')
    new_count = results[0][0]
    # print new_count
    if int(new_count) != int(count):
        update_protein_db()
        config = utils.get_config()
        config["count"] = new_count
        utils.write_to_config_file(config)


def get_pham_no(phage_name, gene_number):
    """
        Gets the pham number of a gene, given the phage name and the gene number
    """
    db = DB()
    gene_number = str(gene_number)
    gene_id_pattern = re.compile(r'^([A-Za-z0-9]*_)*([A-Za-z])*%s$' % re.escape(gene_number))

    def first_match(rows):
        for row in rows:
            if len(row) < 2:
                continue
            pham_id = row[0]
            gene_id = row[1]
            if pham_id is not None and gene_id and gene_id_pattern.match(gene_id):
                return str(pham_id)
        return None

    try:
        # First pass: constrain by phage and select matching gene ids in Python.
        results = db.query(
            "SELECT gene.PhamID, gene.GeneID \n\
             FROM gene JOIN phage ON gene.PhageID = phage.PhageID \n\
             WHERE (phage.Name LIKE %s or phage.PhageID = %s)",
            (phage_name + "%", phage_name),
        )
        pham_no = first_match(results)
        if pham_no is not None:
            return pham_no

        # Fallback pass: match by root of gene id naming.
        results = db.query(
            "SELECT gene.PhamID, gene.GeneID \n\
             FROM gene JOIN phage ON gene.PhageID = phage.PhageID \n\
             WHERE gene.GeneID LIKE %s",
            phage_name + "%",
        )
        pham_no = first_match(results)
        if pham_no is not None:
            return pham_no
    except StarteratorError:
        raise
    except Exception:
        raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))

    raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))


def find_upstream_stop_site(start, stop, orientation, phage_sequence):
    """
        Given the coordinates of a gene, the sequence of the phage it is in, and 
        the orientation of the gene, returns a sequence that contains the gene
        and upstream sequence before a stop site.
    """
    ahead_of_start = 0
    stop_site_found = False
    stop_codons = {'AGT', 'AAT', 'GAT'}
    while not stop_site_found:
        ahead_of_start += 99
        if orientation == 'R':
            if start + ahead_of_start > len(phage_sequence):     # i.e. hit end of phage while looking for stop
                ahead_of_start = len(phage_sequence) - start   # start is zero based counting
                ahead_of_start = ahead_of_start - ahead_of_start % 3
                sequence = Seq(phage_sequence[stop:(start+ahead_of_start)])
                sequence = sequence.reverse_complement()
                return sequence, ahead_of_start

            end = start + ahead_of_start
            if stop <= end:
                frag = phage_sequence[stop:end]
            else:
                # Reverse-strand circular wrap-around: gene crosses the origin
                frag = phage_sequence[stop:] + phage_sequence[:end]
            sequence = Seq(frag).reverse_complement()
        else:
            if start < ahead_of_start:
                ahead_of_start = start - start % 3
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop])
                return sequence, ahead_of_start
            if stop < start:
                end_sequence = phage_sequence[(start-ahead_of_start):]
                start_sequence = phage_sequence[:stop]
                sequence = Seq(end_sequence+start_sequence)
            else:
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop])
        sequence_ahead_of_start = sequence[:ahead_of_start]
        sequence_ahead_of_start = sequence_ahead_of_start[::-1]
        upstream_seq_str = str(sequence_ahead_of_start)
        for index in range(0, len(upstream_seq_str), 3):
            codon = upstream_seq_str[index:index+3]
            if codon in stop_codons:
                new_ahead_of_start = index
                new_sequence = sequence[(ahead_of_start - index):]
                return new_sequence, new_ahead_of_start


class Gene(object):
    def __init__(self, phage, name, start, stop, orientation, db_id=None):
        self.phage = phage
        self.name = name
        self.start = start
        self.stop = stop
        self.orientation = orientation
        self.db_id = db_id
   
    def gene_no(self):
        get_gene_number(self.name)


pham_genes = {}

PRINTED_BAD_STARTS = set()


def new_PhamGene(db_id, start, stop, orientation, phage_id, name, phage_sequence=None):
    if db_id is None:
        return UnPhamGene(db_id, start, stop, orientation, phage_id, phage_sequence)
    if pham_genes.get(db_id, True):
        pham_genes[db_id] = PhamGene(db_id, start, stop, orientation, phage_id, name)
    return pham_genes[db_id]


def get_gene_number(gene_name):
    """ Given a gene_name, returns the number of the gene
    """
    # NAMING IN THIS DATABASE DRIVES ME CRAZY!!!
    # GeneID in database: form of <PhageID>_(<PhageName>([_-]Draft*))*_(gene)*(gp)*<GeneNo>
    #   where PhageID can be a number or the name of the phage (with _Draft perhaps)
    
    # ...and I don't think this function is needed. Precisely because of this!
    # match = re.search(r'^(\w+)([_]*\w*)_([])')
    match = re.search(r'^((\w+)([_-]*\w*)_)*([a-zA-Z]*)([0-9]+)+$', gene_name)
    gene_number = match.groups()[-1]
    try:
        return int(gene_number)
    except:
        # it is one of the 3 horrible genes that do not have an actual number
        # !!! WHAT DO I DO HERE?????  doesn't look like it is terrible for it to be a string?
        # so it will return hypothetical or null - don't ask me why
        return gene_number


class PhamGene(Gene):
    def __init__(self, db_id, start, stop, orientation, phage_id, name, pham_no=None):
        self.db_id = db_id
        self.gene_id = db_id
        self.phage_id = phage_id
        self.start = start
        self.stop = stop
        self.cluster = None
        self.subcluster = None
        self.cluster_hash = None
        self.locustag = None
        self.annot_author = None
        self.gene_no = name
        self.full_name = self.phage_id + "_" + self.gene_no

        self.status = None # 'draft' = auto-annotated, 'final' = final/approved, 'gbk' imported non Pitt phage

        if orientation == 'R':
            self.start_codon_location = stop
            self.stop_codon_location = start + 1
        else:
            self.start_codon_location = start + 1
            self.stop_codon_location = stop

        self.orientation = orientation
        self.pham_no = pham_no
        self.pham_size = None
        # self.translation
        self.ahead_of_start = None
        self.sequence = self.make_gene()
        self.candidate_starts = self.add_candidate_starts()

        #adjacent-start clusters and "bad starts" (all-but-last in each cluster)
        self.adjacent_candidate_start_groups = self._find_adjacent_start_groups()

        # Flatten clusters into a single "bad starts" list:
        # for each adjacent run [a,b,c], treat [a,b] as bad (drop last), then merge across runs.
        bad = []
        for grp in self.adjacent_candidate_start_groups:
            if len(grp) >= 2:
                bad.extend(grp[:-1])
        self.bad_adjacent_candidate_starts = sorted(set(bad))
        self.has_bad_adjacent_candidate_starts = bool(self.bad_adjacent_candidate_starts)


        '''
        if self.has_bad_adjacent_candidate_starts:
            print(
                f"[QC] bad adjacent starts (offsets) gene={getattr(self, 'gene_no', getattr(self, 'number', '?'))}: {self.bad_adjacent_candidate_starts}")
        '''


        self.alignment = None
        self.alignment_start_site = None
        self.alignment_candidate_starts = None
        self.alignment_feature_runs = None
        self.alignment_candidate_start_nums = None
        self.alignment_candidate_start_counts = None
        self.alignment_annot_start_nums = None
        self.alignment_annot_start_counts = None
        self.alignment_annot_start_fraction = None
        self.alignment_annot_counts_by_start = {}
        self.alignment_start_num_called = None
        self.alignment_start_conservation = None
        self.calls_most_annotated = None
        self.has_most_annotated = None
        self.suggested_start = {}

    def make_gene(self):
        """
           makes the gene which is a SeqRecord from Biopython. In this case the "gene" should
           include all the sequence upstream of the annotated start all the way to the first
           in frame stop codon.
        """
        phage = new_phage(phage_id=self.phage_id)
        self.phage_name = phage.get_name()
        self.name = self.phage_name + "_" + self.gene_no
        self.cluster = phage.cluster
        if self.cluster is None:
            self.cluster = "singleton"
        if phage.subcluster:
            self.subcluster = phage.subcluster
        else:
            self.subcluster = self.cluster
        self.cluster_hash = sum([pow(ord(elem), i+1) for i, elem in enumerate(self.subcluster)])
        status = phage.get_status()
        if status == 'final':        # values of 'draft' or 'gbk' considered draft quality by starterator
            self.draftStatus = False
        else:
            self.draftStatus = True

        phage_sequence = phage.get_sequence()
        if self.orientation == 'R':
            temp_start = self.stop
            self.stop = self.start
            self.start = temp_start
        self.genome_length = len(phage_sequence)
        sequence, self.ahead_of_start = find_upstream_stop_site(
                                self.start, self.stop, self.orientation, phage_sequence)
        self.ahead_of_start_coord = self.start - self.ahead_of_start
        gene = SeqRecord(sequence, id=self.gene_id, name=self.name,
                         description="|%i-%i| %s" % (self.start, self.stop, self.orientation))
        return gene

    def add_candidate_starts(self):
        """
            Finds all the possible start site of the gene and returns a list of indexes of start sites
        """
        gene_sequence = self._get_gene_seq_str()
        starts = []
        start_codons = {'ATG', 'GTG', 'TTG'}
        for index in range(0, len(gene_sequence), 3):
            codon = gene_sequence[index:index+3]
            if codon in start_codons:
                starts.append(index)
        return sorted(starts)

    def _get_gene_seq_str(self):
        seq = self.sequence.seq
        seq_key = (id(seq), len(seq))
        cached = getattr(self, "_gene_seq_str", None)
        if cached is None or getattr(self, "_gene_seq_str_key", None) != seq_key:
            cached = str(seq)
            self._gene_seq_str = cached
            self._gene_seq_str_key = seq_key
        return cached

    def _get_alignment_seq_str(self):
        return self._get_alignment_analysis_cache()["sequence"]

    def _get_alignment_analysis_cache(self):
        seq = self.alignment.seq
        seq_len = len(seq)
        seq_key = (id(seq), seq_len)
        cache = getattr(self, "_alignment_analysis_cache", None)
        if cache is not None and cache.get("seq_key") == seq_key:
            return cache

        alignment_sequence = str(seq)
        # Convert the alignment string to a byte array for vectorized ops.
        # Each character becomes its ASCII code (e.g. 'A'->65, '-'->45).
        arr = np.frombuffer(alignment_sequence.encode('ascii'), dtype=np.uint8)

        # Boolean mask: True at every position that is NOT a gap character.
        is_not_gap = arr != ord('-')

        # non_gap_prefix[i] = number of non-gap characters in alignment[:i].
        # This is a running total so that prefix[end] - prefix[start] gives
        # the count of real bases in any alignment slice.
        non_gap_prefix = np.empty(seq_len + 1, dtype=np.int32)
        non_gap_prefix[0] = 0
        if seq_len > 0:
            np.cumsum(is_not_gap, out=non_gap_prefix[1:])

        # Alignment indices where non-gap characters appear.
        # e.g. for "A--TG" this would be [0, 3, 4].
        non_gap_to_alignment_index = np.flatnonzero(is_not_gap)

        # Boolean mask: True at positions containing one of A, G, T, C.
        is_acgt = ((arr == ord('A')) | (arr == ord('G'))
                   | (arr == ord('T')) | (arr == ord('C')))

        # Alignment indices where ACGT characters appear (excludes gaps
        # and any ambiguous/lowercase bases).
        agtc_to_alignment_index = np.flatnonzero(is_acgt)

        # Run-length encoding of contiguous seq/gap segments.
        # Each run is (start_index, end_index, 'seq'|'gap') where 'seq'
        # means the segment contains ACGT bases and 'gap' means it doesn't.
        if seq_len > 0:
            # Cast the boolean is_acgt to 0/1 so np.diff detects transitions.
            segment_types = is_acgt.view(np.uint8)
            # Indices where the segment type changes (0->1 or 1->0).
            boundaries = np.flatnonzero(np.diff(segment_types)) + 1
            # Build parallel start/end arrays from the boundary positions.
            starts = np.empty(len(boundaries) + 1, dtype=np.intp)
            starts[0] = 0
            starts[1:] = boundaries
            ends = np.empty(len(boundaries) + 1, dtype=np.intp)
            ends[:-1] = boundaries
            ends[-1] = seq_len
            # Look up the type (1=seq, 0=gap) at each run's first position.
            types = segment_types[starts]
            base_feature_runs = tuple(
                (int(s), int(e), 'seq' if t else 'gap')
                for s, e, t in zip(starts, ends, types)
            )
        else:
            base_feature_runs = ()
        feature_runs_by_boundary = {0: base_feature_runs, seq_len: base_feature_runs}
        coord_cache = {}

        cache = {
            "seq_key": seq_key,
            "seq_len": seq_len,
            "sequence": alignment_sequence,
            "non_gap_prefix": non_gap_prefix.tolist(),
            "non_gap_to_alignment_index": non_gap_to_alignment_index.tolist(),
            "agtc_to_alignment_index": agtc_to_alignment_index.tolist(),
            "base_feature_runs": base_feature_runs,
            "feature_runs_by_boundary": feature_runs_by_boundary,
            "coord_cache": coord_cache,
        }
        self._alignment_analysis_cache = cache

        # Keep legacy cache attributes in sync for compatibility with existing call sites.
        self._alignment_seq_str = alignment_sequence
        self._alignment_seq_str_key = seq_key
        self._alignment_non_gap_prefix = cache["non_gap_prefix"]
        self._alignment_non_gap_prefix_key = seq_key
        self._alignment_index_coord_cache = coord_cache
        self._alignment_index_coord_cache_key = seq_key
        return cache

    def _compute_alignment_start_and_candidates(self):
        alignment_cache = self._get_alignment_analysis_cache()
        alignment_seq_key = alignment_cache["seq_key"]

        candidate_lookup_key = (id(self.candidate_starts), len(self.candidate_starts))
        if getattr(self, "_candidate_starts_lookup_key", None) != candidate_lookup_key:
            unique_candidates = tuple(sorted(set(self.candidate_starts)))
            self._candidate_starts_lookup_key = candidate_lookup_key
            self._candidate_starts_tuple = unique_candidates
        candidate_starts_tuple = self._candidate_starts_tuple

        scan_key = (alignment_seq_key, self.ahead_of_start, candidate_starts_tuple)
        if getattr(self, "_alignment_start_candidate_scan_key", None) == scan_key:
            start_site = self._alignment_start_site_cached
            aligned_starts = list(self._alignment_candidate_starts_cached)
            return start_site, aligned_starts

        agtc_positions = alignment_cache["agtc_to_alignment_index"]
        ahead_of_start = self.ahead_of_start

        if not agtc_positions:
            start_count = -1
            start_site = 0
        else:
            if ahead_of_start <= 0:
                start_count = 0
                start_site = agtc_positions[0]
            elif ahead_of_start < len(agtc_positions):
                start_count = ahead_of_start
                start_site = agtc_positions[ahead_of_start]
            else:
                start_count = len(agtc_positions) - 1
                start_site = agtc_positions[-1]

        if start_count > ahead_of_start:
            start_site -= 1

        non_gap_positions = alignment_cache["non_gap_to_alignment_index"]
        non_gap_length = len(non_gap_positions)
        aligned_starts = []
        for candidate_start in candidate_starts_tuple:
            if 0 <= candidate_start < non_gap_length:
                aligned_starts.append(non_gap_positions[candidate_start])

        aligned_starts_tuple = tuple(aligned_starts)
        self._alignment_start_candidate_scan_key = scan_key
        self._alignment_start_site_cached = start_site
        self._alignment_candidate_starts_cached = aligned_starts_tuple

        return start_site, list(aligned_starts_tuple)

    def _find_adjacent_start_groups(self):
        """Returns groups of start sites that are adjacent in the same ORF (exactly 3 apart)

        groups neighboring start sites, ignoring ones that aren't adjacent

        bad_starts should NOT be called. Starts where there is at least 1 start codon immediately following it.
        """

        starts_sorted = sorted(self.candidate_starts)
        if not starts_sorted:
            return []
        bad_groups = []
        current = [starts_sorted[0]]

        for s in starts_sorted[1:]:
            if s - current[-1] == 3:
                current.append(s)
            else:
                if len(current) >= 2:
                    bad_groups.append(current)
                current = [s]

        if len(current) >= 2:
            bad_groups.append(current)
        return bad_groups



    def add_alignment_start_site(self):
        """
            Gives the coordinate the called start site in the alignment sequence
        """
        start_site, aligned_starts = self._compute_alignment_start_and_candidates()
        self.alignment_start_site = start_site
        self.alignment_candidate_starts = aligned_starts
        return start_site

    def add_alignment_candidate_starts(self):
        """
            Creates a list of candidate starts of the alignment based on the candidate starts
            of the gene
        """
        start_site, aligned_starts = self._compute_alignment_start_and_candidates()
        self.alignment_start_site = start_site
        self.alignment_candidate_starts = aligned_starts
        return aligned_starts

    def add_alignment_start_stats(self, pham, lookup_cache=None):
        self.alignment_candidate_start_nums = []
        self.alignment_candidate_start_counts = []
        self.alignment_annot_start_nums = []
        self.alignment_annot_start_counts = []
        self.alignment_start_conservation = []
        self.alignment_start_num_called = None

        if self.pham_no is None:
            self.pham_no = pham.pham_no

        self.pham_size = len(pham.genes)

        num_gene_in_pham = len(pham.genes)
        use_lookup = bool(lookup_cache) and all(
            key in lookup_cache for key in (
                "candidate_start_nums_by_gene",
                "called_start_num_by_gene",
                "conservation_count_by_start",
                "annot_count_by_start",
            )
        )

        if use_lookup:
            candidate_start_nums_by_gene = lookup_cache["candidate_start_nums_by_gene"]
            called_start_num_by_gene = lookup_cache["called_start_num_by_gene"]
            conservation_count_by_start = lookup_cache["conservation_count_by_start"]
            annot_count_by_start = lookup_cache["annot_count_by_start"]

            self.alignment_candidate_start_nums = list(candidate_start_nums_by_gene.get(self.full_name, []))
            self.alignment_start_num_called = called_start_num_by_gene.get(self.full_name)

            for num in self.alignment_candidate_start_nums:
                conservation_count = conservation_count_by_start.get(num, 0)
                self.alignment_candidate_start_counts.append(conservation_count)

                if num_gene_in_pham > 0:
                    conserved_fraction = round(float(conservation_count) / float(num_gene_in_pham), 4)
                else:
                    conserved_fraction = 0.0
                self.alignment_start_conservation.append(conserved_fraction)

                annot_count = annot_count_by_start.get(num, 0)
                if annot_count > 0:
                    self.alignment_annot_start_nums.append(num)
                    self.alignment_annot_start_counts.append(annot_count)
        else:
            annotated = [gene.full_name for gene in pham.stats['most_common']['annot_list']]
            for startnum, genelist in pham.stats['most_common']['possible'].items():
                if self.full_name in genelist:
                    self.alignment_candidate_start_nums.append(startnum)

            for num in self.alignment_candidate_start_nums:
                conservation_count = len(pham.stats['most_common']['possible'][num])
                self.alignment_candidate_start_counts.append(conservation_count)

                conserved_fraction = round(float(conservation_count) / float(num_gene_in_pham), 4)
                self.alignment_start_conservation.append(conserved_fraction)

                annot_count = 0
                for gene in pham.stats['most_common']['called_starts'][num]:
                    if gene in annotated:
                        annot_count += 1

                if annot_count > 0:
                    self.alignment_annot_start_nums.append(num)
                    self.alignment_annot_start_counts.append(annot_count)

            for startnum, genelist in pham.stats['most_common']['called_starts'].items():
                if self.full_name in genelist:
                    self.alignment_start_num_called = startnum

        if len(self.alignment_annot_start_counts) > 0:
            most_annot_count = max(self.alignment_annot_start_counts)
            most_annot_index = self.alignment_annot_start_counts.index(most_annot_count)
            most_annotated_start_num = self.alignment_annot_start_nums[most_annot_index]
        else:
            most_annotated_start_num = None

        if most_annotated_start_num is None:
            self.calls_most_annotated = None
            self.has_most_annotated = None
        else:
            if self.alignment_start_num_called == most_annotated_start_num:
                self.calls_most_annotated = True
            else:
                self.calls_most_annotated = False

            if most_annotated_start_num in self.alignment_candidate_start_nums:
                self.has_most_annotated = True
            else:
                self.has_most_annotated = False

        total_annots = sum(self.alignment_annot_start_counts)
        self.alignment_annot_start_fraction = [float(count)/float(total_annots) for count in self.alignment_annot_start_counts]

        self.alignment_annot_counts_by_start = dict(zip(self.alignment_annot_start_nums, self.alignment_annot_start_counts))


        #Adjacent-start checking: determine whether the CALLED start is one of the bad ones
        # bad_adjacent_candidate_starts are OFFSETS (bp) in self.sequence coordinates (0-based into gene sequence)
        # We want to flag only if the called start corresponds to one of those "bad" offsets (NOT the last in a run).

        self.bad_adjacent_start_nums = []
        self.called_start_is_bad = False

        try:
            # offset(bp) -> alignment index
            if self.candidate_starts and self.alignment_candidate_starts:
                offset_to_aln = dict(zip(self.candidate_starts, self.alignment_candidate_starts))
            else:
                offset_to_aln = {}

            start_num_by_alignment_index = {}
            if use_lookup:
                start_num_by_alignment_index = lookup_cache.get("start_num_by_alignment_index", {})
            total_possible = getattr(pham, "total_possible_starts", None)

            if (start_num_by_alignment_index or total_possible) and offset_to_aln and getattr(self, "bad_adjacent_candidate_starts", None):
                bad_nums = set()

                for off in self.bad_adjacent_candidate_starts:
                    aln_idx = offset_to_aln.get(off)
                    if aln_idx is None:
                        continue
                    start_num = start_num_by_alignment_index.get(aln_idx)
                    if start_num is not None:
                        bad_nums.add(start_num)
                    elif total_possible and aln_idx in total_possible:
                        # start num is 1-based index into total_possible
                        bad_nums.add(total_possible.index(aln_idx) + 1)

                self.bad_adjacent_start_nums = sorted(bad_nums)

                # Only flag red if CALLED start num is one of the bad nums.
                self.called_start_is_bad = (self.alignment_start_num_called in bad_nums)

        except Exception:
            # keep defaults if anything goes wrong
            self.bad_adjacent_start_nums = []
            self.called_start_is_bad = False


        return





    def alignment_index_to_coord(self, index):
        """
                Given an index of the alignment
                finds the coordinates of the index on the phage sequence.
                The coordinate is 1 based count, not zero based
        """
        new_start_index = 0
        for i in range(0, index):
            if self.alignment.seq[i] != '-':
                new_start_index += 1
        if self.orientation == 'R':
            new_start_coords = (self.start + self.ahead_of_start - new_start_index - 1) % self.genome_length + 1
        else:
            new_start_coords = (self.start - self.ahead_of_start + new_start_index + 1)
        return new_start_coords

    def alignment_index_to_coord_optimized(self, index):
        seq_len, prefix, coord_cache = self._get_alignment_coord_lookup()
        index = self._normalize_alignment_index(index, seq_len)
        coord = coord_cache.get(index)
        if coord is None:
            coord = self._alignment_coord_from_prefix(prefix[index])
            coord_cache[index] = coord
        return coord

    def alignment_indices_to_coords_optimized(self, indices):
        seq_len, prefix, coord_cache = self._get_alignment_coord_lookup()
        coords = []
        for raw_index in indices:
            index = self._normalize_alignment_index(raw_index, seq_len)
            coord = coord_cache.get(index)
            if coord is None:
                coord = self._alignment_coord_from_prefix(prefix[index])
                coord_cache[index] = coord
            coords.append(coord)
        return coords

    def _normalize_alignment_index(self, index, seq_len):
        index = int(index)
        if index < 0:
            return 0
        if index > seq_len:
            return seq_len
        return index

    def _alignment_coord_from_prefix(self, non_gap_count):
        if self.orientation == 'R':
            return (self.start + self.ahead_of_start - non_gap_count - 1) % self.genome_length + 1
        return self.start - self.ahead_of_start + non_gap_count + 1

    def _get_alignment_coord_lookup(self):
        alignment_cache = self._get_alignment_analysis_cache()
        return (
            alignment_cache["seq_len"],
            alignment_cache["non_gap_prefix"],
            alignment_cache["coord_cache"],
        )

    def _get_feature_runs_for_boundary(self, start_boundary):
        alignment_cache = self._get_alignment_analysis_cache()
        feature_runs_by_boundary = alignment_cache["feature_runs_by_boundary"]
        cached_runs = feature_runs_by_boundary.get(start_boundary)
        if cached_runs is not None:
            return cached_runs

        base_feature_runs = alignment_cache["base_feature_runs"]
        if not base_feature_runs:
            feature_runs = ()
        else:
            split_runs = []
            for segment_start, segment_end, segment_type in base_feature_runs:
                if segment_start < start_boundary < segment_end:
                    split_runs.append((segment_start, start_boundary, segment_type))
                    split_runs.append((start_boundary, segment_end, segment_type))
                else:
                    split_runs.append((segment_start, segment_end, segment_type))
            feature_runs = tuple(split_runs)

        feature_runs_by_boundary[start_boundary] = feature_runs
        return feature_runs

    def add_gaps_as_features(self, feature_template_cache=None):
        alignment_cache = self._get_alignment_analysis_cache()
        sequence = alignment_cache["sequence"]
        sequence_len = alignment_cache["seq_len"]

        start_boundary = self.alignment_start_site
        if start_boundary < 0:
            start_boundary = 0
        elif start_boundary > sequence_len:
            start_boundary = sequence_len

        cache_key = (sequence, start_boundary)
        if feature_template_cache is not None:
            cached_runs = feature_template_cache.get(cache_key)
            if cached_runs is not None:
                self.alignment_feature_runs = cached_runs
                self.alignment.features = []
                return

        if sequence_len == 0:
            self.alignment_feature_runs = ()
            self.alignment.features = []
            if feature_template_cache is not None:
                feature_template_cache[cache_key] = ()
            return

        feature_runs = self._get_feature_runs_for_boundary(start_boundary)
        if feature_template_cache is not None:
            feature_template_cache[cache_key] = feature_runs
        self.alignment_feature_runs = feature_runs
        self.alignment.features = []

    def has_valid_start(self):
        return self.ahead_of_start in self.candidate_starts

    def is_equal(self, other):
        """
            Checks if another PhamGene is equal to this one
            PhamGenes are equal if they have the same called start, if the amount ahead of start
            (amount of sequence before the previous stop site) is the same, if the candidate
            starts of the genes are the same, and if the alignment gaps or not are the same
            (This is essentially, they would look the same on the graph output)

        """
        return self._comparison_signature() == other._comparison_signature()

    def _feature_signature(self):
        if not self.sequence.features:
            return ()
        features = []
        for feature in self.sequence.features:
            features.append((int(feature.location.start), int(feature.location.end), feature.type))
        features.sort()
        return tuple(features)

    def _comparison_signature(self):
        return (
            self.alignment_start_site,
            self.ahead_of_start,
            tuple(self.alignment_candidate_starts),
            self._feature_signature(),
        )

    def get_locustag(self):
        db_return = get_db().get(
                "SELECT phage.annotationauthor, phage.status, gene.locustag from gene JOIN phage on gene.phageid=phage.phageid where gene.geneid = %s",
                self.db_id)

        self.annot_author = db_return[0]
        self.status = db_return[1]
        self.locustag = db_return[2]
        return

    def __repr__(self):
        return 'Phamgene for %s' % self.gene_id


class UnPhamGene(PhamGene):
    def __init__(self, number, start, stop, orientation, phage_name, phage_sequence):
        self.number = number
        self.phage_name = phage_name
        self.gene_id = "%s_%s" % (phage_name, number)
        self.full_name = self.gene_id
        self.name = number
        self.start = start-1
        self.stop = stop
        self.orientation = orientation
        self.pham_size = None
        self.pham_no = None
        self.cluster = "Unassigned"
        self.subcluster = "Unassigned"
        self.cluster_hits = None
        self.subcluster_hits = None
        self.cluster_hash = sum([pow(ord(elem), i + 1) for i, elem in enumerate(self.subcluster)])


        if orientation == 'R':
            self.start_codon_location = stop
            self.stop_codon_location = start
        else:
            self.start_codon_location = start
            self.stop_codon_location = stop

        self.sequence = self.make_gene(phage_sequence)
        self.candidate_starts = self.add_candidate_starts()

        #find all but the last start
        self.adjacent_candidate_start_groups = self._find_adjacent_start_groups()

        bad = []
        for grp in self.adjacent_candidate_start_groups:
            if len(grp) >= 2:
                bad.extend(grp[:-1])
        self.bad_adjacent_candidate_starts = sorted(set(bad))
        self.has_bad_adjacent_candidate_starts = bool(self.bad_adjacent_candidate_starts)





        '''
        if self.has_bad_adjacent_candidate_starts:
            print(
                f"bad adjacent starts (offsets) gene={getattr(self, 'gene_no', getattr(self, 'number', '?'))}: {self.bad_adjacent_candidate_starts}")
        '''


        self.alignment = None
        self.alignment_start = None
        self.alignment_candidate_starts = None
        self.alignment_feature_runs = None
        self.alignment_candidate_start_nums = None
        self.alignment_annot_start_nums = None
        self.alignment_annot_start_counts = None
        self.alignment_start_num_called = None
        self.calls_most_annotated = None
        self.has_most_annotated = None
        self.suggested_start = {}
        self.draftStatus = True

    def make_gene(self, phage_sequence):
        if self.orientation == 'R':
            temp_start = self.stop
            self.stop = self.start
            self.start = temp_start
        sequence, self.ahead_of_start = find_upstream_stop_site(
                                self.start, self.stop, self.orientation, phage_sequence)
        gene = SeqRecord(sequence, id=self.gene_id, name=self.gene_id)
        return gene

    def phambymatch(self):
        #try to must make a perfect match to the translation field
        protein = str(self.sequence[self.ahead_of_start:].seq.translate())
        #repair translations if start codon was TTG or GTG and remove stop codon
        protein = "M" + protein[1:-1]
        db = DB()
        result = db.query("SELECT GeneID FROM gene WHERE gene.Translation = %s", protein)
        if len(result) < 1:
            return None
        else:
            result2 = db.query("SELECT phamid FROM gene WHERE geneid = %s", result[0])
            print("pham %s by exact match to gene %s"%(result2[0],result[0]))
            number, = result2[0]
            self.pham_no = number

            return self.pham_no

    def blast(self):
        """
        Runs BLASTp for this UnPhamGene (if needed) and returns the pham number.
        Uses cached XML if present, but deletes it if it's empty (common point of failure).
        """
        xml_path = os.path.join(utils.INTERMEDIATE_DIR, f"{self.gene_id}.xml")

        # If cached XML exists but is empty, delete it so BLAST reruns
        if os.path.exists(xml_path) and os.path.getsize(xml_path) == 0:
            os.remove(xml_path)

        # If XML doesn't exist, run BLAST
        if not os.path.exists(xml_path):
            protein = SeqRecord(self.sequence[self.candidate_starts[0]:].seq.translate(), id=self.gene_id)

            # short proteins need lower e_value
            query_len = (self.stop - self.start) / 3
            if query_len < 50:
                e_value = math.pow(10, -5)
            else:
                e_value = math.pow(10, -20)

            fasta_path = os.path.join(utils.INTERMEDIATE_DIR, f"{self.gene_id}.fasta")
            SeqIO.write(protein, fasta_path, "fasta")

            db_path = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")  # no quotes
            # Ensure the BLAST database exists (makeblastdb outputs)
            db_required = [db_path + ext for ext in (".pin", ".psq", ".phr")]
            if not all(os.path.exists(p) for p in db_required):
                update_protein_db()

            blast_args = [
                os.path.join(utils.BLAST_DIR, "blastp"),
                "-out", xml_path,
                "-outfmt", "5",
                "-query", fasta_path,
                "-db", db_path,
                "-evalue", str(e_value),
            ]

            try:
                subprocess.check_call(blast_args)
            except Exception as e:
                raise StarteratorError(f"Blast could not run! ({e})")

            # If BLAST produced an empty XML anyway, fail clearly
            if not os.path.exists(xml_path) or os.path.getsize(xml_path) == 0:
                raise StarteratorError("BLAST produced an empty XML output file.")

        return self.parse_blast()

    def parse_blast(self):
        # Always read the BLAST XML from the intermediate directory.
        xml_path = os.path.join(utils.INTERMEDIATE_DIR, f"{self.gene_id}.xml")

        # First try: NCBIXML.read (expects exactly one record)
        try:
            with open(xml_path, "r") as result_handle:
                blast_record = NCBIXML.read(result_handle)
        except:
            # Fallback: NCBIXML.parse (iterator) in case file contains multiple records
            with open(xml_path, "r") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                blast_record = next(blast_records)

        if len(blast_record.descriptions) > 0:
            first_result = blast_record.descriptions[0].title.split(',')[0].split(' ')[-1]
            # print first_result
            if "_" not in first_result:
                first_result = blast_record.descriptions[1].title.split(',')[0].split(' ')[-1]
            # Try to get pham directly from gene name
            db = DB()
            results = db.query("SELECT phamid from gene where geneID = %s", first_result)
            if len(results) == 1:
                number, = results[0]
                self.pham_no = number
                return number
            else:
                phage_name = first_result.split("_")[0]
                #exception for error in locus tags specific to this phage:
                if phage_name == 'FRIAPREACHER':
                    phage_name = 'FRIARPREACHER'
                if phage_name.lower() == "draft":
                    phage_name = first_result.split("_")[-3]
                gene_number = first_result.split("_")[-1]
                # print phage_name, gene_number
                pham_no = get_pham_no(phage_name, gene_number)
                self.pham_no = pham_no
                return pham_no
        else:
            self.pham_no = None
            return None

    def add_cluster_hits(self):
        self.cluster_hits = []
        db = DB()
        result = db.query("SELECT distinct(phage.cluster) FROM phage JOIN gene on gene.phageid = phage.phageid\
                           WHERE gene.phamid = %s", self.pham_no)

        for item in result:
            hit = list(item)[0]
            if hit is not None:
                self.cluster_hits.append(hit)

        self.subcluster_hits = []
        result2 = db.query("SELECT distinct(phage.subcluster) FROM phage JOIN gene on gene.phageid = phage.phageid\
                           WHERE gene.phamid = %s", self.pham_no)
        for item in result2:
            hit2 = list(item)[0]
            if hit2 is not None:
                self.subcluster_hits.append(hit2 )

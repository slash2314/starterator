import gzip
import logging
import math
import time

from starterator.database import DB, get_db
from starterator.phamgene import new_PhamGene
from starterator.phage import new_phage
from Bio import AlignIO
from Bio import SeqIO
from collections import Counter
from . import utils
import subprocess
import multiprocessing
import os
import errno
import sys
from .utils import StarteratorError
import json
import orjson


def get_worker_count():
    """Return the number of parallel workers for pham processing.

    Uses STARTERATOR_WORKERS env var if set, otherwise defaults to half of
    cpu_count() to avoid hyperthreading contention — most CPUs report
    logical cores (2 per physical core), and CPU-bound workloads gain
    little from sibling hyperthreads.
    """
    env = os.environ.get("STARTERATOR_WORKERS")
    if env is not None:
        return max(1, int(env))
    return max(1, multiprocessing.cpu_count() // 2)
from hashlib import sha256



def generate_pham_hashes(output_file=True):
    db = get_db()

    db_version = db.query("SELECT Version from version;")[0][0]

    pham_ids = [row[0] for row in db.query("SELECT PhamID from pham ORDER BY PhamID;")]
    total_phams = len(pham_ids)
    pham_hashes = {}
    BATCH_SIZE = 500
    for i in range(0, len(pham_ids), BATCH_SIZE):
        batch_phams = pham_ids[i:i+BATCH_SIZE]
        print(f"Processing hash batch {math.floor(i / BATCH_SIZE)+1}/{math.ceil(total_phams/BATCH_SIZE)}")
        query = """
        SELECT phage.PhageID, phage.Cluster, phage.Subcluster, phage.Status, pham.PhamID, gene.Translation, gene.Name, gene.LocusTag, gene.Notes, gene.GeneID
        FROM pham
        JOIN gene ON pham.PhamID = gene.PhamID
        JOIN phage ON gene.PhageID = phage.PhageID
        WHERE pham.PhamID IN %s
        ORDER BY pham.PhamID, gene.GeneID
        """
        db_results = db.query(query, (tuple(batch_phams),))
        pham_groups = {}
        for row in db_results:
            if row[4] not in pham_groups:
                pham_groups[row[4]] = []
            pham_groups[row[4]].append(row)
        # Generate hashes for each pham in the batch
        for pham_id, genes in pham_groups.items():
            sorted_genes = sorted(genes, key=lambda x: x[-1])  # Sort by GeneID
            membership_data = '|'.join([g[-1] for g in sorted_genes])
            sequences_data = '|'.join([utils.decode_if_bytes(g[5]) for g in sorted_genes])
            annotations_data = '|'.join([g[6] for g in sorted_genes])
            # LocusTag, Notes, Cluster, Subcluster, Status
            metadata_data = '|'.join([f"{g[7]}|{g[8]}|{g[1]}|{g[2]}|{g[3]}" for g in sorted_genes])
            membership_hash = sha256(membership_data.encode()).hexdigest()
            sequences_hash = sha256(sequences_data.encode()).hexdigest()
            annotations_hash = sha256(annotations_data.encode()).hexdigest()
            metadata_hash = sha256(metadata_data.encode()).hexdigest()
            combined_hash = sha256(f"{membership_hash}|{sequences_hash}|{annotations_hash}|{metadata_hash}".encode()).hexdigest()
            pham_hashes[pham_id] = {
                "membership_hash": membership_hash,
                "sequences_hash": sequences_hash,
                "annotations_hash": annotations_hash,
                "metadata_hash": metadata_hash,
                "combined_hash": combined_hash,
                "gene_count": len(sorted_genes)
            }
        pham_groups.clear()
    sorted_pham_ids = sorted(pham_hashes, key=lambda pham_id: str(pham_id))
    all_combined_hashes = '|'.join(
        [pham_hashes[pham_id]["combined_hash"] for pham_id in sorted_pham_ids]
    )
    overall_hash = sha256(all_combined_hashes.encode()).hexdigest()
    hash_output = {
        "database_version": db_version,
        # Iso time
        "generated_timestamp": time.ctime(),
        "overall_hash": overall_hash,
        "total_phams": len(pham_hashes),
        "phams": pham_hashes
    }
    if output_file:
        with open(os.path.join(utils.FINAL_DIR, f"pham_hashes_v{db_version}.json"), 'wb') as f:
            f.write(orjson.dumps(hash_output, option=orjson.OPT_INDENT_2 | orjson.OPT_NON_STR_KEYS))
        print(f"Output written to {utils.FINAL_DIR}/pham_hashes_v{db_version}.json")
    return hash_output

def compare_hash_snapshots(hash1, hash2):
    """Compare two in-memory pham hash snapshots."""
    hash_report = {
        "overall_hash_matches": hash1["overall_hash"] == hash2["overall_hash"],
    }
    phams1 = {str(pham_id): pham_info for pham_id, pham_info in hash1["phams"].items()}
    phams2 = {str(pham_id): pham_info for pham_id, pham_info in hash2["phams"].items()}
    modified_phams = []
    phams1_ids = set(phams1.keys())
    phams2_ids = set(phams2.keys())
    for pham_id, pham_info in phams1.items():
        if pham_id in phams2 and pham_info != phams2[pham_id]:
            modified_phams.append(pham_id)
    added_phams = phams2_ids - phams1_ids
    removed_phams = phams1_ids - phams2_ids
    hash_report["phams_added"] = list(added_phams)
    hash_report["phams_removed"] = list(removed_phams)
    hash_report["phams_modified"] = modified_phams
    return hash_report

def compare_hashes_current(hash_file_name1, hash_file_name2):
    """
    Compares two hash file outputs generated from generate_pham_hashes.
    Outputs a report with the phams with differences: add, removed, modified.
    """
    with open(hash_file_name1, 'r') as f:
        hash1 = json.load(f)
    with open(hash_file_name2, 'r') as f:
        hash2 = json.load(f)
    return compare_hash_snapshots(hash1, hash2)


def _validated_genes_for_pham(pham_no):
    """
    Build validated gene membership for a pham using Starterator's
    existing filtering behavior.
    """
    results = get_db().query(
        "SELECT `gene`.`GeneID`, `gene`.`phageID`, "
        " `gene`.`Length`, `gene`.`Start`, `gene`.`Stop`, `gene`.`Orientation`, `gene`.`name`"
        " FROM `gene`"
        " JOIN `phage` ON `gene`.`PhageID` = `phage`.`PhageID`"
        " WHERE `gene`.`PhamID` = %s"
        " AND `phage`.`Status` != 'unknown'; ", pham_no)

    # --- Batch phage prefetch ---
    # Collect unique phage IDs from results, filter to those not already cached,
    # and load them all in one query to avoid N+1 individual get_name() calls.
    from . import phage as _phage_mod
    unique_phage_ids = {row[1] for row in results}
    uncached_ids = [pid for pid in unique_phage_ids if pid not in _phage_mod.phage_list]
    if uncached_ids:
        phage_rows = get_db().query(
            "SELECT PhageID, Name, Cluster, Sequence, Status, AnnotationAuthor, Subcluster"
            " FROM phage WHERE PhageID IN %s",
            (tuple(uncached_ids),)
        )
        from .phage import Phage
        for row in phage_rows:
            phage_id = row[0]
            p = Phage(phage_id=phage_id)
            p.name = row[1]
            p.cluster = row[2]
            p.sequence = row[3]
            p.status = row[4]
            p.annot_author = row[5]
            p.subcluster = row[6]
            _phage_mod.phage_list[phage_id] = p

    genes = {}
    phage_has_n = {}
    skipped_n = []
    skipped_start = []

    for gene_info in results:
        gene_id = gene_info[0]
        phage_id = gene_info[1]
        start = gene_info[3]
        stop = gene_info[4]
        orientation = gene_info[5]
        name = gene_info[6]
        gene = new_PhamGene(gene_id, start, stop, orientation, phage_id, name)

        if phage_id not in phage_has_n:
            phage = new_phage(phage_id=phage_id)
            genome_seq = utils.decode_if_bytes(phage.get_sequence())
            phage_has_n[phage_id] = "N" in genome_seq
        if phage_has_n[phage_id]:
            skipped_n.append(gene_id)
            continue

        if gene.has_valid_start():
            genes[gene.gene_id] = gene
        else:
            skipped_start.append(gene_id)

    return {
        "genes": genes,
        "total_count": len(results),
        "skipped_n": skipped_n,
        "skipped_start": skipped_start,
    }


def get_all_phams():
    """Get relevant pham numbers from the database.

    Uses a single SQL query to pre-filter likely candidate phams, then
    applies Starterator's per-gene validation so only phams with at least
    two retained genes are emitted.
    """
    db = get_db()
    candidate_rows = db.query(
        "SELECT gene.PhamID, COUNT(*) as cnt"
        " FROM gene JOIN phage ON gene.PhageID = phage.PhageID"
        " WHERE phage.Status != 'unknown' AND phage.Sequence NOT LIKE '%N%'"
        " GROUP BY gene.PhamID HAVING cnt >= 2"
        " ORDER BY cnt DESC"
    )
    relevant_phams = []
    for pham_no, _count in candidate_rows:
        validated = _validated_genes_for_pham(pham_no)
        if len(validated["genes"]) >= 2:
            relevant_phams.append(pham_no)
    return relevant_phams

def _process_pham_inprocess(args):
    """In-process worker: call PhamReport.final_report() directly, no subprocess."""
    pham_no, no_pdfs, compress_json = args
    start_time = time.time()
    logging.info(f"Started processing Pham {pham_no}")
    # Spawned worker processes start with a fresh module state where
    # utils.INTERMEDIATE_DIR / FINAL_DIR are empty strings (module-level defaults).
    # Call get_config() to initialize them from the config file / env vars,
    # matching what main() does in the subprocess path.
    utils.get_config()
    # Clear module-level caches to prevent cross-pham contamination within this worker.
    # These caches (phage_list, pham_genes) are singletons that would carry stale data
    # from the previous pham processed by this same worker process.
    from . import phage as _phage_mod
    from . import phamgene as _phamgene_mod
    _phage_mod.phage_list.clear()
    _phamgene_mod.pham_genes.clear()
    try:
        from .report import PhamReport
        report = PhamReport(pham_no)
        report.final_report(save_json=True, no_pdfs=no_pdfs, compress_json=compress_json)
    except Exception as e:
        # Match subprocess behavior: log the error and continue to the next pham
        # instead of killing the entire batch.
        logging.error(f"Pham {pham_no} failed: {e}")
        return
    elapsed = time.time() - start_time
    if elapsed < 1:
        logging.info(f"Finished processing Pham {pham_no} in {elapsed*1000:.0f}ms")
    else:
        logging.info(f"Finished processing Pham {pham_no} in {elapsed:.1f}s")


def start_pham_job(pham_no, no_pdfs=False, compress_json=False):
    """
    Start a subprocess for the given pham number.
    Kept for single-pham CLI invocation and harness profiling.
    """
    start_time = time.time()
    logging.info(f"Started processing Pham {pham_no}")
    command = [sys.executable, '-m', 'starterator.starterate', '-n', str(pham_no), '-j', 'True']
    if no_pdfs:
        command.append('--no-pdfs')
    if compress_json:
        command.append('--compress-json')
    subprocess.call(command)
    elapsed = time.time() - start_time
    if elapsed < 1:
        logging.info(f"Finished processing Pham {pham_no} in {elapsed*1000:.0f}ms")
    else:
        logging.info(f"Finished processing Pham {pham_no} in {elapsed:.1f}s")

def process_pham_list(phams, no_pdfs=False, compress_json=False):
    """Process an explicit iterable of pham IDs using in-process workers."""
    phams = list(phams)
    workers = get_worker_count()
    logging.info(f"Found {len(phams)} phams to process using {workers} workers")
    jobs = [(pham_no, no_pdfs, compress_json) for pham_no in phams]
    with multiprocessing.Pool(processes=workers) as pool:
        list(pool.imap_unordered(_process_pham_inprocess, jobs, chunksize=1))
    logging.info("Batch pham processing complete!")

def process_all_phams(no_pdfs=False, compress_json=False):
    """
    Process all phams in the database using in-process workers.
    """
    process_pham_list(get_all_phams(), no_pdfs=no_pdfs, compress_json=compress_json)

def process_phams_from_file(filepath, no_pdfs=False, compress_json=False):
    """Process phams listed in a line-delimited file using in-process workers."""
    with open(filepath, 'r') as f:
        phams = [line.strip() for line in f if line.strip()]
    logging.info(f"Loaded {len(phams)} phams from {filepath}")
    process_pham_list(phams, no_pdfs=no_pdfs, compress_json=compress_json)

def get_pham_number(phage_name, gene_number):
    try:
        db = DB()
        results = db.query("SELECT pham.Name \n\
            FROM gene JOIN pham ON gene.GeneID = pham.Gene \n\
            JOIN phage ON gene.PhageID = phage.PhageID \n\
            WHERE phage.Name LIKE %s AND gene.Name LIKE %s \n\
            ESCAPE '!'", (phage_name+"%", '%'+str(gene_number)))
        row = results[0]
        pham_no = row[0]
        return str(pham_no)
    except:
        raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))


def get_pham_colors(phams=None):
    db = DB()
    results = db.query("SELECT `PhamID`, `Color` from `pham`;")
    pham_colors = {}
    if phams:
        for row in results:
            if row[0] in phams:
                pham_colors[str(row[0])] = row[1]
    else:
        for row in results:
            pham_colors[str(row[0])] = row[1]
    return pham_colors

def get_version():
    results = get_db().query("SELECT Version from version;")
    return int(results[0][0])


def _seq_runs_from_feature_runs(feature_runs):
    """Convert internal feature runs to compact JSON-friendly sequence spans."""
    if not feature_runs:
        return []
    return [[start, end] for start, end, feature_type in feature_runs if feature_type == "seq"]


class Pham(object):
    def __init__(self, pham_no, genes=None):
        self.pham_no = pham_no
        self.stats = {}
        self.file = ""
        self.count = 0
        self.genes = self.get_genes()
        self.color = self.get_color()
        if genes:
            for gene in genes:
                self.add(gene)
            whole = "All" if len(genes) > 1 else "One"
            self.file = "%s%s" % (genes[0].phage_name, whole)
        self.aligner = None

    def get_genes(self):
        """
            Get the genes of the Phamily
        """
        validated = _validated_genes_for_pham(self.pham_no)
        genes = validated["genes"]
        skipped_n = validated["skipped_n"]
        skipped_start = validated["skipped_start"]
        self.count = validated["total_count"]
        if len(genes) < 1:
            if self.count == 0:
                raise StarteratorError("Pham %s not found in database." % self.pham_no)
            reasons = []
            if skipped_n:
                reasons.append("%d gene(s) skipped due to N's in genome: %s" % (len(skipped_n), ", ".join(skipped_n)))
            if skipped_start:
                reasons.append("%d gene(s) skipped due to invalid start codon: %s" % (len(skipped_start), ", ".join(skipped_start)))
            raise StarteratorError(
                "Pham %s has %d gene(s) but all failed validation. %s" % (self.pham_no, self.count, "; ".join(reasons))
            )
        return genes

    def get_phage_genes(self):
        pass

    def add(self, gene):
        """
            Add an unphameratored gene to the pham
        """
        self.genes[gene.gene_id] = gene

    def get_color(self):
        """
            Get the color of the phamily from the database
        """
        try:
            result = get_db().get("SELECT `phamid`, `color`\n\
                FROM `pham` WHERE `phamid` = %s;", self.pham_no)
            return result[1]
        except:
            raise StarteratorError("Pham number %s not found in database!" % self.pham_no)

    def add_alignment(self, alignment):
        """
            Using the alignment, add the alignment to the each gene in the pham
        """
        feature_template_cache = {}
        for record in alignment:
            gene = self.genes[record.id]
            gene.alignment = record
            # Call the underlying computation once instead of through two wrappers.
            start_site, aligned_starts = gene._compute_alignment_start_and_candidates()
            gene.alignment_start_site = start_site
            gene.alignment_candidate_starts = aligned_starts
            gene.add_gaps_as_features(feature_template_cache=feature_template_cache)

    def call_clustal(self, fasta_file):
        # self.aligner = 'ClustalO'
        self.aligner = 'MAFFT'
        #self.aligner = 'ClustalW'

        if self.aligner == 'ClustalO':
            outfile = fasta_file.replace(".fasta", ".aln")
            subprocess.check_call(['clustalo', '--infile=%s' % fasta_file, '--outfile=%s' % outfile, '--outfmt=clu', '--threads', str(multiprocessing.cpu_count())])
            # subprocess.check_call(['clustalo', '--infile=%s' % fasta_file, '--outfile=%s' % outfile, '--outfmt=clu'])
        elif self.aligner == 'MAFFT':
            outfile = fasta_file.replace(".fasta", ".aln")
            mafft_fasta_out = outfile + ".fasta"
            # Keep MAFFT single-threaded because phams are processed in parallel workers.
            mafft_threads = 1
            with open(mafft_fasta_out, 'w') as alignment_file:
                try:
                    # subprocess.check_call(
                    #     ['mafft', '--retree', '2', '--quiet', '--thread', str(multiprocessing.cpu_count()), fasta_file],
                    #     stdout=alignment_file
                    # )
                    subprocess.check_call(
                        ['mafft', '--retree', '2', '--quiet', '--thread', str(mafft_threads), fasta_file],
                        stdout=alignment_file
                    )
                except OSError as e:
                    if e.errno == errno.ENOENT:
                        raise StarteratorError("MAFFT aligner selected but 'mafft' is not installed or not on PATH.")
                    raise StarteratorError("MAFFT alignment failed: %s" % e)
                except subprocess.CalledProcessError as e:
                    raise StarteratorError("MAFFT alignment failed with exit code %s" % e.returncode)

            try:
                alignment = AlignIO.read(mafft_fasta_out, "fasta")
                # Downstream start-site logic checks for uppercase A/G/T/C only,
                # so normalize MAFFT FASTA output before clustal conversion.
                for record in alignment:
                    record.seq = record.seq.upper()
                AlignIO.write(alignment, outfile, "clustal")
            except Exception as e:
                raise StarteratorError("Failed to convert MAFFT FASTA alignment to clustal: %s" % e)
            finally:
                if os.path.exists(mafft_fasta_out):
                    os.remove(mafft_fasta_out)
            # Return the in-memory alignment directly — skip re-reading the Clustal file.
            # The .aln file was still written above for caching by align().
            return alignment
        elif self.aligner == 'ClustalW':
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call(['clustalw', '-infile=%s' % (fasta_file), '-quicktree', '-quiet'], stderr=devnull, stdout=devnull)
        else:
            raise StarteratorError("Unknown aligner '%s'. Supported aligners are ClustalW, ClustalO, and MAFFT." % self.aligner)

        aln_file = fasta_file.replace(".fasta", ".aln")
        alignment = AlignIO.read(aln_file, "clustal")
        return alignment

    def make_fasta(self, file_name=None):
        if file_name is None:
            file_name = os.path.join(utils.INTERMEDIATE_DIR, "%sPham%s" % (self.file, self.pham_no))
        sorted_genes = sorted(self.genes.values(), key=lambda gene: gene.gene_id)
        genes = [gene.sequence for gene in sorted_genes]
        count = SeqIO.write(genes, "%s.fasta" % file_name, "fasta")

    def align(self):
        """
            Makes a fasta file of the genes in the Pham
            if the alignment already exists, uses that .aln as the alignment
            Otherwise, calls Clustalw from the command line and creates alignment
        """
        file_name = os.path.join(utils.INTERMEDIATE_DIR, "%sPham%s" % (self.file, self.pham_no))
        fasta_path = file_name + ".fasta"
        aln_path = file_name + ".aln"
        sorted_genes = sorted(self.genes.values(), key=lambda gene: gene.gene_id)
        genes = [gene.sequence for gene in sorted_genes]

        if len(self.genes) == 1:
            alignment = [genes[0]]
        else:
            try:
                alignment = AlignIO.read(aln_path, "clustal")
            except:
                SeqIO.write(genes, fasta_path, "fasta")
                alignment = self.call_clustal(fasta_path)
        self.add_alignment(alignment)

    def add_total_possible_starts(self):
        """ Returns a list of all the candidate starts from the alignment
        """
        unique_sites = set()
        for gene in self.genes.values():
            for site in gene.alignment_candidate_starts:
                unique_sites.add(site)
        self.total_possible_starts = sorted(unique_sites)
        return self.total_possible_starts

    def add_alignment_stats_to_phamgenes(self, lookup_cache=None):
        for gene in self.genes.values():
            gene.add_alignment_start_stats(self, lookup_cache=lookup_cache)
        return

    def group_similar_genes(self, start_with=None):
        """
            Groups genes that have the same called start site, the same candidate starts
            and the same alignment (gaps are the same) together
        start_with: phage to be first item of first list
        """
        genes = list(self.genes.values())
        groups = []
        groups_by_signature = {}

        for gene in genes:
            signature = gene._comparison_signature()
            group = groups_by_signature.get(signature)
            if group is None:
                group = []
                groups_by_signature[signature] = group
                groups.append(group)
            group.append(gene)

        if start_with:
            split_gene_lists_on = None
            start_with_name = start_with.lower()
            for i, gene_list in enumerate(groups):
                if len(gene_list) == 1:
                    if start_with_name == gene_list[0].phage_name.lower():
                        split_gene_lists_on = i
                else:
                    for j, gene in enumerate(gene_list):
                        if start_with_name == gene.phage_name.lower():
                            split_gene_lists_on = i
                            groups[i] = groups[i][j:] + groups[i][:j]
            if split_gene_lists_on is not None:
                groups = groups[split_gene_lists_on:] + groups[:split_gene_lists_on]

                if groups[0][0].subcluster != 'Unassigned':
                    sort_from = groups[0][0].cluster_hash
                else:
                    sort_from = 1

                if len(groups) > 1:
                    remaining = groups[1:]
                    remaining.sort(key=lambda x: abs(x[0].cluster_hash-sort_from))
                    groups[1:] = remaining

        else:
            groups.sort(key=lambda x: x[0].subcluster)

        return groups

    def find_most_common_start(self, ignore_draft=False):
        """
            From the total candidate strats of each gene in the pham and all the start
            called in each gene, finds the start that is most commonly called.
            Returns a dictionary containing:
                "most_called" : a list of genes currently call the "most common start"
                "most_not_called" : a list of genes that have the "most common start" but do not call it
                "no_most_called" : a list of genes that do no have the "most called start"
                "possible" : a list containing lists of genes with the start of the index of self.total_possible_starts
                "called_start: a list containing lists of genes with the called start of the index of self.total_possible_starts

            Also, for each gene in the pham, a suggested start is given, gene.suggested_start["most_commom"]
            For genes that have the most common start called (or not) a tuple containing the index of the
            most common start and the coordinate of the sequence is given.
            For genes that do not have the most common start, a list of all possible starts, containing the index
            (useful when looking at the graphical output), and the coordinate is given.
        """
        # TODO:
        # add functionality for ignoring DRAFT phages?
        # use term Called_start for all genes irrespective of method to determine location of start codon
        # use term Annotated_start for genes in which manual annotation was used to determine start codon
        # use term predicted_start for gene in which computational prediction was used to determine start codon
        genes = list(self.genes.values())
        all_start_sites = []
        all_annotated_start_sites = []
        all_predicted_start_sites = []
        possible_by_site = {}
        called_by_site = {}

        for gene in genes:
            start_site = gene.alignment_start_site
            gene_name = gene.full_name

            all_start_sites.append(start_site)
            if gene.draftStatus:
                all_predicted_start_sites.append(start_site)
            else:
                all_annotated_start_sites.append(start_site)

            if start_site in called_by_site:
                called_by_site[start_site].append(gene_name)
            else:
                called_by_site[start_site] = [gene_name]

            # Preserve per-site gene ordering while ensuring each gene contributes at most once per site.
            for site in dict.fromkeys(gene.alignment_candidate_starts):
                if site in possible_by_site:
                    possible_by_site[site].append(gene_name)
                else:
                    possible_by_site[site] = [gene_name]

        annotated_start_site_set = set(all_annotated_start_sites)
        start_stats = {}
        # creates two lists each containing a list of gene ids
        # for each candidate start of the pham:
        # start_stats["possible"] contains a list of genes with the candidate starts
        # for the index of each start in the pham
        # start_stats["called_starts"] contains of list of the genes that have the site
        #   of the index called as their start
        start_stats["possible"] = {}
        start_stats["called_starts"] = {}
        # start_stats["most_called"] = {}
        self.add_total_possible_starts()
        start_num_by_alignment_index = {}
        for i, site in enumerate(self.total_possible_starts):
            start_num = i + 1
            start_num_by_alignment_index[site] = start_num
            start_stats["possible"][start_num] = possible_by_site.get(site, [])
            # start_stats["most_called"][i+1] = []
            start_stats["called_starts"][start_num] = called_by_site.get(site, [])

        all_starts_count = Counter(all_start_sites)
        all_annot_count = Counter(all_annotated_start_sites)
        all_predicted_count = Counter(all_predicted_start_sites)

        called_starts_count = all_starts_count.most_common()
        annot_starts_count = all_annot_count.most_common()
        annot_count_by_alignment_index = dict(annot_starts_count)
        predicted_starts_count = all_predicted_count.most_common()

        most_called_start_index = None
        for start_site, _count in called_starts_count:
            if start_site in start_num_by_alignment_index:
                most_called_start_index = start_num_by_alignment_index[start_site]
                break
        if most_called_start_index is None:
            raise StarteratorError(
                "Unable to map called starts to candidate starts for pham %s." % self.pham_no
            )

        most_annot_start_index = None
        if len(annot_starts_count) > 0:  # i.e. at least 1 annotated gene
            for start_site, _count in annot_starts_count:
                if start_site in start_num_by_alignment_index:
                    most_annot_start_index = start_num_by_alignment_index[start_site]
                    break

        genes_start_most_called = start_stats["called_starts"][most_called_start_index]
        genes_start_most_called_set = set(genes_start_most_called)
        possible_most_called_set = set(start_stats["possible"][most_called_start_index])
        start_stats["most_called_start"] = most_called_start_index
        start_stats["most_annotated_start"] = most_annot_start_index

        if most_annot_start_index is not None:
            genes_start_most_annot_set = set(start_stats["called_starts"][most_annot_start_index])
            possible_most_annot_set = set(start_stats["possible"][most_annot_start_index])
        else:
            genes_start_most_annot_set = set()
            possible_most_annot_set = set()

        # start_stats["most_called"] = start_stats["called_starts"][most_called_start_index]
        start_stats["most_called"] = []
        start_stats["most_not_called"] = []
        start_stats["no_most_called"] = []
        start_stats["most_annotated"] = []
        start_stats["most_not_annotated"] = []
        start_stats["no_most_annot"] = []
        start_stats["annot_list"] = [g for g in genes if not g.draftStatus]
        start_stats["draft_list"] = [g for g in genes if g.draftStatus]
        annotated_set = {gene.full_name for gene in start_stats["annot_list"]}

        candidate_start_nums_by_gene = {}
        conservation_count_by_start = {}
        for start_num, gene_names in start_stats["possible"].items():
            conservation_count_by_start[start_num] = len(gene_names)
            for gene_name in gene_names:
                if gene_name in candidate_start_nums_by_gene:
                    candidate_start_nums_by_gene[gene_name].append(start_num)
                else:
                    candidate_start_nums_by_gene[gene_name] = [start_num]

        start_stats['called_counts'] = {}
        called_start_num_by_gene = {}
        annot_count_by_start = {}
        for start_num, gene_names in start_stats['called_starts'].items():
            if len(gene_names) > 0:
                start_stats['called_counts'][start_num] = len(gene_names)
            annot_count = 0
            for gene_name in gene_names:
                called_start_num_by_gene[gene_name] = start_num
                if gene_name in annotated_set:
                    annot_count += 1
            if annot_count > 0:
                annot_count_by_start[start_num] = annot_count

        start_stats['annot_counts'] = {}
        for gene in start_stats["annot_list"]:
            start_number = called_start_num_by_gene.get(gene.full_name)
            if start_number is None:
                continue
            if start_number not in start_stats['annot_counts']:
                start_stats['annot_counts'][start_number] = 0
            start_stats['annot_counts'][start_number] += 1
             #   start_stats['annot_counts'][k] = len(annotated)

        lookup_cache = {
            "candidate_start_nums_by_gene": candidate_start_nums_by_gene,
            "called_start_num_by_gene": called_start_num_by_gene,
            "conservation_count_by_start": conservation_count_by_start,
            "annot_count_by_start": annot_count_by_start,
            "start_num_by_alignment_index": start_num_by_alignment_index
        }

        most_called_alignment_index = self.total_possible_starts[most_called_start_index - 1]
        if most_annot_start_index is not None:
            most_annot_alignment_index = self.total_possible_starts[most_annot_start_index - 1]
        else:
            most_annot_alignment_index = None

        # print "phams.find_most_common_start: genes_start_most_called " + str(genes_start_most_called)
        for gene in genes:
            gene_name = gene.full_name
            # check if the gene even has the most called start
            if gene_name in possible_most_called_set:
                if gene_name in genes_start_most_called_set:
                    if gene.orientation == 'F':   # only +1 for forward genes
                        # genes where most called start is present and it is called as the start are "most_called"
                        gene.suggested_start["most_called"] = (most_called_start_index, gene.start+1)
                    else:
                        gene.suggested_start["most_called"] = (most_called_start_index, gene.start)
                    start_stats["most_called"].append(gene_name)
                else:
                    # genes where most called start is present but it's not the called start are "most_not_called
                    start_stats["most_not_called"].append(gene_name)
                    suggested_start = gene.alignment_index_to_coord_optimized(most_called_alignment_index)
                    gene.suggested_start["most_called"] = (most_called_start_index, suggested_start)

            else:
                # genes where the most called start is NOT even present are no_most_called
                start_stats["no_most_called"].append(gene_name)
                possible_start_nums = [
                    start_num_by_alignment_index[start] for start in gene.alignment_candidate_starts
                ]
                possible_start_coords = gene.alignment_indices_to_coords_optimized(gene.alignment_candidate_starts)
                possible_starts_coords = []
                for index, coord in zip(possible_start_nums, possible_start_coords):
                    possible_starts_coords.append((index, coord + 1))
                gene.suggested_start["most_called"] = possible_starts_coords

            if most_annot_start_index is not None:
                if gene_name in possible_most_annot_set:
                    if gene_name in genes_start_most_annot_set:
                        # code below used for deprecated "suggested starts" list
                        # if gene.orientation == 'F':  # only +1 for forward genes
                        #     # genes where most annotated start is present and it is called as the start are "most_annotated"
                        #     gene.suggested_start["most_annotated"] = (most_annot_start_index, gene.start + 1)
                        # else:
                        #     gene.suggested_start["most_annotated"] = (most_annot_start_index, gene.start)
                        start_stats["most_annotated"].append(gene_name)
                    else:
                        # genes where most annotated start is present but it's not the called start are "most_not_annotated"
                        start_stats["most_not_annotated"].append(gene_name)
                        # code below used for deprecated "suggested starts" list
                        suggested_start = gene.alignment_index_to_coord_optimized(
                            most_annot_alignment_index)  # +1 issue dealt with in function
                        gene.suggested_start["most_called"] = (most_annot_start_index, suggested_start)

                else:
                    # genes where the most annotated start is NOT even present are no_most_annot
                    start_stats["no_most_annot"].append(gene_name)
                    # Code below used for deprecated "suggested starts" list
                    # possible_starts_coords = []
                    # for start in gene.alignment_candidate_starts:
                    #     index = self.total_possible_starts.index(start) + 1
                    #     new_start = gene.alignment_index_to_coord(start) + 1
                    #     possible_starts_coords.append((index, new_start))
                    # gene.suggested_start["most_called"] = possible_starts_coords

            # section to add summary of annotations based only on set of starts found in the gene
            # start by looking through all possible starts in this particular gene and see if
            # there are any annotations that call that start

            alignment_start_coord_with_annotations = [
                start for start in gene.alignment_candidate_starts if start in annotated_start_site_set
            ]

            gene.suggested_start["alignment_start_coord_with_annotations"] = alignment_start_coord_with_annotations

            alignment_start_indices_with_annotations = [
                start_num_by_alignment_index[start] for start in alignment_start_coord_with_annotations
            ]

            gene.suggested_start["alignment_start_indices_with_annotations"] = alignment_start_indices_with_annotations

            alignment_start_counts_with_annotations = [
                annot_count_by_alignment_index[start]
                for start in alignment_start_coord_with_annotations
                if start in annot_count_by_alignment_index
            ]

            gene.suggested_start["alignment_start_counts_with_annotations"] = alignment_start_counts_with_annotations
            try:  # this happens when annotated start of gene is not one of the three typical start codons (ATG,GTG,TTG)
                gene.suggested_start["current_start_number"] = start_num_by_alignment_index[gene.alignment_start_site]
            except:
                gene.suggested_start["current_start_number"] = None

        self.stats["most_common"] = start_stats

        # now update genes based on start analysis
        self.add_alignment_stats_to_phamgenes(lookup_cache=lookup_cache)
        self.add_cluster_stats(start_stats)
        return start_stats

    def annot_summary(self):
        genes_missing_locustag = [gene for gene in self.genes.values() if gene.locustag is None and gene.db_id]
        if genes_missing_locustag:
            gene_ids = [gene.db_id for gene in genes_missing_locustag]
            placeholders = ",".join(["%s"] * len(gene_ids))
            results = get_db().query(
                "SELECT gene.geneid, phage.annotationauthor, phage.status, gene.locustag "
                "from gene JOIN phage on gene.phageid=phage.phageid "
                f"where gene.geneid IN ({placeholders})",
                tuple(gene_ids)
            )
            by_gene_id = {row[0]: row for row in results}
            for gene in genes_missing_locustag:
                row = by_gene_id.get(gene.db_id)
                if row is None:
                    continue
                gene.annot_author = row[1]
                gene.status = row[2]
                gene.locustag = row[3]

        summary_dict = {}
        summary_dict['Name'] = self.pham_no
        summary_dict['MemberCount'] = self.count
        summary_dict['AnnotCount'] = len(self.stats['most_common']['annot_list'])
        summary_dict['TotalStarts'] = len(self.total_possible_starts)
        summary_dict['DbVersion'] = get_version()
        summary_dict['Aligner'] = self.aligner
        summary_dict['TotalPossibleStartAlignmentIndices'] = list(self.total_possible_starts)

        genelist = []
        for gene in self.genes.values():
            gene_dict = {}
            gene_dict['GeneID'] = gene.gene_id
            gene_dict['Start'] = gene.start if gene.orientation == "R" else gene.start + 1 #switch to 1 based counting
            gene_dict['Stop'] = gene.stop  if gene.orientation == "F" else gene.stop + 1 #switch to 1 based counting
            gene_dict['Orientation'] = gene.orientation
            gene_dict['AvailableStarts'] = gene.alignment_candidate_start_nums
            gene_dict['AvailableCoord'] = gene.alignment_indices_to_coords_optimized(gene.alignment_candidate_starts)
            gene_dict['AvailableAlignmentIndices'] = list(gene.alignment_candidate_starts)
            gene_dict['DraftStatus'] = "Draft" if gene.draftStatus else "Final"

            if gene.locustag is None or gene.annot_author is None:
                gene.get_locustag()
            if gene.locustag != "":
                gene_dict['locustag'] = gene.locustag
            gene_dict['seaphage'] = "True" if gene.annot_author == 1 else "False"

            # Alignment-level data for visualization
            # gene_dict['AlignmentSequence'] = str(gene.alignment.seq)
            gene_dict['SeqRuns'] = _seq_runs_from_feature_runs(gene.alignment_feature_runs)
            gene_dict['CalledStartNum'] = gene.alignment_start_num_called
            if isinstance(gene.alignment_start_num_called, int) and 1 <= gene.alignment_start_num_called <= len(self.total_possible_starts):
                gene_dict['CalledStartAlignmentIndex'] = self.total_possible_starts[gene.alignment_start_num_called - 1]
            else:
                gene_dict['CalledStartAlignmentIndex'] = None
            gene_dict['AnnotatedStartNums'] = gene.alignment_annot_start_nums
            gene_dict['AnnotatedStartCounts'] = gene.alignment_annot_start_counts

            # Quality flags
            gene_dict['CallsMostAnnotated'] = gene.calls_most_annotated
            gene_dict['HasMostAnnotated'] = gene.has_most_annotated
            gene_dict['CalledStartIsBad'] = gene.called_start_is_bad
            gene_dict['BadAdjacentStartNums'] = gene.bad_adjacent_start_nums

            # Cluster/subcluster
            gene_dict['Cluster'] = gene.cluster
            gene_dict['Subcluster'] = gene.subcluster

            genelist.append(gene_dict)

        summary_dict['Genes'] = genelist
        annotlist = {}
        annotlist['Starts'] = list(self.stats['most_common']['annot_counts'].keys())
        annotlist['Counts'] = list(self.stats['most_common']['annot_counts'].values())
        summary_dict['Annots'] = annotlist

        conservationdict = {}
        for i in range(0, summary_dict['TotalStarts']):
            conservationdict[i+1] = round(float(len(self.stats["most_common"]['possible'][i+1]))/self.count, 4)
        summary_dict['Conservation'] = conservationdict

        # Gene groups: ordered list of lists of GeneIDs, matching PDF render order
        groups = self.group_similar_genes()
        summary_dict['GeneGroups'] = [[g.gene_id for g in group] for group in groups]

        # Pham-level most called and most annotated start numbers
        summary_dict['MostCalledStartNum'] = self.stats['most_common']['most_called_start']
        summary_dict['MostAnnotatedStartNum'] = self.stats['most_common']['most_annotated_start']

        # Alignment length (all genes share the same length after MAFFT)
        first_gene = next(iter(self.genes.values()))
        summary_dict['AlignmentLength'] = len(first_gene.alignment)

        return summary_dict

    def export_json(self, filename, compress=False):
        blob = self.annot_summary()
        data = orjson.dumps(blob, option=orjson.OPT_NON_STR_KEYS)
        if compress:
            with gzip.open(filename + '.gz', 'wb') as f:
                f.write(data)
        else:
            with open(filename, 'wb') as f:
                f.write(data)

    def add_cluster_stats(self, self_stats):
        clusters_present = set()
        genes_by_cluster = {}
        for g, pg in self.genes.items():
            clusters_present.add(pg.subcluster)
            if pg.subcluster in genes_by_cluster.keys():
                genes_by_cluster[pg.subcluster].append(g)
            else:
                genes_by_cluster[pg.subcluster] = [g]

        self_stats['clusters_present'] = [str(t) for t in clusters_present]
        self_stats['genes_by_cluster'] = genes_by_cluster

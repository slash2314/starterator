#!/usr/bin/env python
# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Starterate function 

import argparse
import json
import sys
from multiprocessing import Pool, Process, Queue, Semaphore
from .phams import (
    compare_hash_snapshots,
    compare_hashes_current,
    generate_pham_hashes,
    get_all_phams,
    process_all_phams,
    process_pham_list,
    process_phams_from_file,
)
from . import utils
from .utils import clean_up_files, StarteratorError
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter, A4
#from gi.repository import Gtk, Gdk, GObject
from . import report
from . import phamgene

"""
def gui():
    GObject.threads_init()
    Gdk.threads_init()
    win = StarteratorWindow()
    win.connect('delete-event', Gtk.main_quit)
    win.show_all()
    Gdk.threads_enter()
    Gtk.main()
    Gdk.threads_leave()
"""
def get_output_one_pham(pham, pham_no, config):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statisitics
    """
    output_dir = config['intermediate_file_dir']
    doc = SimpleDocTemplate("%s%sPham%sText.pdf" % (output_dir, phage+one_or_all, pham_no), pagesize=letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Pham %s Report </font>' % pham_no
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    output = pham.output_start_sites()
    for line in output:
        if line == '':
            story.append(Spacer(1, 12))
        text = '<font size=12> %s </font>' % line
        story.append(Paragraph(text, styles['Normal']))
    suggested_start = pham.output_suggested_starts()
    story.append(Spacer(1, 12))
    for line in suggested_start:
        text = '<font size=12>%s</font>' % line
        story.append(Paragraph(text, styles["Normal"]))
    doc.build(story)

def build_parser():
    parser = argparse.ArgumentParser(prog='starterate.py', usage='Phameratored Phage Report: %(prog)s -p {Phage Name}\n'
            + 'One Gene of Phameratored Phage Report:  %(prog)s -p {Phage Name} -n {Pham Number}\n'
            + 'Unphameratored Phage Report:  %(prog)s -p {Phage Name} -u True -f {Path to DNAMaster profile file}\n'
            + 'One Gene of Unphameratored Phage Report:  %(prog)s -p {Phage Name} -u True -s {Start of Gene} '
            + '-t {Stop of Gene} -o {Orientation of Gene} -g {Number of Gene}')
    parser.add_argument('-n', '--pham_no', default = -1,
                        help='Number of the Pham. For case when want report of a phameratored phage gene.')
    parser.add_argument('-p' , '--phage', default=None,
                        help='The Phamerator database Phage Name. Always needed')
    parser.add_argument('-u', '--unphamed', type=bool, default=False,
                        help='Boolean. If phage has been phameratored: False.'
                        +' If phage is unphameratored: True. For use when want report with an unphameratored phage')
    parser.add_argument('-s', '--given_start', type=int, default=-1,
                        help= 'The start of a gene. Use to report on one gene of an unphameratored phage.')
    parser.add_argument('-t', '--given_stop', type=int,
                        help= 'The stop of a gene. Use to report on one gene of an unphameratored phage.')
    parser.add_argument('-o', '--given_orientation',
                        help='The orientation of a gene. Use to report on one gene of an unphameratored phage.')
    parser.add_argument('-g', '--gene_number', default=-1,
                        help='The number of a gene. Use to report on one gene of an unphameratored phage.')
    parser.add_argument('-d', '--profile',
                        help='Path to a DNAMaster profile. For case when want whole report of an unphameratored phage')
    parser.add_argument('-f', '--fasta', help='Path to Fasta File')
    parser.add_argument('-j', '--save_json', type=bool, default=False,
                        help='Boolean, use with -n to save json file describing complete results.')
    parser.add_argument('--compress-json', action='store_true',
                        help='Compress JSON output as .json.gz (default JSON output is uncompressed).')
    parser.add_argument('--no-pdfs', action='store_true',
                        help='Skip PDF generation and only emit JSON output (pham report mode).')
    parser.add_argument('--all-phams', action='store_true',
                        help='Batch process all phams in the database')
    parser.add_argument('--with-phams', type=str, default=None,
                        help='Process phams listed in a line-delimited file (one pham ID per line)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('--get-phams', action='store_true',
                        help='Get all pham ids (one per line)')
    parser.add_argument('--get-pham-hashes', action='store_true',
                        help='Get all pham hashes')
    parser.add_argument('--compare-hash-files', nargs=2, metavar=('FILE1', 'FILE2'),
                        help='Compare two hash files')
    parser.add_argument('--delta-phams', type=str, default=None,
                        help='Process only added or modified phams compared to a previous hash file')
    return parser


def _get_selected_primary_modes(args):
    selected_modes = []
    if args.delta_phams:
        selected_modes.append('--delta-phams')
    if args.compare_hash_files:
        selected_modes.append('--compare-hash-files')
    if args.get_phams:
        selected_modes.append('--get-phams')
    if args.get_pham_hashes:
        selected_modes.append('--get-pham-hashes')
    if args.all_phams:
        selected_modes.append('--all-phams')
    if args.with_phams:
        selected_modes.append('--with-phams')
    if args.phage is not None:
        selected_modes.append('--phage')
    elif args.pham_no != -1:
        selected_modes.append('--pham_no')
    return selected_modes


def validate_primary_mode(args, parser):
    selected_modes = _get_selected_primary_modes(args)
    if len(selected_modes) > 1:
        parser.error("primary operation modes are mutually exclusive: %s" % ", ".join(selected_modes))


def get_arguments(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_primary_mode(args, parser)
    return args


def _load_previous_hash_snapshot(previous_hash_file):
    try:
        with open(previous_hash_file, 'r') as f:
            return json.load(f)
    except OSError as exc:
        raise StarteratorError(
            "Could not read previous pham hash file %s: %s" % (previous_hash_file, exc)
        )
    except ValueError as exc:
        raise StarteratorError(
            "Invalid JSON in previous pham hash file %s: %s" % (previous_hash_file, exc)
        )



def starterate(info, gui=None, event=None):
    global one_or_all, phage, protein_db, output_dir, final_dir
    one_or_all = 'All' if info['all'] else 'One'
    phage = info['phage']
    # print phage
    protein_db = utils.PROTEIN_DB + 'ProteinsDB'
    output_dir = utils.INTERMEDIATE_DIR
    final_dir = utils.FINAL_DIR
    if info['all'] and info['phamerated']:
        phage = report.PhageReport(phage, gui=gui, event=event)
        final_file, short_final = phage.final_report()
    elif info['all'] and not info['phamerated']:
        # clear intermediate files if doing this
        clean_up_files(utils.INTERMEDIATE_DIR)

        # phams.update_protein_db(db, config)
        phage = report.UnPhamPhageReport(phage, fasta_file=info['fasta'], profile_file=info['profile'],
                                         gui=gui, event=event)
        final_file, short_final = phage.final_report()
    elif not info['all'] and info['phamerated'] and not info['pham']:
        gene = report.GeneReport(info['phage'], number=info['gene_no'])
        gene.get_pham()
        gene.make_report()
        final_file, s = gene.merge_report()
    elif not info['all'] and not info['phamerated'] and not info['pham']:
        # phams.update_protein_db(db, config)
        gene = report.GeneReport(info['phage'], info["gene_no"], fasta_file=info['fasta'])
        gene.make_unpham_gene(int(info['start']), int(info['stop']), info['orientation'])
        gene.make_report()
        final_file, s = gene.merge_report()
    else:
        # Pham without phage'
        pham = report.PhamReport(info['pham'])
        final_file, s = pham.final_report()
    return final_file


import os
import logging
from . import utils

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main(argv=None):
    # Use the same configuration approach as utils.py
    config_path = os.path.abspath(os.path.join(
        os.getenv("STARTERATOR_CONFIG_DIR", os.path.join(os.environ["HOME"], ".starterator")), 
        "starterator.config"))
    
    # Configuration will be created automatically if it doesn't exist
    config = utils.get_config()
    args = get_arguments(argv)
    if args.verbose:
        # Log the content of the configuration file
        with open(config_path, 'r') as config_file:
            logging.info("Configuration file {} content:\n{}".format(config_path, config_file.read()))

    # Existing code
    config = utils.get_config()

    if args.get_phams:
        phams = get_all_phams()
        # phams returned one per line
        for pham in phams:
            print(pham)
        return
    
    if args.compare_hash_files:
        # Returns changed phams one per line
        results = compare_hashes_current(args.compare_hash_files[0], args.compare_hash_files[1])
        if not results["overall_hash_matches"]:
            for pham_id in results["phams_added"]:
                print(f"+{pham_id}")
            for pham_id in results["phams_modified"]:
                print(f"~{pham_id}")
            for pham_id in results["phams_removed"]:
                print(f"-{pham_id}")
        else:
            print("No differences found!")
        return

    if args.get_pham_hashes:
        generate_pham_hashes()
        return

    phamgene.check_protein_db(config["count"])

    if args.delta_phams:
        try:
            previous_hash = _load_previous_hash_snapshot(args.delta_phams)
            current_hash = generate_pham_hashes(output_file=True)
            results = compare_hash_snapshots(previous_hash, current_hash)
        except (KeyError, TypeError, StarteratorError) as e:
            logging.error(str(e))
            sys.exit(1)

        candidate_phams = set(results["phams_added"]) | set(results["phams_modified"])
        eligible_phams = [
            pham_id for pham_id in get_all_phams()
            if str(pham_id) in candidate_phams
        ]
        if not eligible_phams:
            print("No eligible changed phams found; nothing to process.")
            return
        process_pham_list(
            eligible_phams,
            no_pdfs=args.no_pdfs,
            compress_json=args.compress_json,
        )
        return

    if args.all_phams:
        process_all_phams(no_pdfs=args.no_pdfs, compress_json=args.compress_json)
        return

    if args.with_phams:
        process_phams_from_file(
            args.with_phams,
            no_pdfs=args.no_pdfs,
            compress_json=args.compress_json,
        )
        return


    try:
        # --Phamerated and only one gene
        if args.gene_number != -1 and args.phage is not None and args.unphamed is False:
            gene = report.GeneReport(args.phage, args.gene_number, True)
            # print gene
            gene.get_pham()
            gene.make_report()
            final_file, s = gene.merge_report()

        # --Unphameratored Phage with only one gene
        elif args.given_start > -1 and args.phage is not None and args.unphamed is True:
            # given start and stop coordinates and orientation
            one_or_all = 'One'
            given_start = args.given_start
            given_stop = args.given_stop
            given_orientation = args.given_orientation
            gene_name = args.phage + '_' + str(args.gene_number)
            gene = report.GeneReport(args.phage, args.gene_number, fasta_file=args.fasta)
            gene.make_unpham_gene(given_start, given_stop, given_orientation)
            gene.make_report()
            final_file, s = gene.merge_report()

        # --Phameratored or Unphameratored Phages with all genes
        elif args.pham_no == -1 and args.phage is not None and args.unphamed is False:
            phage = report.PhageReport(args.phage, gui=None)
            final_file, short_final = phage.final_report()

        elif args.pham_no == -1 and args.phage is not None and args.unphamed is True:
            phage = report.UnPhamPhageReport(args.phage, fasta_file=args.fasta, profile_file=args.profile, gui=None)
            final_file, short_final = phage.final_report()
        elif args.phage is None:
            pham = report.PhamReport(args.pham_no)
            save_json = args.save_json is True or args.no_pdfs
            final_file, short_final = pham.final_report(
                save_json=save_json,
                no_pdfs=args.no_pdfs,
                compress_json=args.compress_json,
            )
    except StarteratorError as e:
        logging.error(str(e))
        sys.exit(1)

    # clean_up_files()
    # email_final_report(args.email, short_final)



if __name__ == "__main__":
    main()

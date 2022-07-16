#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2021-12-30 -> 2022-07-15
Purpose: Performs BiSearch primer design for a given fasta or all .fasta files in a given directory.
"""

# TODO: Continue to refine the rate-limiter for BiSearch. Also improve documentation and continue refactoring.

import argparse
from typing import NamedTuple

import pathlib
import re
import random
import time
import pandas as pd
from Bio import SeqIO

try:
    import mechanicalsoup
except ImportError as e:
    print('MechanicalSoup is required for BiSearch WebScraping.')

class Args(NamedTuple):
    """ Command-line arguments """

    input_path: pathlib.Path
    output_path: pathlib.Path

    bis: bool
    strand: str
    max: int

    sub: int
    tries: int

    ignore_rate_lim: bool
    verbose: bool
    seed: int
# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Program parses .FASTA files and generates primers for bisulfite-converted template via BiSearch.',
        usage='%(prog)s [options] [input_path] [output_path]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'input_path',
        metavar='input_path',
        type=pathlib.Path,
        help='path of .FASTA file or directory containing .FASTA files')

    parser.add_argument(
        'output_path',
        metavar='output_path',
        type=pathlib.Path,
        help="path of output directory")

    group_BiSearch = parser.add_argument_group('BiSearch parameters')

    group_BiSearch.add_argument(
        '--bis',
        help="bisulfite-converted",
        action='store_true')

    group_BiSearch.add_argument(
        '--strand',
        help="specify target strand after bisulfite-conversion {sense= 'p', antisense = 'm'}",
        choices=['p', 'm'],
        type=str,
        default=None)

    group_BiSearch.add_argument(
        '--max',
        help="max length for PCR",
        metavar='<n>',
        type=int,
        default=400)

    group_subsequence = parser.add_argument_group('subsequence search options')

    group_subsequence.add_argument(
        '-s',
        '--sub',
        help='divide the sequence into <n> subsequences to pull tries from',
        metavar='<n>',
        type=int,
        default=1)

    group_subsequence.add_argument(
        '-t',
        '--tries',
        help='number of tries to attempt in each subsequences',
        metavar='<n>',
        type=int,
        default=1)

    group_debug = parser.add_argument_group('debug')

    group_debug.add_argument(
        '--ignore_rate_lim',
        help='ignore rate-limiters',
        action='store_true')

    group_debug.add_argument('-v',
        dest='verbose',
        action='count',
        default=0,
        help='verbose output for debugging')

    group_debug.add_argument('--seed',
        help='seed random',
        metavar='<n>',
        type=int,
        default=990318)

    args = parser.parse_args()

    if args.strand and not args.bis:
        parser.error('--strand is only applicable if bisulfite-converted.')

    if args.ignore_rate_lim:
        print(print_warning('Ignoring the rate-limiter risks blacklisting! And it will be completely your fault.'))

    return Args(args.input_path, args.output_path, args.bis, args.strand, args.max, args.sub, args.tries, args.ignore_rate_lim, args.verbose, args.seed)
# --------------------------------------------------
class BiSearch_condition:
    """
    A class to represent a set of conditions to be passed to BiSearch Primer Design. 

    Attributes (that you care about):
        args (Args): set of arguments to pass for handling verbosity and logging
        run_info (str): information stored to handle logging
        sequence (str): query sequence to search for primers
        bisulfite (bool): bisulfite conversion status, True = Bisulfite converted
        strand (str): only applicable after bisulfite conversion, 'p' = (+) strand, 'm' = (-) strand
        forward_start (int): forward primer, start of query region
        forward_end (int): forward primer, end of query region
        reverse_start (int): reverse primer, start of query region
        reverse_end (int): reverse primer, end of query region
        PCR_max_length (int): maximum length (nucleotides) of PCR amplicon to query
        primer_conc (float): concentration of primers (uM)
        potassium_conc (float): concentration of potassium (mM)
        magnesium_conc (float): concentration of magnesium (mM)
        results_list (int): maximum number of top results to show, best keep this at 10 unless you're absolutely desperate
        genome (str): genome to query for ePCR

    BiSearch_primer:
        A class to represent a primer retrieved from BiSearch Primer Design.

        do_BiSearch_ePCR():
            A function to perform BiSearch ePCR with a given set of primer information.
    
    do_BiSearch_primer_design():
        A function to perform BiSearch Primer Design with a given condition.
    """

    BiSearch_form_mapping = (
        'seq',      # query sequence
        'bis',      # bisulfite status (True = bisulfite-converted)
        'sens',     # strand (only available if bisulfite-converted)
        'fbeg',     # forward primer, beginning of primer query
        'fend',     # forward primer, end of primer query
        'rbeg',     # reverse primer, beginning of primer query
        'rend',     # reverse primer, end of primer query
        'len',      # maximum PCR length

        'pc',       # primer concentration (umol)
        'kc',       # potassium concentration (mmol)
        'mgc',      # magnesium concentration (mmol)

        'nstore',   # results list size (still sorted by score, not necessary)

        'db')       # database for ePCR
    def __init__(
            self,
            args: Args,                     # list of arguments to easily handle verbosity and logging
            run_info: str,                  # run_info (str) to handle logging
            sequence: str = '',             # query sequence to search for primers
            bisulfite: bool = False,        # bisulfite conversion status, True = bisulfited converted
            strand = None,                  # only applicable after bislfite conversion, 'p' = (+) strand, 'm' = (-) strand
            forward_start: int = None,      # forward primer, start of query region
            forward_end: int = None,        # forward primer, end of query region
            reverse_start: int = None,      # reverse primer, start of query region
            reverse_end: int = None,        # reverse primer, end of query region
            PCR_max_length: int = 400,      # maximum length (nucleotides) of PCR amplicon to query

            primer_conc: float = 2,         # concentration of primers (uM)
            potassium_conc: float = 1.0,    # concentration of potassium (mM)
            magnesium_conc: float = 1.5,    # concentration of magnesium (mM)

            results_list = 10,              # maximum number of top results to show, best keep this at 10

            genome = 'Bos taurus'):         # genome to query for ePCR

        self.args = args
        self.run_info = run_info

        self.sequence = str(sequence.seq)
        self.bisulfite = bisulfite

        # strand is only important if bisulfite converted. If bisulfite converted and no strand specified, default to 'p' strand
        if bisulfite and strand:
            self.strand = strand
        elif bisulfite and not strand:
            self.strand = 'p'
            # logging information
            current_action = print_runtime(f'No strand was given, defaulted to positive strand.')
            log_and_print(current_action, 2, args, run_info)
        elif not bisulfite:
            self.strand = None

        # everything else
        self.forward_start = forward_start if forward_start else 0
        self.forward_end = forward_end
        self.reverse_start = reverse_start
        self.reverse_end = reverse_end if reverse_end else len(self.sequence)
        self.PCR_max_length  = PCR_max_length

        self.primer_conc = primer_conc
        self.potassium_conc = potassium_conc
        self.magnesium_conc = magnesium_conc

        self.results_list = results_list

        self.genome = genome
        self.summary_conditions = (
            self.sequence,
            self.bisulfite,
            self.strand,

            self.forward_start,
            self.forward_end,
            self.reverse_start,
            self.reverse_end,

            self.PCR_max_length,
            self.primer_conc,
            self.potassium_conc,
            self.magnesium_conc,

            self.results_list,

            self.genome)

        self.primers_list = []
        
        # logging information
        current_action = print_runtime(f'Initialized BiSearch conditions for web service.')
        log_and_print(current_action, 2, args, run_info)
    def __repr__(self):
        return f'< {self.strand} (start-{self.forward_end}) -- ({self.reverse_start}-end) >'
    class BiSearch_primer:
        """
        A class to represent a BiSearch primer as retrieved by BiSearch Primer Design.

        Attributes (that you care about):
            basically everything from before
        
        do_BiSearch_ePCR():
            A function to perform BiSearch ePCR with a given set of primer information.
        """

        ePCR_form_mapping = (
            'fp',       # forward primer sequence
            'rp',       # reverse primer sequence
            'bis',      # bisulfite status
            'db',       # genome database
            'fpcrlen'   # max ePCR length
            )

        def __init__(
                self,
                args,
                run_info,
                strand,
                score,
                f_seq, f_length, f_GC, f_Tm,
                r_seq, r_length, r_GC, r_Tm,
                amp_start, amp_end, amp_length,
                CpG_count,
                f_self_anneal, f_self_end_anneal,
                r_self_anneal, r_self_end_anneal,
                pair_anneal, pair_end_anneal,
                primer_info_link
                ):

            self.args = args
            self.run_info = run_info
            self.strand = strand
            self.score = score
            self.f_seq = f_seq
            self.f_length = f_length
            self.f_GC = f_GC
            self.f_Tm = f_Tm
            self.r_seq = r_seq
            self.r_length = r_length
            self.r_GC = r_GC
            self.r_Tm = r_Tm
            self.amp_start = amp_start
            self.amp_end = amp_end
            self.amp_length = amp_length
            self.CpG_count = CpG_count
            self.f_self_anneal = f_self_anneal if f_self_anneal else 0
            self.f_self_end_anneal = f_self_end_anneal if f_self_end_anneal else 0
            self.r_self_anneal = r_self_anneal if r_self_anneal else 0
            self.r_self_end_anneal = r_self_end_anneal if r_self_end_anneal else 0
            self.pair_anneal = pair_anneal if pair_anneal else 0
            self.pair_end_anneal = pair_end_anneal if pair_end_anneal else 0
            self.primer_info_link = primer_info_link if primer_info_link else 0

            self.primer_binds = []
        def __repr__(self) -> None:
            return f'< str: {self.strand}, f: {self.f_seq}, r: {self.r_seq}, score: {self.score}, length: {self.amp_length}, CpG count: {self.CpG_count}>'
        def do_BiSearch_ePCR(self, max_pcr_len: int = 1000) -> tuple:
            """
            Function to automatically perform BiSearch ePCR and retrieve results.

            Parameters:
                self: the instance of the BiSearch primer
                max_pcr_len (int): maximum length to perform in silico PCR
            
            Returns:
                (tuple) of primer bind information and intended PCR product 
            """

            if not self.args.ignore_rate_lim:
                current_action = print_runtime(f'Performing necessary rate-limiting for BiSearch ePCR ...')
                log_and_print(current_action, 3, self.args, self.run_info)
                time.sleep(15)
                current_action = print_runtime(f'Rate-limiting done.')
                log_and_print(current_action, 3, self.args, self.run_info)

            # define the ePCR conditions
            # initialize browser instance for BiSearch ePCR
            browser = mechanicalsoup.StatefulBrowser(soup_config={'features': 'lxml'})
            browser.open("http://bisearch.enzim.hu/?m=genompsearch")
            browser.select_form('form[action="?run"]')

            # deals with the degenerate primers
            for r in (('Y', 'T'),
                    ('R', 'A')):
                self.f_seq = self.f_seq.replace(*r)
            for r in (('Y', 'T'),
                    ('R', 'A')):
                self.r_seq = self.r_seq.replace(*r)

            list(map(browser.form.set, self.ePCR_form_mapping, (self.f_seq, self.r_seq, True, 'Bos taurus', max_pcr_len)))
            current_action = print_runtime(f'Submitting request for BiSearch ePCR ...')
            log_and_print(current_action, 2, self.args, self.run_info)
            browser.submit_selected()

            #  find primer binds based on their position on the website
            # --------------------------------------------------

            # initialize primer binds
            self.primer_binds = [
                0,	# forward primer binds on the plus strand
                0,	# reverse primer binds on the plus strand
                0,	# forward primer binds on the minus strand
                0]	# reverse primer binds on the minus strand

            # define relative position of primer binds on strands
            pos_primer_matches = browser.get_current_page().find_all(text=re.compile('found based'))

            try:
                # for each in primer binds, find the values
                for i, _ in enumerate(self.primer_binds):
                    self.primer_binds[i] = pos_primer_matches[i].strip()[:-55] if len(pos_primer_matches[i].strip()[:-55]) > 0 else 0
                current_action = print_runtime(f'Retrieved primer bind values.')
                log_and_print(current_action, 1, self.args, self.run_info)
            except IndexError:
                # logging information
                current_action = print_warning(f'An error occurred in finding primer binds!\nTried to use  <F:{self.f_seq}> and <R:{self.r_seq}>.')
                log_and_print(current_action, 0, self.args, self.run_info)
                return None

            #  find  number of PCR products on the sense and antisense strand
            # --------------------------------------------------

            # define relative position of the PCR product number depending on how many PCR products are generated in the sense and antisense strand
            self.pos_PCR_prod_p = 1 if str(browser.get_current_page().find_all('h3')[3].find_next_siblings()[2]).strip() == '<br/>' else 2
            self.pos_PCR_prod_n = 1 if str(browser.get_current_page().find_all('h3')[6].find_next_siblings()[2]).strip() == '<br/>' else 2

            # number of PCR products is equal to the number of elements in those tables.
            self.PCR_products_p = len(browser.get_current_page().find_all('h3')[3].find_next_siblings()[self.pos_PCR_prod_p].find_all(text=re.compile('len')))
            self.PCR_products_n = len(browser.get_current_page().find_all('h3')[6].find_next_siblings()[self.pos_PCR_prod_n].find_all(text=re.compile('len')))


            #  find the preferred amplicon on either the sense or antisense strand
            # --------------------------------------------------

            # initialize  preferred amplicon's chromosome number and genomic position
            self.amp_chromosome_num = self.amp_genomic_start = self.amp_genomic_end = 'N/A'

            # define  relative web position of the preferred amplicon depending on the strand
            pos_amp = None
            if self.strand == 'p' and self.PCR_products_p > 0:
                pos_amp = 3
                pos_amp_str = self.pos_PCR_prod_p
            elif self.strand == 'm' and self.PCR_products_n > 0:
                pos_amp = 6
                pos_amp_str = self.pos_PCR_prod_n

            if pos_amp:
                try:
                    i_containing_amplicon = [i for i, s in enumerate(browser.get_current_page().find_all('h3')[pos_amp].find_next_siblings()[pos_amp_str].find_all(text=re.compile('len'))) if str(self.amp_length) in s][0]
                    web_container = str(browser.get_current_page().find_all('h3')[pos_amp].find_next_siblings()[pos_amp_str].find_all('a')[i_containing_amplicon]).split('?')[1].split(';')

                    self.amp_chromosome_num = web_container[0][4:-4]
                    self.amp_genomic_start = web_container[2][6:-4]
                    self.amp_genomic_end = web_container[3].split('"')[0][4:]
                
                    current_action = print_runtime(f'Retrieved amplicon coordinates ...')
                    log_and_print(current_action, 2, self.args, self.run_info)
                except:
                    current_action = print_warning(f'Error occurred in retrieving amplicon coordinates!')
                    log_and_print(current_action, 0, self.args, self.run_info)

            self.num_degen_bases = sum(map(self.f_seq.count, ['R','Y'])) + sum(map(self.r_seq.count, ['R','Y']))
            self.num_repeats = sum(map(self.f_seq.count, ['AA','TT', 'CC', 'GG'])) + sum(map(self.r_seq.count, ['AA','TT', 'CC', 'GG']))
            
            current_action = print_runtime(f'Received ePCR results for <F: {self.f_seq}> and <R: {self.r_seq}>.')
            log_and_print(current_action, 0, self.args, self.run_info)

            browser.close()
            current_action = print_warning(f'Closed browser instance of BiSearch ePCR.')
            log_and_print(current_action, 2, self.args, self.run_info)
            return (
                self.strand,
                self.score,

                self.f_seq,
                self.f_length,
                self.f_GC,
                self.f_Tm,

                self.r_seq,
                self.r_length,
                self.r_GC,
                self.r_Tm,

                self.amp_start,
                self.amp_end,
                self.amp_length,
                self.CpG_count,

                self.f_self_anneal,
                self.f_self_end_anneal,

                self.r_self_anneal,
                self.r_self_end_anneal,

                self.pair_anneal,
                self.pair_end_anneal,

                self.primer_info_link,

                self.primer_binds[0],
                self.primer_binds[1],
                self.PCR_products_p,
                self.primer_binds[2],
                self.primer_binds[3],
                self.PCR_products_n,

                self.amp_chromosome_num,
                self.amp_genomic_start,
                self.amp_genomic_end,

                self.num_degen_bases,
                self.num_repeats
            )
    def do_BiSearch_primer_design(self) -> list:
        """
        Function to automatically perform BiSearch Primer design and retrieve results.

        Parameters:
            self: the instance of the BiSearch conditions
        
        Returns:
            (list) of primers designed by BiSearch Primer design 
        """

        # initialize browser instance for BiSearch ePCR, select form, fill in form according to class variables
        browser = mechanicalsoup.StatefulBrowser(soup_config={'features': 'lxml'})
        
        # Maybe wait a bit before starting a new instance
        if not self.args.ignore_rate_lim:
            current_action = print_runtime(f'Performing necessary rate-limiting for initial primer design ...')
            log_and_print(current_action, 3, self.args, self.run_info)
            time.sleep(60)
            current_action = print_runtime(f'Rate-limiting done.')
            log_and_print(current_action, 3, self.args, self.run_info)

        current_action = print_runtime(f'Opened BiSearch browser instance.')
        log_and_print(current_action, 2, self.args, self.run_info)

        browser.open("http://bisearch.enzim.hu/?m=search")
        browser.select_form('form[action="?run"]')
        list(map(browser.form.set, self.BiSearch_form_mapping, self.summary_conditions, (True,)*len(self.summary_conditions)))
        current_action = print_runtime(f'Submitting request to BiSearch ...')
        log_and_print(current_action, 0, self.args, self.run_info)
        browser.submit_selected()

        # translate HTML into something that pandas can read -- essentially gets rid of unnecessary form settings
        BiSearchTable = browser.get_current_page().find_all('table')[11]
        HTML_string = str(BiSearchTable)
        for r in (('<form action="?run" method="post">', ''),
                ('</form>', ''),
                ('<input name="prg" type="hidden" value="cgi/fpcr.cgi"/>', ''),
                ('<tr class="r1">', '<tr>'),
                ('<tr class="r2">', '<tr>')):
            HTML_string = HTML_string.replace(*r)

        current_action = print_runtime(f'Received table of primer results from BiSearch.')
        log_and_print(current_action, 1, self.args, self.run_info)

        list_of_primer_info_links = []
        for tr in BiSearchTable.findAll('tr'):
            trs = tr.findAll('td')
            for each in trs:
                try:
                    if str(each.find('a')['href']).startswith('index.php'):
                        primer_info_link_str = 'http://bisearch.enzim.hu/'+ str(each.find('a')['href'])
                        list_of_primer_info_links.append(primer_info_link_str)
                        list_of_primer_info_links.append(' ')
                    else:
                        pass
                except:
                    pass

        current_action = print_runtime(f'Processed hyperlinks for primer-dimer information.')
        log_and_print(current_action, 1, self.args, self.run_info)

        BiSearch_Results = pd.read_html(HTML_string, skiprows=1)[0]

        current_action = print_runtime(f'Processing table of primer results ...')
        log_and_print(current_action, 1, self.args, self.run_info)

        # check if the results page is a legitimate table of primer results -- if there aren't 14 columns, there weren't any results
        if BiSearch_Results.shape[1] > 1:
            for row_i, _ in BiSearch_Results[::2].iterrows():
                x = BiSearch_Results.iloc

                try:
                    # This is just some webscraping laziness
                    score = x[row_i, 1]                     # primer score, as scored by BiSearch
                    strand = self.strand                    # strand
                    
                    f_seq = x[row_i, 2]                     # sequence of the foward primer
                    f_length = len(x[row_i, 2])             # length of the foward primer
                    f_GC = x[row_i, 5]                      # %GC of the foward primer
                    f_Tm = x[row_i, 6]                      # T_m of the foward primer

                    r_seq = x[row_i + 1, 2]                 # sequence of the reverse primer
                    r_length = len(x[row_i + 1, 2])         # length of the reverse primer
                    r_GC = x[row_i + 1, 5]                  # %GC of the reverse primer
                    r_Tm = x[row_i + 1, 6]                  # T_m of the reverse primer

                    amp_start = x[row_i, 3]                 # start of amplicon
                    amp_end = (x[row_i + 1, 3]) + 2         # end of amplicon
                    amp_length = (amp_end - amp_start) - 1  # length of amplicon
                    CpG_count = int(x[row_i, 8])            # number of CpG

                    try:
                        f_self_anneal = int(str(x[row_i, 9]).strip())           # forward primer self annealing
                        f_self_end_anneal = int(str(x[row_i, 10]).strip())      # forward primer self end-annealing

                        r_self_anneal = int(str(x[row_i + 1, 9]).strip())       # reverse primer self annealing
                        r_self_end_anneal = int(str(x[row_i + 1, 10]).strip())  # reverse primer self end-annealing
                    except:
                        current_action = print_warning(f'Failed to retrieve primer self-annealing information!')
                        log_and_print(current_action, 0, self.args, self.run_info)

                    try:
                        pair_anneal = int(str(x[row_i, 11]).strip())        # primer pair self annealing
                        pair_end_anneal = int(str(x[row_i, 12]).strip())    # primer pair self end-annealing
                    except:
                        current_action = print_warning(f'Failed to retrieve primer pair self-annealing information!')
                        log_and_print(current_action, 0, self.args, self.run_info)

                    try:
                        primer_info_link = list_of_primer_info_links[row_i]
                    except:
                        current_action = print_warning(f'Failed to process hyperlink information!')
                        log_and_print(current_action, 0, self.args, self.run_info)
                except:
                    current_action = print_warning(f'An error occurred in table processing!')
                    log_and_print(current_action, 0, self.args, self.run_info)

                    current_action = print_warning(
                        x[row_i, 1],
                        self.strand,
                        x[row_i, 2],
                        len(x[row_i, 2]),
                        x[row_i, 5],
                        x[row_i, 6],
                        x[row_i + 1, 2],
                        len(x[row_i + 1, 2]),
                        x[row_i + 1, 5],
                        x[row_i + 1, 6],
                        x[row_i, 3],
                        (x[row_i + 1, 3]) + 2,
                        (amp_end - amp_start) - 1,
                        int(x[row_i, 8]),
                        int(str(x[row_i, 9]).strip()),
                        int(str(x[row_i, 10]).strip()),
                        int(str(x[row_i + 1, 9]).strip()),
                        int(str(x[row_i + 1, 10]).strip()),
                        int(str(x[row_i, 11]).strip()),
                        int(str(x[row_i, 12]).strip()),
                    )
                    log_and_print(current_action, 0, self.args, self.run_info)

                try:
                    self.primers_list.append(
                        self.BiSearch_primer(self.args, self.run_info, strand, score, f_seq, f_length, f_GC, f_Tm, r_seq, r_length, r_GC, r_Tm, amp_start, amp_end, amp_length, CpG_count, f_self_anneal, f_self_end_anneal, r_self_anneal, r_self_end_anneal, pair_anneal, pair_end_anneal, primer_info_link)
                    )
                except:
                    current_action = print_warning(f'An error occurred in appending the primer to the list!')
                    log_and_print(current_action, 0, self.args, self.run_info)

        current_action = print_warning(f'Processed table of primer results.')
        log_and_print(current_action, 0, self.args, self.run_info)
        browser.close()
        current_action = print_warning(f'Closed browser instance of BiSearch Primer Design.')
        log_and_print(current_action, 2, self.args, self.run_info)
        return self.primers_list
# --------------------------------------------------
def main() -> None:
    """ Do the thing """

    args = get_args()

    # generate an output file
    output_path = pathlib.Path.resolve(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # run info
    run_info = time.strftime("%Y_%m_%d_%H%M%S", time.localtime(time.time()))
 
    # generate params file
    generate_log_file(output_path, run_info, args)

    # generate a list of files to iterate through
    if args.input_path.is_file():
        list_of_fasta_files = [args.input_path]
        # logging information
        current_action = print_runtime(f'Input path to file was used and found: {args.input_path.name}')
        log_and_print(current_action, 1, args, run_info)
    elif args.input_path.is_dir():
        list_of_fasta_files = [fasta_file for fasta_file in args.input_path.glob('*.fasta')]
        # logging information
        current_action = print_runtime(f'Input path to directory was used and found: {[fasta_file.name for fasta_file in list_of_fasta_files]}')
        log_and_print(current_action, 1, args, run_info)

    compiled_primers_total = {}
    for fasta_file in list_of_fasta_files:
        # logging information
        current_action = print_runtime(f'Generating BiSearch conditions for {fasta_file.name} ...')
        log_and_print(current_action, 0, args, run_info)

        query_sequence = SeqIO.read(fasta_file, "fasta")
        fasta_subseqs = generate_subseqs(query_sequence, args.sub, args.max+200)
        
        # logging information
        current_action = print_runtime(f'Generating subsequences to try for {fasta_file.name}...')
        log_and_print(current_action, 1, args, run_info)

        subsequence_list = []
        for subseq_i in fasta_subseqs.keys():
            tries_arg = args.tries
            if tries_arg > len(fasta_subseqs[subseq_i].values()):
                # logging information
                current_action = print_warning(f'Splitting the FASTA into {tries_arg} tries in {len(fasta_subseqs.keys())} x {args.max+200} bp subsequences isn\'t useful. The query region will overlap into the other regions anyway.')
                log_and_print(current_action, 2, args, run_info)
                
                while tries_arg > len(fasta_subseqs[subseq_i].values()):
                    tries_arg -= 1
                current_action = print_runtime(f'Rounded down to {tries_arg} tries in {len(fasta_subseqs.keys())} regions to avoid this.')
                log_and_print(current_action, 2, args, run_info)
            subsequence_list += random.choices(list(fasta_subseqs[subseq_i].values()), k=tries_arg)
        # logging information
        current_action = print_runtime(f'Generated {len(subsequence_list)} stepwise subsequence(s) to try.')
        log_and_print(current_action, 1, args, run_info)
    
        list_of_conditions = []
        for subsequence in subsequence_list:
            if int(subsequence['start']) < 20:
                subsequence['start'] = subsequence['start'] + 20
                subsequence['end'] = subsequence['end'] + 20
            if args.bis and not args.strand:
                try:
                    pos_strand_conditions = BiSearch_condition(
                        args=args,
                        run_info=run_info,
                        sequence = query_sequence,
                        reverse_start = subsequence['start'],
                        forward_end = subsequence['end'],
                        bisulfite = args.bis,
                        strand = 'p',
                        )
                    neg_strand_conditions = BiSearch_condition(
                        args=args,
                        run_info=run_info,
                        sequence = query_sequence,
                        reverse_start = subsequence['start'],
                        forward_end = subsequence['end'],
                        bisulfite = args.bis,
                        strand = 'm',
                        )
                    current_action = print_runtime(f'Generated BiSearch conditions for the (+) and (-) strands.')
                    log_and_print(current_action, 0, args, run_info)

                    current_action = print_runtime(f'Trying the conditions for the (+) strand ...')
                    log_and_print(current_action, 0, args, run_info)
                    condition = pos_strand_conditions.do_BiSearch_primer_design()
                    list_of_conditions.append(condition)
                    
                    current_action = print_runtime(f'Trying the conditions for the (-) strand ...')
                    log_and_print(current_action, 0, args, run_info)
                    condition = neg_strand_conditions.do_BiSearch_primer_design()
                    list_of_conditions.append(condition)

                except IndexError:
                    current_action = print_warning(f'An error occurred when processing BiSearch conditions for (+) and (-) strands!')
                    log_and_print(current_action, 0, args, run_info)
            elif args.bis and args.strand:
                try:
                    condition = BiSearch_condition(
                        args=args,
                        run_info=run_info,
                        sequence = query_sequence,
                        reverse_start = subsequence['start'],
                        forward_end = subsequence['end'],
                        bisulfite = args.bis,
                        strand = args.strand,
                        )
                    current_action = print_runtime(f'Generated BiSearch conditions for the {args.strand} strand.')
                    log_and_print(current_action, 0, args, run_info)

                    current_action = print_runtime(f'Trying the current BiSearch condition ...')
                    log_and_print(current_action, 0, args, run_info)
                    list_of_conditions.append(condition.do_BiSearch_primer_design())
        
                except IndexError:
                    current_action = print_warning(f'An error occurred when processing BiSearch conditions for {args.strand} strand!')
                    log_and_print(current_action, 0, args, run_info)
    
        compiled_primers_for_fasta = []
        for primer_condition in list_of_conditions:
            for BiSearch_primer in primer_condition:
                if args.verbose > 2 : print_runtime(f'Performed BiSearch ePCR with {BiSearch_primer}')
                compiled_primers_for_fasta.append(BiSearch_primer.do_BiSearch_ePCR())

        compiled_primers_total.update({fasta_file.stem: compiled_primers_for_fasta})

        generate_primers_csv(compiled_primers_for_fasta, fasta_file.stem, output_path)
        current_action = print_runtime(f'Generated .csv files containing primers for {fasta_file}.fasta .')
        log_and_print(current_action, 0, args, run_info)
    
    current_action = print_runtime(f'Finished job.')
    log_and_print(current_action, 0, args, run_info)
def generate_primers_csv(data, csv_name, output_path) -> None:
    """
    Function outputs BiSearch-generated primers to .csv
    """
    pd.DataFrame(data,
            columns=[
                'strand',           # primer strand, 'p' = sense, 'm' = antisense
                'score',            # primer score, as scored by BiSearch

                'f_seq',            # sequence of forward primer
                'f_length',         # length of forward primer
                'f_GC',             # %GC of forward primer
                'f_Tm',             # T_m of forward primer

                'r_seq',            # sequence of reverse primer
                'r_length',         # length of reverse primer
                'r_GC',             # %GC of reverse primer
                'r_Tm',             # T_m of reverse primer

                'amp_start',        # start of amplicon
                'amp_end',          # end of amplicon
                'amp_size',         # size of amplicon
                'CpG_count',        # number of CpG

                'f_self_anneal',    # self annealing of the forward primer
                'f_self_end_anneal',# self end-annealng of the forward primer

                'r_self_anneal',    # self annealing of the reverse primer
                'r_self_end_anneal', # self end-annealing of the reverse primer

                'pair_anneal',      # pair annealing
                'pair_end_anneal',  # pair end-annealing

                'primer_info_link', # the link containing the primer info

                'f_bind (+)',       # number of forward primer binds on (+) strand
                'r_bind (+)',       # number of reverse primer binds on (+) strand
                'PCR_products (+)', # number of PCR products on (+) strand

                'f_bind (-)',       # number of forward primer binds on (-) strand
                'r_bind (-)',       # number of reverse primer binds on (-) strand
                'PCR_products (-)', # number of PCR products on (-) strand

                'chr',              # chromosome number
                'start',            # start pos in chromosome
                'end',              # end pos in chromosome

                'degen_bases',      # number of degenerate bases in the primer sequences
                'repeats'           # number of repeats in the primer sequences
                ]).to_csv(output_path.joinpath(f'{csv_name}.csv'))
def generate_subseqs(sequence, subseq_num, frame_length) -> dict:
    """
    Function will find subsequences of a given length within a fasta file.
    """
    subseq_len = int(len(sequence)/subseq_num)
    subseq_data = {}
    frame_length = subseq_len if frame_length > subseq_len else frame_length

    for subseq_count in range(subseq_num):
        subseq_list = {}

        subseq_start = subseq_count * subseq_len
        subseq_end = (subseq_count + 1) * subseq_len

        i = subseq_start
        while (i + frame_length) <= subseq_end:
            i += 1
            subseq_list.update({i: {
            'start': i,
            'end': i + frame_length,
            'CpG_count': sequence.seq[i:i+frame_length].count('CG')}})

        subseq_data.update({subseq_count: subseq_list})

    return subseq_data
def log_and_print(action_arg, verbosity_arg, args, run_info) -> None:
    """ Log the current action and print it if the verbosity is high enough """
    with open(pathlib.Path.resolve(args.output_path).joinpath(f'log-{run_info}.txt'), 'a', encoding='UTF-8') as log_file:
        log_file.write(f'\n{action_arg}')
    if args.verbose > verbosity_arg - 1 : print(action_arg)
def generate_log_file(output_path, run_info, args) -> None:
    """ Generates log file. """
    with open(output_path.joinpath(f'log-{run_info}.txt'), 'w', encoding='UTF-8') as log_file:
        output_vars = [
            f'input={args.input_path.resolve()}',
            f'output_dir={args.output_path.resolve()}',
            f'bis={args.bis}',
            f'strand={args.strand}',
            f'max_PCR_len={args.max}',
            f'n_subseq={args.sub}',
            f'n_tries={args.tries}',
            f'ignore_rate_lim={args.ignore_rate_lim}',
            f'verbose={args.verbose}',
            f'seed={args.seed}']
        log_file.write('\n'.join(output_vars))
        log_file.write('\n')
    print(print_runtime(f"Created log file in {output_path} ."))
def print_warning(action) -> None:
    """ Return the time and a warning. """
    return f'\n[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] !!! WARNING !!! :\n{action}'
def print_runtime(action) -> None:
    """ Return the time and some defined action. """
    return f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}'
# --------------------------------------------------
if __name__ == '__main__':
    main()

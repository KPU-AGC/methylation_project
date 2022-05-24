#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-05-20 -> 2022-05-23
Purpose: To manually screen JSON files produced from tracy variant calling.
"""


# TODO: SNV are still the most reliable to screen. Correctly highlighting indels is still difficult.
#       Screening indels is even more difficult. A system exists to compare it to the reference sequence
#       below the plot, but this implementation is currently flawed. There does not seem to be a reliable
#       way to get the flanking region from the reference without a substring search.


import argparse
from typing import NamedTuple, TextIO

import pathlib
import json

import pprint
import time

import matplotlib.pyplot as plt

from Bio.Seq import Seq

class Args(NamedTuple):
    """ Command-line arguments """
    target_dir: pathlib.Path

    qc_flag: bool
    all_vars: bool

    hide_nuc: bool
    hide_tracy: bool

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Open a .JSON file produced from Tracy variant calling and screen the variants.',
        epilog='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'target_dir',
        type=pathlib.Path,
        metavar='path',
        help='directory containing .JSON files/standalone .JSON file')

    parser.add_argument(
        '-q',
        '--qc_flag',
        action='store_true',
        help='include LowQual variants')

    parser.add_argument(
        '-a',
        '--all_vars',
        action='store_true',
        help='include all types of variants')

    group_UX_options = parser.add_argument_group('UX options (to remove bias?)')

    group_UX_options.add_argument(
        '-n',
        '--hide_nuc',
        action='store_true',
        help='hide nucleotide labels')

    group_UX_options.add_argument(
        '-t',
        '--hide_tracy',
        action='store_true',
        help='hide tracy labels')

    args = parser.parse_args()

    return Args(args.target_dir, args.qc_flag, args.all_vars, args.hide_nuc, args.hide_tracy)

# --------------------------------------------------
class JSON_data():
    """
    A class to represent data from the JSON files.

    Attributes (that you care about):
        orig_JSON_data (dict): original parsed JSON data for comparison purposes, not to be updated
        JSON_data (dict): parsed JSON data for updating variants
        sample_name (str): input file used for tracy variant calling as listed in JSON file
        seq (str): input
    """

    nucleotide_data = {
        'A': {'color': 'green', 'complement': 'T'},
        'T': {'color': 'red', 'complement': 'A'},
        'C': {'color': 'blue', 'complement': 'G'},
        'G': {'color': 'black', 'complement': 'C'},
    }


    def __init__(self, path_arg, output_arg, hide_nuc_arg, hide_tracy_arg):
        """
        Construct a JSON data class.

        Parameters:
            path_arg (pathlib.Path): path of the JSON file
            output_arg (pathlib.Path): path of the output file for this JSON file
            hide_nuc_arg (bool): if True, hide nucleotides when displaying displaying peaks
            hide_tracy_arg (bool): if True, hide tracy data when displaying peaks
        """
        self.hide_tracy_arg = hide_tracy_arg
        self.hide_nuc_arg = hide_nuc_arg
        self.output_path = output_arg

        with open(path_arg, 'r', encoding='UTF-8') as input_JSON:
            # raw JSON data
            self.orig_JSON_data = json.load(input_JSON)
        with open(path_arg, 'r', encoding='UTF-8') as input_JSON:
            # raw JSON data
            self.JSON_data = json.load(input_JSON)

        # sample name
        self.sample_name = self.JSON_data['meta']['arguments']['input']

        # seq
        self.seq = self.JSON_data['primarySeq']

        # channels
        self.channels = {}
        for nuc in self.nucleotide_data:
            self.channels[nuc] = self.JSON_data['peak' + nuc]

        # basecall positions
        self.basecall_pos = self.JSON_data['basecallPos']


    def __repr__(self) -> str:
        return f'<JSON data for: {self.sample_name}>'


    def get_variants(self) -> list:
        """
        Returns a list of list of variants in the JSON data.

        Parameters:
            None

        Returns:
            List of list of variants in JSON data
        """
        return self.JSON_data['variants']['rows']


    def display_peaks(self, variant_value_arg: list, nuc_flank: int = 5) -> None:
        """
        Display the peaks of a given position using matplotlib graph.

        Parameters:
            variant_value_arg (list): list containing parameters for the query variant
            nuc_flank (int): the number of nucleotides to flank the current position.

        Returns:
            None
        """

        #basecall_pos_text = variant_value_arg[9]
        basepos = variant_value_arg[9]
        rev_com_var = 1
        genotype_text = variant_value_arg[8][:4].upper()
        signalpos = variant_value_arg[10]
        tracy_variant_allele = variant_value_arg[4]


        # dynamic plot options
        list_of_peak_pos = [int(i) for i in self.JSON_data['basecalls'].keys()]
        average_peak_distance = sum([list_of_peak_pos[i + 1] - list_of_peak_pos[i] for i in range(len(list_of_peak_pos) - 1)])/len(list_of_peak_pos)
        #max_peak_height = max([max(peak_heights)for peak_heights in [self.channels[channel] for channel in self.channels]])
        max_peak_height = 1500

        # display params
        display_params = {
            'lower buffer': max_peak_height,
            'nucleotide': max_peak_height + 200,
            'basepos': max_peak_height + 200,
            'genotype': max_peak_height + 300,
            'buffer': max_peak_height + 600
        }

        # fix the spacing depending on hidden vars
        if self.hide_nuc_arg: display_params['buffer'] = display_params['buffer'] - 100
        if self.hide_tracy_arg: display_params['buffer'] = display_params['buffer'] - 300


        # plot options
        plt.title(self.sample_name)
        plt.rcParams['toolbar'] = 'None'

        #plt.subplots_adjust(bottom = 0.2)

        ax = plt.gca()
        plt.cla()
        ax.axis('off')
        ax.axes.xaxis.set_ticks([])
        ax.axes.yaxis.set_ticks([])


        # check reverse complement
        if self.JSON_data['ref1forward'] == 1:
            reverse_complement = False
        elif self.JSON_data['ref1forward'] == 0:
            reverse_complement = True


        # deal with reverse complementation
        if reverse_complement:

            # nucleotide channels are set to complementary nucleotide channel and reversed
            for nuc in self.nucleotide_data:
                self.JSON_data[f'peak{nuc}'] = self.orig_JSON_data[f"peak{self.nucleotide_data[nuc]['complement']}"][::-1]


            # sequence reverse complemented
            self.JSON_data['primarySeq'] = str(Seq(self.orig_JSON_data['primarySeq']).complement())

            # re-number the basecall positions
            for i, _ in enumerate(self.orig_JSON_data['basecallPos']):
                self.JSON_data['basecallPos'][i] = (max(self.orig_JSON_data['pos']) - self.orig_JSON_data['basecallPos'][i])

            signalpos = self.JSON_data['basecallPos'][basepos - 1]

            rev_com_var = -1


        # plot the channels
        for nuc in self.nucleotide_data:
            ax.plot(self.JSON_data[f'peak{nuc}'], color=self.nucleotide_data[nuc]['color'], label=nuc, linewidth=0.8)
            ax.legend(ncol=int(len(self.nucleotide_data)/2))


        # plot bases above channels
        if not self.hide_nuc_arg:
            for nuc_i, nuc_val in enumerate(self.JSON_data['primarySeq']):
                if nuc_val in self.nucleotide_data: nuc_col = self.nucleotide_data[nuc_val]['color']
                else: nuc_col ='purple'

                # plot the nucleotide text above
                plt.text(self.JSON_data['basecallPos'][nuc_i], display_params['nucleotide'], nuc_val, horizontalalignment='center', color=nuc_col, clip_on=True).set_in_layout(False)


        # highlight variant nucleotide position
        if reverse_complement: left_offset = 0
        else: left_offset = 2
        left_nuc_peak = self.JSON_data['basecallPos'][basepos - left_offset]
        left_distance = (self.JSON_data['basecallPos'][basepos - left_offset] - self.JSON_data['basecallPos'][basepos - 1])/2
        left_edge = left_nuc_peak - left_distance

        if reverse_complement: right_offset = 2
        else: right_offset = 0
        right_count = len(tracy_variant_allele) - 1
        right_nuc_peak = self.JSON_data['basecallPos'][basepos - right_offset - rev_com_var*(right_count)]
        right_distance = (self.JSON_data['basecallPos'][basepos - right_offset - rev_com_var*(right_count)] -  self.JSON_data['basecallPos'][basepos - 1])/2
        right_edge = right_nuc_peak - right_distance
        bar_width = right_edge - left_edge

        ax.bar(left_edge, display_params['nucleotide']*2, width=bar_width, align='edge', color = 'c', alpha=0.1)


        # plot tracy information above
        #plt.text(signalpos, display_params['basepos'], basepos, horizontalalignment='center')

        plt.text(signalpos, display_params['basepos']-100, basepos, horizontalalignment='center')
        if not self.hide_tracy_arg:

            plt.text(signalpos, display_params['genotype'], genotype_text, horizontalalignment='center')
            plt.text(signalpos, display_params['genotype']+100, variant_value_arg[7], horizontalalignment='center')


            # show reference sequence below
            sample_subseq = self.seq[basepos - 1 - nuc_flank:basepos - 1]
            ref_offset =  self.JSON_data['ref1align'].find(sample_subseq)
            str_subseq = self.JSON_data['ref1align'][ref_offset : ref_offset + (nuc_flank * 2)]
            if reverse_complement:
                str_subseq = Seq(str_subseq).reverse_complement()
            reference_text = ''.join([str(i + ' ') for i in str_subseq])
            plt.xlabel(f"Reference: {reference_text}")

        # limit view to current nucleotide window
        plt.xlim(signalpos - (average_peak_distance * nuc_flank), signalpos + (average_peak_distance * nuc_flank))
        plt.ylim(0, display_params['buffer'])

        plt.tight_layout()
        plt.show(block=False)


    def validate_nucleotide_arg(self, variant_allele_arg: str = '') -> list:
        """
        Validate the nucleotide that was passed into this function.

        Parameters:
            variant_allele_arg (str): argument that was passed into this function

        Returns:
            (list-ish): return the list of alleles if valid, otherwise False
        """

        # special cases:
        # if nothing was entered, skip
        if variant_allele_arg == '':
            return 'SKIP'

        # if * was entered, flag for resequencing
        if variant_allele_arg == '*':
            return 'SKIP_ERROR'

        # ex: deal with A/ or /A or A cases and return as just [A]
        variant_alleles = variant_allele_arg.upper().split('/')
        if '' in variant_alleles: variant_alleles.remove('')

        # triple alleles not allowed
        # ex: A/A/A is disallowed
        if len(variant_alleles) > 2: return False

        # invalid characters not allowed
        # ex: APOM is disallowed
        for allele in variant_alleles:
            for nucleotide in allele:
                if nucleotide not in ['A', 'T', 'C', 'G']: return False

        # homozygous in heterozygous format, somewhat allowed
        # ex: A/A = A
        if len(variant_alleles) == 2 and variant_alleles[0] == variant_alleles[1]:
            return [variant_alleles[0]]

        # allowed, return alleles
        if len(variant_alleles) == 2:
            return [variant_alleles[0], variant_alleles[1]]
        return variant_alleles


    def update_variant(self, variant_i_arg: int, variant_nucleotide_arg: list, update_type: str) -> None:
        """
        Update the list of variants in the working JSON file.

        Parameters:
            variant_i_arg (int): the integer position of the variant in the list of variants in the JSON
            variant_nucleotide_arg (list): list of alleles that passed validate_nucleotide_arg()
            update_type (str): the type of update, see cases below

        Returns:
            None

        --------------------------------------------------
        Update types:
            PASS: the screening passed tests, update directly
            SKIP_ERROR: the screening was unsuccessful, flag for resequencing
        """

        # deal with the * sequence flag
        if update_type == 'SKIP_ERROR':
            self.JSON_data['variants']['rows'][variant_i_arg][4] = '*'
            self.JSON_data['variants']['rows'][variant_i_arg][5] = '-'
            self.JSON_data['variants']['rows'][variant_i_arg][6] = 'ERR_SCREENED'
            self.JSON_data['variants']['rows'][variant_i_arg][7] = '-'
            self.JSON_data['variants']['rows'][variant_i_arg][8] = '-'
            prettied_json =  pprint.pformat(self.JSON_data, compact=False, sort_dicts=False).replace("'", '"')
            with open(self.output_path.joinpath(f'{self.sample_name}_s.json'), 'w') as new_JSON:
                new_JSON.write(prettied_json)
            return None

        # handle heterozygous input
        if isinstance(variant_nucleotide_arg, list):
            genotype='het.'
            for allele in variant_nucleotide_arg:
                if allele == self.JSON_data['variants']['rows'][variant_i_arg][3]:
                    variant_nucleotide_arg.remove(allele)
        else:
            genotype='hom. ALT'

        # update JSON data
        if variant_nucleotide_arg:
            self.JSON_data['variants']['rows'][variant_i_arg][4] = variant_nucleotide_arg[0]
            self.JSON_data['variants']['rows'][variant_i_arg][6] = 'SCREENED'
        else:
            self.JSON_data['variants']['rows'][variant_i_arg][4] = self.JSON_data['variants']['rows'][variant_i_arg][3]
            self.JSON_data['variants']['rows'][variant_i_arg][6] = 'ERR_SCREENED'
        self.JSON_data['variants']['rows'][variant_i_arg][5] = '-'
        self.JSON_data['variants']['rows'][variant_i_arg][8] = genotype

        # write to screened JSON file
        prettied_json =  pprint.pformat(self.JSON_data, compact=False, sort_dicts=False).replace("'", '"')
        with open(self.output_path.joinpath(f'{self.sample_name[:-4]}_s.json'), 'w') as new_JSON:
            new_JSON.write(prettied_json)


    def translate_variant(self, variant_value_arg: list) -> str:
        """
        Check the variant genotype and translate the alleles into allele/allele format for easier viewing.
        I have decided the convention is reference/variant for heterozygotes.

        Parameters:
            variant_value_arg (list): list of variant data as listed in JSON data

        Output:
            (str): translated string in allele/allele format

        --------------------------------------------------
        Example:
            If 'C' was the reference allele,
                a het. sample with the 'A' variant allele would return C/A.
                a hom. sample with the 'A' variant allele would return A/A.
        """
        if variant_value_arg[8][:4] == 'het.':
            return f'{variant_value_arg[3]}/{variant_value_arg[4]}'
        elif variant_value_arg[8][:4] == 'hom.':
            return f'{variant_value_arg[4]}/{variant_value_arg[4]}'


    def screen_variants(self, qc_flag_arg: bool, all_vars_arg: bool) -> None:
        """
        Main function for screening variants.

        Parameters:
            qc_flag_arg (bool): If True, also show LowQual variants
            all_vars_arg (bool): If True, show every kind of variant for screening, more info below

        Returns:
            None

        --------------------------------------------------
        all_vars_arg:
            Only SNV are shown unless the all_vars_arg flag is show. The validity of discerning indels
            given the electropherogram is still up for debate. The current system allows for rough inspection
            via the reference sequence at the bottom of the display. This requires further discussion.
        """

        for variant_i, variant_value in enumerate(self.get_variants()):
            if (variant_value[7] == 'SNV')  or (all_vars_arg):
                if (variant_value[6] == 'PASS') or (variant_value[6] == 'LowQual' and qc_flag_arg):
                    self.display_peaks(variant_value)

                    check_validated = False
                    while not check_validated:
                        plt.draw()
                        if not self.hide_tracy_arg: input_nucleotide = input(f'{variant_value[9]} {variant_value[7]} ({self.translate_variant(variant_value)}): ')
                        else: input_nucleotide = input(f'{variant_value[9]}: ')
                        if input_nucleotide == 'i': print_instructions()
                        check_validated = self.validate_nucleotide_arg(input_nucleotide)
                        if not check_validated: print('Enter a valid nucleotide.')
                    if check_validated == 'SKIP':
                        pass
                    elif check_validated == 'SKIP_ERROR':
                        self.update_variant(variant_i, self.validate_nucleotide_arg(input_nucleotide), update_type='SKIP_ERROR')
                    else:
                        self.update_variant(variant_i, self.validate_nucleotide_arg(input_nucleotide), update_type='PASS')
                else:
                    pass


def print_instructions() -> None:
    """ Print instructions for manual base-calling"""

    print(\
    """
    IMPORTANT INSTRUCTIONS:
    To perform manual base-calling, enter the nucleotide peaks that are present at that variant position.

    Allowed nucleotides: A, T, C, G

    Examples:
        a/a = homozygous A
        a/t = heterozygous A/T
        *   = ambiguous, impossible to call

    Press enter without entering a nucleotide to skip this position.
    Enter 'i' to bring up these instructions again.
    """)

def generate_reseq_list(output_path_arg: pathlib.Path) -> None:
    """
    Iterate through the generated ouput and create a list in a .txt file for variants with the re-sequence flag (*)

    Parameters:
        output_path_arg (pathlib.Path): the path of the output directory

    Returns: None
    """

    resequencing_list = []
    for JSON_file in output_path_arg.glob('*.json'):
        with open(JSON_file, 'r') as input_JSON:
            query_JSON_data = json.load(input_JSON)
            for variant_data in query_JSON_data['variants']['rows']:
                if variant_data[4] == '*': resequencing_list.append(JSON_file.stem)

    resequencing_list = sorted(set(resequencing_list))
    with open(output_path_arg.joinpath('reseq_report.txt'), 'w', encoding='UTF-8') as text_file:
        text_file.write(f'Good luck, future sequencing monkey!\n\n')
        text_file.write(f'To be re-sequenced:\n')
        for i in resequencing_list:
            text_file.write(f'{i}\n')
        text_file.write(f'\nNote: This is not a conclusive list. This does not include failed sequences.\n')
    print_runtime('Generated .txt report for re-sequencing')

def print_runtime(action) -> None:
    """Prints the runtime."""
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

# --------------------------------------------------
def main() -> None:
    """ docstring """

    args = get_args()
    target_path = pathlib.Path.resolve(args.target_dir)

    print_instructions()

    # create output folder
    if target_path.is_file(): output_path = target_path.parent.joinpath('tracy-variant-screening')
    elif target_path.is_dir(): output_path = target_path.joinpath('tracy-variant-screening')
    output_path.mkdir(parents=True, exist_ok=True)

    if target_path.exists() and target_path.is_file():
        print_runtime(f'Performing Tracy screening on {target_path.name}')
        query_JSON = JSON_data(target_path, output_path, args.hide_nuc, args.hide_tracy)
        query_JSON.screen_variants(args.qc_flag, args.all_vars)
        print_runtime(f'Finished Tracy screening on {target_path.name}')

    if target_path.exists() and target_path.is_dir():
        for JSON_file in target_path.glob('*.json'):
            print_runtime(f'Performing Tracy screening on {JSON_file.name}')
            query_JSON = JSON_data(JSON_file, output_path, args.hide_nuc, args.hide_tracy)
            query_JSON.screen_variants(args.qc_flag, args.all_vars)

    if target_path.is_dir():
        generate_reseq_list(output_path)

# --------------------------------------------------
if __name__ == '__main__':
    main()
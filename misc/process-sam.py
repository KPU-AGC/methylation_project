from pathlib import Path
import subprocess
import pandas as pd
import sys
from collections import Counter

primers_dict = {
    "DMAP1": ("NC_037330.1:101832510-101832748", True),
    "DNMT3A": ("NC_037338.1:74014960-74015346", True),
    "DNMT3B": ("NC_037340.1:62166550-62166873", True),
    "GNAS": ("NC_037340.1:57520605-57520891", False),
    "H19": ("NC_037356.1:49504607-49504946", True),
    "IGF2R": ("NC_037336.1:96223067-96223367", False),
    "KCNQ1": ("NC_037356.1:48907634-48907882", True),
    "LIF": ("NC_037344.1:69264000-69264354", False),
    "LIFR": ("NC_037347.1:35949053-35949420", True),
    "MEST": ("NC_037331.1:94249893-94250283", False),
    "NNAT": ("NC_037340.1:66465753-66466057", True),
    "PEG10": ("NC_037331.1:12063279-12063673", True),
    "PEG3": ("NC_037345.1:64120688-64120941", False),
    "PLAGL1": ("NC_037336.1:81418732-81419041", False),
    "RTL1": ("NC_037348.1:65778329-65778707", False),
    "SLC2A8": ("NC_037338.1:98154376-98154714", True),
    "SNRPN": ("NC_037348.1:1937284-1937513", True),
    "SUV39H1": ("NC_037357.1:86781615-86781929", False),
    "TXNIP": ("NC_037330.1:21383423-21383784", True),
    "XIST": ("NC_037357.1:77162729-77163121", True),
}

def reverse_complement(input_seq: str) -> str:
    """
    """

    nucleotides = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N'
    }

    return ''.join([nucleotides[nuc] for nuc in input_seq.upper()[::-1]])

def find_CpG(input_seq: str) -> list:
    """
    """

    positions = []
    for i, _ in enumerate(input_seq[:-1]):
        if input_seq[i] in 'Cc' and input_seq[i+1] in 'Gg':
            positions.append(i)
    return positions

def find_C(input_seq: str) -> list:

    positions = []
    for i, _ in enumerate(input_seq[:-1]):
        if input_seq[i] in 'Cc':
            positions.append(i)
    return positions

def find_non_CpG(input_seq: str) -> list:

    positions = []
    for i, _ in enumerate(input_seq[:-1]):
        if input_seq[i] in 'Cc' and input_seq[i+1] not in 'Gg':
            positions.append(i)
    return positions

def iter_print_meth(sample_path: Path, region_str: str, strand: bool, print_profiles: bool = False):

    REFERENCE_PATH: Path = Path("GCF_002263795.2_ARS-UCD1.3_genomic.fasta")

    direct_output = subprocess.run(['samtools', 'view', sample_path, region_str], capture_output=True)

    profiles = []
    all_positions = []
    methylation_per_read = []
    for line in direct_output.stdout.decode().split('\n'):
        if len(line.split('\t')) < 9: continue

        read_id, _, chrom, start_pos, _, _, _, _, _, seq, *_ = line.split('\t')

        ref_region_str: str = f"{chrom}:{start_pos}-{int(start_pos)+len(seq)-1}"

        ref_output = subprocess.run(['samtools', 'faidx', REFERENCE_PATH, ref_region_str], capture_output=True)
        ref_seq = ''.join([line for line in ref_output.stdout.decode().split('\n') if not line.startswith('>')])


        seq = seq if strand else reverse_complement(seq)
        ref_seq = ref_seq if strand else reverse_complement(ref_seq)

        C_positions = find_non_CpG(ref_seq)
        CpG_positions = find_CpG(ref_seq)

        C_pos_methylated = [1 if seq[pos] == 'C' else 0 for pos in C_positions]
        CpG_pos_methylated = [1 if seq[pos] == 'C' else 0 for pos in CpG_positions if seq[pos] in 'CT']
        if not CpG_pos_methylated: continue

        read_methylation_perc: float = ''.join([str(i) for i in CpG_pos_methylated]).count('1')/len(CpG_pos_methylated)


        mod = 1 if strand else -1
        adjusted_CpG_pos = [int(start_pos) + pos*mod for pos in CpG_positions]
        all_positions += adjusted_CpG_pos
        read_profile = {key: value for key, value in zip(adjusted_CpG_pos, CpG_pos_methylated)}        
        profiles.append(read_profile)

        if not print_profiles: 
            methylation_per_read.append(read_methylation_perc)
            #print(read_methylation_perc)

    if print_profiles:
        positions_to_exclude = Counter(all_positions).items()
        #print(Counter(all_positions).values(), 0.1*max(Counter(all_positions).values()))
        positions_to_exclude = list(filter(lambda x: x[1] < 0.1*max(Counter(all_positions).values()), positions_to_exclude))
        positions_to_exclude = [i[0] for i in positions_to_exclude]

        profiles_df = pd.DataFrame(profiles)
        profiles_translated = []
        for index, row in profiles_df.iterrows():
            genotype = ""
            for col_name in sorted(profiles_df.columns):
                if col_name in positions_to_exclude: continue
                i = row[col_name]
                if i == 1: genotype += "1"
                elif i == 0: genotype += "0"
                else: genotype += "-"
            profiles_translated.append(genotype)

        for key, value in sorted(Counter(profiles_translated).items(), key=lambda x: x[1]):
            print(key, value)
    
    return methylation_per_read

def main():
    #IMPRINTED_GENES: list = ("GNAS", "H19", "IGF2R", "KCNQ1", "MEST", "NNAT", "PEG10", "PEG3", "PLAGL1", "RTL1", "SNRPN", "XIST")
    input_path = list(Path('bam').glob(f'{sys.argv[1]}_*.bam'))[0]
    input_gene = sys.argv[2]
    sample_name = input_path.stem.split('_')[0]
    
    #samples_list = sorted(set([file.stem.split('_')[0] for file in Path('bam').glob('*.bam')]))

    methylation_per_read = iter_print_meth(input_path, *primers_dict[input_gene], False)

    output_path = Path('histo').joinpath(input_gene)
    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path.joinpath(f"{sample_name}.txt"), 'w') as output_file:
        for read in methylation_per_read:
            output_file.write(f"{read}\n")

if __name__=="__main__":
    main()
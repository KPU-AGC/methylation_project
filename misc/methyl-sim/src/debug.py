from Bio import Align
from Bio import SeqIO
import pickle
import pathlib
import pandas as pd

REF_POSITIONS: dict = {
    "DMAP1-BP1X": "NC_037330.1:101832510-101832748", 
    "DNMT3A-BP1X": "NC_037338.1:74014960-74015346", 
    "DNMT3B-BP2X": "NC_037340.1:62166550-62166873", 
    "EHMT2-BN1X": "NC_037350.1:27462619-27462871", 
    "GNAS-BN1X": "NC_037340.1:57520605-57520891", 
    "H19-BP1X": "NC_037356.1:49504607-49504946", 
    "IGF2R-BN2X": "NC_037336.1:96223067-96223367", 
    "KCNQ1-BP1X": "NC_037356.1:48907634-48907882", 
    "LIF-BN1X": "NC_037344.1:69264000-69264354", 
    "LIFR-BP1X": "NC_037347.1:35949053-35949420", 
    "MEST-BN1X": "NC_037331.1:94249893-94250283", 
    "NNAT-BP1X": "NC_037340.1:66465753-66466057", 
    "PEG10-BP1X": "NC_037331.1:12063279-12063673", 
    "PEG3-BN1X": "NC_037345.1:64120688-64120941", 
    "PLAGL1-BN1X": "NC_037336.1:81418732-81419041", 
    "RTL1-BN1X": "NC_037348.1:65778329-65778707", 
    "SLC2A8-BP1X": "NC_037338.1:98154376-98154714", 
    "SNRPN-BP1X": "NC_037348.1:1937284-1937513", 
    "SUV39H1-BN1X": "NC_037357.1:86781615-86781929", 
    "TXNIP-BP1X": "NC_037330.1:21383423-21383784", 
    "XIST-BP1X": "NC_037357.1:77162729-77163121", 
}

def _get_adjusted_pos(input_seq: str, start: int, end: int) -> dict:
    output_dict = {}
    running_count = 0
    for i, nuc in enumerate(input_seq):
        if nuc == '-': continue
        output_dict[i] = running_count+start
        running_count+=1
    return output_dict

with open(pathlib.Path('pkl/new_weights.pkl'), mode='rb') as input_file:
    data_df = pickle.load(input_file)

target_region = 'DNMT3A-BP1X'

refseq = data_df['templates'][target_region]
chromosome = REF_POSITIONS['DNMT3A-BP1X'].split(':')[0]
ref_start = int(REF_POSITIONS['DNMT3A-BP1X'].split(':')[1].split('-')[0])
ref_end = int(REF_POSITIONS['DNMT3A-BP1X'].split(':')[1].split('-')[1])


aligner = Align.PairwiseAligner(mode='local', match_score=2, mismatch_score=1)

methylated_CpG = []
unmethylated_CpG = []
for i in SeqIO.parse('SIMULATED.extendedFrags.fastq', 'fastq'):
    alignments = aligner.align(refseq, i.seq)
    adjusted_pos =_get_adjusted_pos(alignments[0][0], ref_start, ref_end)
    CpG_pos = [i for i, nuc in enumerate(alignments[0][0]) if nuc == 'C']
    meth_CpG_pos = [adjusted_pos[i] for i in CpG_pos if alignments[0][1][i] == 'C']
    unmeth_CpG_pos = [adjusted_pos[i] for i in CpG_pos if alignments[0][1][i] == 'T']

    methylated_CpG += meth_CpG_pos
    unmethylated_CpG += unmeth_CpG_pos

from collections import Counter

meth_pos = Counter(methylated_CpG)
uneth_pos = Counter(unmethylated_CpG)
for i in sorted(set(list(meth_pos.keys()) + list(uneth_pos.keys()))):
    val_1 = meth_pos.get(i, 0)
    val_2 = uneth_pos.get(i, 0)
    print(f"{chromosome}\t{i}\t{i}\t{val_1/(val_1+val_2)}\t{val_1}\t{val_2}")
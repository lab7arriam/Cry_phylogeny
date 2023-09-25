#/usr/bin/python3.7
import argparse
from collections import defaultdict
import os.path
import os
from Bio import pairwise2 as pw2
from Bio import SeqIO
from itertools import combinations
from annotate_clusters_consistency import write_csv

def find_sequence_identity(seq_tuple: tuple, ret_mode='ident', gap_mode='no_pen', get_alns='no') -> float:
    """
      Performes global pairwise alignmet of two sequences
      Returns blast-like identity between sequences
    """
    if gap_mode=='no_pen':
        global_align = pw2.align.globalxx(seq_tuple[0], seq_tuple[1])
    elif gap_mode=='local':
        global_align = pw2.align.localxs(seq_tuple[0], seq_tuple[1], -10.0, -0.5 )
    else:
        global_align = pw2.align.globalxs(seq_tuple[0], seq_tuple[1], -1, -0.1)
    
    seq_length = min(len(seq_tuple[0]), len(seq_tuple[1]))

    matches = global_align[0][2]
    percent_match_loc = (matches / seq_length) * 100

    al_length = global_align[0][4]
    max_seq_length = max(len(seq_tuple[0]), len(seq_tuple[1]))

    percent_match_glob = (matches / al_length) * 100
    percent_match_max_seq = (matches / max_seq_length) * 100
 
    if get_alns=='yes':
        return([global_align[0][0], global_align[0][1]])
    if ret_mode=='ident':
        return([percent_match_loc, percent_match_glob, al_length])
    elif ret_mode=='matches':
        return([round(matches), max_seq_length])
    elif ret_mode=='mathes_with_aln':
        return([round(matches), max_seq_length, al_length])

def make_sequence_dict(seq_dir: str) -> defaultdict(list):
    """
    Finds fasta files with domains and agregates them in a dictionary
    Returns a defaultidct(list) object with sequences of all 3 domains  
    """

    seq_dict=defaultdict(list)

    for filename in sorted(os.listdir(seq_dir)):
        #print(filename)
        if 'domain_1' in filename or 'domain1' in filename:
            for record in SeqIO.parse(os.path.join(seq_dir, filename),"fasta"):
                seq_dict[record.id]=[]
                seq_dict[record.id].append(record.seq.upper())

        elif 'domain_2' in filename or 'domain2' in filename:
            for record in SeqIO.parse(os.path.join(seq_dir, filename),"fasta"):
                seq_dict[record.id].append(record.seq.upper())

        elif 'domain_3' in filename or 'domain3' in filename:
            for record in SeqIO.parse(os.path.join(seq_dir, filename),"fasta"):
                seq_dict[record.id].append(record.seq.upper())

    return(seq_dict)

def agregate_identity(prot_dict:defaultdict(list), nucl_dict:defaultdict(list)) -> list:
    """
    Agregates domain's identity for nucleotide and protein sequences
    Returns a tsv-table with identity for each domain and delta between nucleotides and proteins  
    """
    ret_list_full=[['tox1','tox2','dom1_prot_loc','dom1_prot_glob','dom1_nucl_loc','dom1_nucl_glob','dom1_delta_loc','dom1_delta_glob'
                    ,'dom2_prot_loc','dom2_prot_glob','dom2_nucl_loc','dom2_nucl_glob','dom2_delta_loc','dom2_delta_glob',
                   'dom3_prot_loc','dom3_prot_glob','dom3_nucl_loc','dom3_nucl_glob','dom3_delta_loc','dom3_delta_glob']]
    tox_list=set()
    for pair_tup in list(combinations(list(prot_dict.keys()),2)):
        if pair_tup[0] not in tox_list:
            print(pair_tup[0])
            tox_list.add(pair_tup[0])
        add_list=[pair_tup[0], pair_tup[1]]
        for dom_ind in range(3):
            prot_pair=find_sequence_identity((prot_dict[pair_tup[0]][dom_ind],prot_dict[pair_tup[1]][dom_ind]))
            nucl_pair=find_sequence_identity((nucl_dict[pair_tup[0]][dom_ind],nucl_dict[pair_tup[1]][dom_ind]))
            add_list.extend(prot_pair)
            add_list.extend(nucl_pair)
            add_list.extend([abs(prot_pair[0]-nucl_pair[0]),abs(prot_pair[1]-nucl_pair[1])])
        if abs(prot_pair[1]-nucl_pair[1])>30:
            print(abs(prot_pair[1]-nucl_pair[1]))
        ret_list_full.append(add_list)

    return(ret_list_full)
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a table with identities of the domains for nucleotide and protein sequences')
    parser.add_argument('-o', '--out', dest='out_dir', help='the path to the output directory',
                        type=str)
    parser.add_argument('-p', '--prot', dest='prot_dir', help='the path to the proteins\'directory',
                        type=str)
    parser.add_argument('-n', '--nuc', dest='nuc_dir', help='the path to the nucleotides\'directory',
                        type=str)
    args = parser.parse_args()
    prot_dict = make_sequence_dict(args.prot_dir)
    nucl_dict = make_sequence_dict(args.nuc_dir)
    write_csv(agregate_identity(prot_dict, nucl_dict),os.path.join(os.path.realpath(args.out_dir),'pairs_identity_no_pat.tsv'))

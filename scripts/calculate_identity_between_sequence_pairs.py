import argparse
import csv
from collections import defaultdict
import os.path
from Bio import SeqIO, Entrez
from Bio import pairwise2 as pw2
from itertools import combinations, chain
from annotate_clusters_consistency import write_csv, create_dict_from_list, read_csv_to_list
from make_data_for_heatmap_universal import make_sequence_dict, find_sequence_identity
from itertools import combinations
import multiprocessing


def run_id_seqs(i) -> float:
    if i%1000==0:
        print(i)
    id1=seq_list[i][0][0]
    id2=seq_list[i][1][0]

    seq1=seq_list[i][0][1]
    seq2=seq_list[i][1][1]

    seq_ident = find_sequence_identity((seq1,seq2))[1]
    return((id1, seq_ident), (id2, seq_ident))


def find_mulitclust_similarity(fasta_file, out_file, threads) :
    """
        Parses fasta file and calculates mean identity between sequences
    """
    global seq_list
    seq_list=[]
    id_list=[]
    id_tuple_dict=defaultdict(list)

    for record in SeqIO.parse(fasta_file,"fasta"):
        seq_list.append((record.id, record.seq))

    seq_list=list(combinations(seq_list,2))

    pool = multiprocessing.Pool(threads)

    id_tuple_list = pool.map(run_id_seqs, range(0,len(seq_list)))
    for seq_tuple in id_tuple_list:
        if int(seq_tuple[0][1])==0 or int(seq_tuple[1][1])==0:
            pass
        else:
            id_tuple_dict[seq_tuple[0][0]].append(seq_tuple[0][1])
            id_tuple_dict[seq_tuple[1][0]].append(seq_tuple[1][1])

    for tox in id_tuple_dict:
        id_list.append([tox, sum(id_tuple_dict[tox])/len(id_tuple_dict[tox])])
        

    write_csv(id_list, out_file) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reads fasta-file and calculates identity between groups of sequences')
    parser.add_argument('-f', '--fast_f', dest='fasta_file', help='path to the fasta file',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_name', help='the name of the output file',
                        type=str)
    parser.add_argument('-t', '--th', dest='threads', help='the number of threads',
                        type=int, default=5)
    args = parser.parse_args()

    find_mulitclust_similarity(args.fasta_file, args.out_name, args.threads)

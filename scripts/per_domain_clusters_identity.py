
import argparse
import csv
from collections import defaultdict
import os.path
from Bio import SeqIO, Entrez
from Bio import pairwise2 as pw2
from itertools import combinations, chain
from annotate_clusters_consistency import write_csv, create_dict_from_list, read_csv_to_list
from make_data_for_heatmap_universal import make_sequence_dict, find_sequence_identity


def find_mulitclust_similarity(clust_dict : defaultdict(list), nucl_dict : defaultdict(list)) -> list:
    """
        Parses table with Cd-hit clusters and find unique sequences for each domains
        Returns a table with remained sequences and similarities
    """
    ret_rows=list()

    for clust in clust_dict:
        if clust!='NOV_Cry8Bb1_98.4_CS131028.1' and len(clust_dict[clust][0].split(';')) >1:
            for i in range(3):
                filt_inter=set()
                ref_seq=nucl_dict[clust][i]
                non_dups=dict()
                for rec in clust_dict[clust][0].split(';'):
                    if rec!=clust:
                        try:
                            id_thr=find_sequence_identity((ref_seq,nucl_dict[rec][i]))[0]
                            if id_thr<100.0:
                                non_dups[rec]=id_thr
                        except:
                            pass
                    
                filt_inter=set()
            
                if len(non_dups)>1:
                    comb_list=list(combinations(non_dups.keys(),2))
                    non_passed=set()
                    filt_inter.add(comb_list[0][0])
                    for pair in list(combinations(non_dups.keys(),2)):
                        th2=find_sequence_identity((nucl_dict[pair[0]][i],nucl_dict[pair[1]][i]))[0]
                        if pair[0] not in non_passed:
                            if th2==100.0:
                                non_passed.add(pair[1])
                            if pair[0] not in non_passed and th2!=100.0:
                                filt_inter.add(pair[0])

                if len(filt_inter)>0:
                    ret_rows.append([clust, 'domain' +str(i+1), ';'.join(filt_inter),';'.join([str(round(non_dups[tox],2)) for tox in filt_inter])])
                    print(clust, 'domain' +str(i+1),min([round(non_dups[tox],2) for tox in filt_inter]))

    return(ret_rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses CD-HT clusterins results and returns a table with unique per-domain sequences')
    parser.add_argument('-cf', '--cl_f', dest='clust_file', help='path to the CD-HIT cluster file',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_dir', help='the path to the output directory',
                        type=str)
    parser.add_argument('-n', '--nucl', dest='nucl_seqs', help='the path to the nucleotide sequences',
                        type=str)
    args = parser.parse_args()

    clust_dict = create_dict_from_list(read_csv_to_list(args.clust_file, headless=True,delim='\t'),0,9)
    nucl_dict = make_sequence_dict(args.nucl_seqs)
    
    write_csv(find_mulitclust_similarity(clust_dict,nucl_dict),os.path.join(os.path.realpath(args.out_dir),'unique_multiclusters.tsv'))

   

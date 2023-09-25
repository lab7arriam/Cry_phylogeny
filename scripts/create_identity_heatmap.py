from Bio import pairwise2 as pw2
from Bio import SeqIO
import sys
import csv
from collections import defaultdict
import argparse
import pandas as pd
import numpy as np
import os.path
from annotate_clusters_consistency import write_csv, read_csv_to_list
from build_recombination_graph import get_order_from_tree

def build_heatmap(tree_file, id_table, output, domain):
    cry_order=get_order_from_tree(tree_file)
    matrix=np.zeros((len(cry_order), len(cry_order)))
    heatmap_df=pd.DataFrame(matrix, columns=cry_order, index=cry_order)
    
    id_table = read_csv_to_list(os.path.realpath(id_table), headless=False)
    id_df = pd.DataFrame(id_table[1:], columns=id_table[0])
    #dom1_nucl_glob

    for tox1 in cry_order:
        print(tox1)
        for tox2 in cry_order:
            if tox1==tox2:
                heatmap_df.loc[tox1, tox2]=100
            else:
                row_heatmap=id_df.loc[(id_df['tox1'] == tox1) & (id_df['tox2'] == tox2)]
                if len(row_heatmap)==1:
                    pass
                else:
                    row_heatmap=id_df.loc[(id_df['tox2'] == tox1) & (id_df['tox1'] == tox2)]

                if domain=='domain3':
                    heatmap_df.loc[tox1, tox2]=float(row_heatmap['dom3_nucl_glob']) #dom3_nucl_loc
                elif domain=='domain2':
                    heatmap_df.loc[tox1, tox2]=float(row_heatmap['dom2_nucl_glob']) #dom2_nucl_loc
                elif domain=='domain1':
                    heatmap_df.loc[tox1, tox2]=float(row_heatmap['dom1_nucl_glob']) #dom1_nucl_loc

    heatmap_df.to_csv(output, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='builds identity heatmap based on a tree')
    parser.add_argument('-i', '--id', dest='id_tab', help='the path to the table with identity',
                        type=str)
    parser.add_argument('-t', '--tree', dest='tree_file', help='the path to the tree of the toxins',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_file', help='the path to the output directory',
                        type=str)
    parser.add_argument('-d', '--dom', dest='dom_type', help='speficied domein name',
                        type=str)
    args = parser.parse_args()

    build_heatmap(args.tree_file, args.id_tab, args.out_file, args.dom_type)

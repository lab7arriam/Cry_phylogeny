import os
import os.path
import click
import pandas as pd
from ete3 import Tree
from collections import defaultdict
from annotate_clusters_consistency import write_csv, read_csv_to_list
from get_breakpoints_from_alignments import  get_major_and_minor_parents
from Overall_toxicity_summary import get_parent_names
from make_data_for_heatmap_universal import make_sequence_dict, find_sequence_identity
from itertools import combinations
from statistics import mean



def make_clusters_stat(clust_tab, nucl_dict):
    cd_hit_clust_tab = read_csv_to_list(os.path.realpath(clust_tab), headless=True)
    clust_new_dict = defaultdict(list)

    for line in cd_hit_clust_tab:

        clust_list_uniq = line[4].split(';')
        clust_new_dict[line[0]]= clust_list_uniq

    for_clust_table = [['cluster_ref', 'num_tox', 'dom1_id','dom2_id','dom3_id']]

    for cluster in clust_new_dict:
        cluster_prots = clust_new_dict[cluster]
        cluster_combinations = list(combinations(cluster_prots, 2))
        if len(cluster_prots) == 1:
            ids = [100,100,100]

        else:
            ids = []

            for dom_ind in range(3):
                id_dom = []
                for pair_tup in cluster_combinations:
                    print(pair_tup)
                    pair_id = find_sequence_identity((nucl_dict[pair_tup[0]][dom_ind],nucl_dict[pair_tup[1]][dom_ind]))
                    id_dom.append(pair_id[0])
                ids.append(round(mean(id_dom),4))

        for_clust_table.append([cluster, len(cluster_prots)] + ids)

    
    write_csv(for_clust_table, 'tab_S2.csv')



            
@click.command()            
@click.option('--clust_tab', '-c', help="the file with clusters", 
              type=click.Path(exists=True), metavar='<PATH>')  
@click.option('--nucl_dir', '-n', help="the directory with domain sequences", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(clust_tab, nucl_dir): 

    nucl_dict = make_sequence_dict(nucl_dir)
    make_clusters_stat(clust_tab, nucl_dict)

if __name__ == '__main__':
   main() 

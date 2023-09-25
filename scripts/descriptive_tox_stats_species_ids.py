#/usr/bin/python3.7
import click
import os
from Bio import SeqIO
from collections import defaultdict
from annotate_clusters_consistency import read_csv_to_list, write_csv, create_dict_from_list
from make_data_for_heatmap_universal import make_sequence_dict, find_sequence_identity
from itertools import combinations, chain



def parse_annot_table(annot_table: str, tox_list: list) -> dict:
    """
        Parses annotation table and returns the list of accession numbers as well as missed toxins
    """

    ret_dict=dict()

    with open(annot_table,'r') as csv_file:
        my_reader = csv.reader(csv_file, delimiter='\t')
        next(my_reader)
        for row in my_reader:
            if row[0]!='--':
                tox=row[0]
                if tox in tox_list:
                    ret_dict[tox]=row[6].strip().replace(' ','').replace('\r','')

    return(ret_dict)


def summarize_species_for_toxins(annot_tab, clusters_all):
    annot_tab_stat=read_csv_to_list(os.path.realpath(annot_tab),delim='\t' ,headless=True)
    clusters_all_list=read_csv_to_list(os.path.realpath(clusters_all),delim='\t' ,headless=True)
    clusters_all_list=[el[0] for el in clusters_all_list]

    species_dict=defaultdict(list)
    for row in annot_tab_stat:
        if row[0]!='--':
            tox=row[0]
            try:
                if 'Bacillus thuringiensis' in row[8]:
                    species_dict[tox].append('Bacillus thuringiensis')
                else:
                    species_dict[tox].append(row[8])
            except:
                pass

    species_dict_ret=dict()

    for clust in clusters_all_list:
        if clust not in species_dict:
            print(clust)


    for tox in species_dict:
        if tox in clusters_all_list:
            if 'Bacillus thuringiensis' in species_dict[tox] or 'Unknown' in species_dict[tox] :
                species_dict_ret[tox]='Bacillus thuringiensis'
            elif 'Bacillus cereus' in species_dict[tox]:
                species_dict_ret[tox]='Bacillus cereus'
            elif 'Paenibacillus popilliae' in species_dict[tox]:
                species_dict_ret[tox]='Paenibacillus popilliae'
            elif 'Brevibacillus brevis' in species_dict[tox]:
                species_dict_ret[tox]='Brevibacillus brevis'
            elif 'Brevibacillus brevis X23' in species_dict[tox]:
                species_dict_ret[tox]='Brevibacillus brevis'
            elif 'Paenibacillus lentimorbus' in species_dict[tox]:
                species_dict_ret[tox]='Paenibacillus lentimorbus'
            elif 'Paraclostridium bifermentans' in species_dict[tox]:
                species_dict_ret[tox]='Paraclostridium bifermentans'
            elif 'Paenibacillus popilliae ATCC 14706' in species_dict[tox]:
                species_dict_ret[tox]='Paenibacillus popilliae'
            elif 'Bacillus mycoides KBAB4' in species_dict[tox]:
                species_dict_ret[tox]='Bacillus mycoides'
            else:
                species_dict_ret[tox]=species_dict[tox][0]            

    return(species_dict_ret)

    

@click.command()           
@click.option('--annot_tab', '-a', help="the table with IPG-base annotations", 
              type=click.Path(exists=True), metavar='<PATH>')  
@click.option('--nuc_dir', '-n', help="the path to the nucleotide sequences' directory", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--clusters_all', '-c', help="the table with names from all 733 clusters", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--prop_tab', '-p', help="the path to the file with domains' properties", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--clust_tab', '-l', help="the path to the file with clusterizations", 
              type=click.Path(exists=True), metavar='<PATH>')

def main(annot_tab, nuc_dir, prop_tab, clust_tab, clusters_all):
    nucl_dict = make_sequence_dict(nuc_dir)
    nov_id_dict=defaultdict(list)

    raw_prop_stat=read_csv_to_list(os.path.realpath(prop_tab),delim='\t' ,headless=True)
    properties_tab=[['Cry','Type','Domain','Length','ID_pair','ID_NOV','NOM_group', 'Species']]

    species_for_all_clusters = summarize_species_for_toxins(annot_tab, clusters_all)

    for tox in nucl_dict:
        if tox in species_for_all_clusters:
            if 'NOV' in tox:
                bt_ref='BT_'+tox.split('_')[1]
                print(bt_ref, tox)
                dom1_id=find_sequence_identity((nucl_dict[bt_ref][0],nucl_dict[tox][0]))[1]
                dom2_id=find_sequence_identity((nucl_dict[bt_ref][1],nucl_dict[tox][1]))[1]
                dom3_id=find_sequence_identity((nucl_dict[bt_ref][2],nucl_dict[tox][2]))[1]
                nov_id_dict[tox] = [dom1_id,dom2_id,dom3_id]

    for row in raw_prop_stat:
        NOV_flag=row[0].split('_')[0]
        if NOV_flag=='BT':
            seq_id=100
        else:
            if row[2]=='domain1':
                seq_id=nov_id_dict[row[0]][0]
            if row[2]=='domain2':
                seq_id=nov_id_dict[row[0]][1]
            if row[2]=='domain3':
                seq_id=nov_id_dict[row[0]][2]
        if row[0]=='BT_Cry10Aa4':
            species_for_all_clusters[row[0]]='Bacillus thuringiensis'
        properties_tab.append([row[0],row[3],row[2],row[5], row[1], seq_id, NOV_flag, species_for_all_clusters[row[0]]])

    write_csv(properties_tab,'mearged_domain_props.csv')

    raw_clust_stat=read_csv_to_list(os.path.realpath(clust_tab),delim='\t' ,headless=True)
    new_clust_tab=[['cluster', 'Num_tox','domain','ID']]

    clust_iter=1
    for cluster_row in raw_clust_stat:
        clust_iter += 1
        clust_seqs = cluster_row[9].split(';')

        comb_clusters=list(combinations(clust_seqs,2))

        dom1_id=100
        dom2_id=100
        dom3_id=100
        if len(clust_seqs)>1:
            clust_ids1=[]
            clust_ids2=[]
            clust_ids3=[]

            for pair in comb_clusters:
                dom1_id=find_sequence_identity((nucl_dict[pair[0]][0], nucl_dict[pair[1]][0]))[1]
                dom2_id=find_sequence_identity((nucl_dict[pair[0]][1], nucl_dict[pair[1]][1]))[1]
                dom3_id=find_sequence_identity((nucl_dict[pair[0]][2], nucl_dict[pair[1]][2]))[1]

                clust_ids1.append(dom1_id)
                clust_ids2.append(dom2_id)
                clust_ids3.append(dom3_id)

                dom1_id = sum(clust_ids1)/len(clust_ids1)
                dom2_id = sum(clust_ids2)/len(clust_ids2)
                dom3_id = sum(clust_ids3)/len(clust_ids3)
            

        new_clust_tab.append([clust_iter, cluster_row[8], 'domain1', dom1_id])
        new_clust_tab.append([clust_iter, cluster_row[8], 'domain2', dom2_id])
        new_clust_tab.append([clust_iter, cluster_row[8], 'domain3', dom3_id])
        print([clust_iter, cluster_row[8], 'domain1', dom1_id])

    write_csv(new_clust_tab,'Clusters_len_and_ids.csv')

if __name__ == '__main__':
   main()

import os
import click
import os.path
import pandas as pd
from statistics import mean
from collections import defaultdict
from ete3 import Tree
from create_evol_data_on_fixed_recombinants import gca, get_node_names
from annotate_clusters_consistency import write_csv, read_csv_to_list
from get_breakpoints_from_alignments import  get_major_and_minor_parents
from Overall_toxicity_summary import get_parent_names


#Check subtrees from the evol article:
#cry_list_from_article_full = ['BT_Cry1Aa10','BT_Cry1Ab2','BT_Cry1Ac39','BT_Cry1Ea10','BT_Cry1Ba1','BT_Cry9Ec1','BT_Cry9Ea10','BT_Cry9Eb1','BT_Cry9Da1','BT_Cry9Db1','BT_Cry9Ba1', 'BT_Cry8Aa1', 'BT_Cry8Da1', 'BT_Cry7Aa1','BT_Cry7Ab1']
#cry_list_from_article_no_cry7 = ['BT_Cry1Aa10','BT_Cry1Ab2','BT_Cry1Ac39','BT_Cry1Ea10','BT_Cry1Ba1','BT_Cry9Ec1','BT_Cry9Ea10','BT_Cry9Eb1','BT_Cry9Da1','BT_Cry9Db1','BT_Cry9Ba1', 'BT_Cry8Aa1', 'BT_Cry8Da1']


#for tree_ind, tree in zip([ 'dom1', 'dom2', 'dom3', 'full'], trees):
#    print(tree_ind, len(gca(tree, cry_list_from_article_no_cry7)), 'no_cry7', 'non-collapsed', sep ='\t')
#    print(tree_ind, len(gca(tree, cry_list_from_article_full)), 'cry7', 'non-collapsed', sep ='\t')

#testing_tree = gca(trees[0], ['BT_Cry29Aa1'])



def get_tree_list(tree_path):
    trees = [] #full, dom1, dom2, dom3
    for tf in sorted(os.listdir(os.path.realpath(tree_path))):
        tf_name = os.path.join(os.path.realpath(tree_path), tf)
        t = Tree(tf_name)
        trees.append(t)
    return(trees)


def get_node_depth(tree, label):
    path = []
    tree = tree.copy()

    D = tree&label
    while D.up:
        path.append(D)
        D = D.up
    return(len(path))


def find_the_longest_depth(tree):
    analyzed_tree = tree.copy()
    labels = get_node_names(analyzed_tree)

    depth_list = []
    for label in labels:
        depth_list.append(get_node_depth(analyzed_tree,label))
    return(max(depth_list))



def remove_uknonwn_from_tree(tree, exclude_list):
    analyzed_tree = tree.copy()


    for node in exclude_list:
        node_to_remove = analyzed_tree.search_nodes(name=node)[0]
        node_to_remove.detach()

    return(analyzed_tree)


def classify_nodes_in_subtree(node, tox_dict_events, parent_flag):
    toxins_in_node = get_node_names(node)

    events = []
    recs = set()
    pars = set()

    rec_node_flag = 0
    for toxin in toxins_in_node:
        for tox_type in ['rec']+parent_flag:
            if toxin in tox_dict_events[tox_type]:
                if tox_type=='rec':
                    recs.add(toxin)

                else:
                    pars.add(toxin)
                events.extend(tox_dict_events[tox_type][toxin])

    p_num = 0
    if recs.union(pars) == set(toxins_in_node):
        p_num = len(pars)
        rec_node_flag = 1
    return((toxins_in_node, rec_node_flag, p_num))


def analyze_tree_depth(tree, exclude_list, tox_dict_events, depth, parent_flag):

    tree_without_unknowns = remove_uknonwn_from_tree(tree, exclude_list)
    node_content = tree.get_cached_content(store_attr='name')
    
    passed_tox = list()

    nodes_with_rec = []
    nodes_without_recs = []
    num_pars_list = []

    for node in node_content:
        if len(node)>50:
            continue

        depth_node = find_the_longest_depth(node)


        if depth_node==depth:
            if  len(set(get_node_names(node)).intersection(set(passed_tox)))!=0:
                continue
            tox_list, rec_node_flag, num_pars = classify_nodes_in_subtree(node, tox_dict_events,parent_flag)

            if rec_node_flag==1:
                nodes_with_rec.append(node)
                num_pars_list.append(num_pars)
            else:
                nodes_without_recs.append(node)

            passed_tox.extend(tox_list)

    print(round(len(nodes_with_rec)/len(nodes_without_recs),3), round(mean(num_pars_list),3))

def analyze_unknowns(trees, rec_df):

    unknowns_dict={'domain1':[], 'domain2':[], 'domain3':[]}
    d_list = [1,2]

    tox_dict_events = defaultdict(dict)
    for domain in ['domain1', 'domain2', 'domain3']:
        tox_dict_events[domain] = defaultdict(dict)
        for tox_type in ['rec', 'parent_1', 'parent_2', 'parent_3']:
            tox_dict_events[domain][tox_type]= defaultdict(list)

    for line in rec_df.iterrows():
        
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']

        if rec_id in ['22','45','107']:
            continue

        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        tox_all = [el for el in major_all_pars + minor_pars+recs]

        dom1_pars =  get_parent_names(line[1]['parent_1'])
        dom2_pars =  get_parent_names(line[1]['parent_2'])
        dom3_pars =  get_parent_names(line[1]['parent_3'])

        per_dom_pars = dom1_pars+dom2_pars+dom3_pars

        unknown_flag = 0

        if 'unknown' in tox_all or 'unknown' in per_dom_pars:
            unknown_flag = 1
            unknowns_dict[rec_type].extend(tox_all)

        if unknown_flag==0:
            for rec in recs:
                if rec not in tox_dict_events[rec_type]['rec']:
                    tox_dict_events[rec_type]['rec'][rec] = [rec_id]
                else:
                    tox_dict_events[rec_type]['rec'][rec].append(rec_id)

            for par in dom1_pars:
                if par not in tox_dict_events[rec_type]['parent_1']:
                    tox_dict_events[rec_type]['parent_1'][par] = [rec_id]
                else:
                    tox_dict_events[rec_type]['parent_1'][par].append(rec_id)

            for par in dom2_pars:
                if par not in tox_dict_events[rec_type]['parent_2']:
                    tox_dict_events[rec_type]['parent_2'][par] = [rec_id]
                else:
                    tox_dict_events[rec_type]['parent_2'][par].append(rec_id)


            for par in dom3_pars:
                if par not in tox_dict_events[rec_type]['parent_3']:
                    tox_dict_events[rec_type]['parent_3'][par] = [rec_id]
                else:
                    tox_dict_events[rec_type]['parent_3'][par].append(rec_id)


               
    for dom in unknowns_dict:
        tox_list=list(set([el for el in unknowns_dict[dom] if el!='unknown']))
        unknowns_dict[dom]=tox_list

    #Minors:
    analyze_tree_depth(trees[0], unknowns_dict['domain1'], tox_dict_events['domain1'], 1, ['parent_1'])
    analyze_tree_depth(trees[1], unknowns_dict['domain2'], tox_dict_events['domain2'], 1, ['parent_2'])
    analyze_tree_depth(trees[2], unknowns_dict['domain3'], tox_dict_events['domain3'], 1, ['parent_3'])

    #Majors:
    analyze_tree_depth(trees[1], unknowns_dict['domain3'], tox_dict_events['domain3'], 1, ['parent_2'])

    print(round(1*1.2307692307692308/0.138297)) #minor D1
    ##print(round(1*1/0.00943)) #minor D2 
    print(round(3*1.3461538461/0.3333333333333333)) #minor D3
    print(round(5*1.04347826086956/0.2738095238095)) #major D3 + minor D2


#Minor(d1-d3)

#0.138 1.231
#0.009 1
#0.333 1.346

#Major(d2)

#0.274 1.043

#U(d):

#9
#106
#12
#15

#U(d) = n(d)*p(d)/f(d)


@click.command() 
@click.option('--rec_tab', '-r', help="the path with recombination results", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--tree_path', '-t', help="the path to renamed trees", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(rec_tab, tree_path): 
    file_reader = read_csv_to_list(os.path.realpath(rec_tab), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])
    
    trees = get_tree_list(tree_path)

    analyze_unknowns(trees, rec_df)


if __name__ == '__main__':
   main()





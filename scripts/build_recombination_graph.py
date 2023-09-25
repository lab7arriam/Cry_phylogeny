import sys
import csv
import os
from collections import defaultdict, Counter
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pygraphviz
import pandas as pd
import itertools
import argparse
import pandas as pd
from count_host_changes import get_parent_names
from annotate_clusters_consistency import write_csv, read_csv_to_list
from get_breakpoints_from_RDP import get_major_and_minor_parents

#Old color
# 1 - '#EEE8AA' - light yellow
# 2 - "#CD853F" - brown
# 3 -"#DC143C" - red
# 1 & 2 -"#556B2F" - green
# 1 & 3 - '#5F9EA0' - light blue
# 2 & 3 -'#DDA0DD' - pink

#New color
# 1 - #a60b0b
# 2  - #2980b9
# 3 - #bdbd00
# 1+2 - #556B2F
# 1+3 - #71C685
# 2+3 - #DDA0DD

#89D5D2 - both
#9A73BE - pars
#CD5922 - recs


def make_events(rec_filt):
    filt_res = read_csv_to_list(os.path.realpath(rec_filt), headless=False)
    filt_df = pd.DataFrame(filt_res[1:], columns=filt_res[0])

    filt_events_dict = defaultdict(dict)
    for line in filt_df.iterrows(): 
        ID = line[1]['ID']

        for index in ['ID', 'Rec','Min_par','Maj_par','Unknown_flag','Type','start', 'stop', 'd1_start', 'd1_stop', 'd2_start', 'd2_stop', 'd3_start', 'd3_stop','parent_1','parent_2','parent_3']:
            filt_events_dict[ID][index]= line[1][index]

        filt_events_dict[ID]['intersect_maj'] = get_major_and_minor_parents(line, mode='intersect')[2]
        filt_events_dict[ID]['new_all_maj'] = get_major_and_minor_parents(line, mode='intersect')[1]
        filt_events_dict[ID]['new_min_list'] = get_major_and_minor_parents(line, mode='intersect')[0]

    return(filt_events_dict)



def transform_nov_names(name):
    if 'NOV' in name:
        list_for_name = name.split('_')
        name = list_for_name[0]
        for i in range(1, len(list_for_name) - 1):
            name = name + '_' + list_for_name[i]
    return(name.replace('_WP',''))

def return_domain_names(domain: str) -> list:
    rtype = domain
    return ['domain1','domain2','domain3']


def make_multigraph(graph_dict, set_of_recs, set_of_parents):
    list_for_graph = []
    G = nx.MultiDiGraph()
    for node in graph_dict:
        for neigh in graph_dict[node]:
            for i in range(len(graph_dict[node][neigh])):
                col = graph_dict[node][neigh][i][0]
                connection_1 = [node, neigh, col]
                connection_2 = [neigh, node, col]
                if connection_1 not in list_for_graph and connection_2 not in list_for_graph:
                    if col == color_rec:
                        G.add_edge(node, neigh, color=col, dir='none') #, style='dotted'
                    else:
                        if node in set_of_parents and neigh in set_of_recs:
                            G.add_edge(node, neigh, color=col)
                        else:
                            G.add_edge(neigh, node, color=col)
                    list_for_graph.append(connection_1)
                    list_for_graph.append(connection_2)
    return(G)


def get_order_from_tree(tree):
    cry_order=[]
    with open (tree, 'r',newline='') as csvfile2:
        my_reader3 = csv.reader(csvfile2, delimiter='\t')
        for row in my_reader3:
            for cry in row[0].replace(')','').replace('(','').replace(';','').split(','):
                cry_order.append(cry.split(':')[0])
    return(cry_order)

def get_cry_rank(cry):
    if 'like' in cry:
        return([cry,cry])
    else:
        ret_ind=len(cry)
        third_rank=''

        for sybm in reversed(cry):
            if sybm.isnumeric():
                ret_ind-=1
            else:
                break

            third_rank=cry[:ret_ind]
   
        two_rank=''
        one_rank=''
        ret_ind=len(third_rank)

        for sybm in reversed(third_rank):
            if sybm.islower():
                ret_ind-=1
            else:
                break

            two_rank=third_rank[:ret_ind]
        

        ret_ind=len(two_rank)

        for sybm in reversed(two_rank):
            if not sybm.isnumeric():
                ret_ind-=1
            else:
                break

            one_rank=two_rank[:ret_ind]
        return([two_rank,one_rank])

   
def prepare_graph_dict(link_to_recombs, alpha_chanel='', rec_data = 'initial'):
    global rec_nodes
    global color_rec
    rec_nodes = '#CD5922' + alpha_chanel
    color_rec = '#000000' + alpha_chanel
    if rec_data =='initial':
        ind_dict = {'domain1':51, 'domain2':52, 'domain3':53} 
    else:
        ind_dict = {'domain1':21, 'domain2':22, 'domain3':23} #for combined graph



    global dict_of_NOV_names
    global set_of_recs
    global set_of_parents

    set_of_recs = set()
    set_of_parents = set()

    dict_of_NOV_names=dict()

    with open(link_to_recombs, 'r', newline='') as recombs:
        recombs_reader = csv.reader(recombs, delimiter='\t')
        next(recombs_reader)
        graph_dict = {}
        node_color_dict = {}
        for line in recombs_reader:
            id = line[0]
            recs = line[1].replace(' ', '').split(';')
            domain = line[5]

            pars1 = sorted(list(set(get_parent_names(line[ind_dict[return_domain_names(domain)[0]]]))))
            pars2 = sorted(list(set(get_parent_names(line[ind_dict[return_domain_names(domain)[1]]]))))
            pars3 = sorted(list(set(get_parent_names(line[ind_dict[return_domain_names(domain)[2]]]))))
            pars12=[]
            pars23=[]
            pars13=[]	

            for name in recs+pars1+pars2+pars3:
                dict_of_NOV_names[transform_nov_names(name)]=name

            recs=[transform_nov_names(name) for name in recs]
            pars1=[transform_nov_names(name) for name in pars1]
            pars2=[transform_nov_names(name) for name in pars2]
            pars3=[transform_nov_names(name) for name in pars3]
            
            col1='#a60b0b'+alpha_chanel
            col2='#2980b9'+alpha_chanel
            col3='#bdbd00'+alpha_chanel
            col12='#556B2F'+alpha_chanel
            col23='#71C685'+alpha_chanel
            col13='#DDA0DD'+alpha_chanel
            global col_to_links_dict
            col_to_links_dict={col1:'domain1', col2:'domain2', col3:'domain3',col12:'domains12', col23:'domain23',col13:'domains13', color_rec:'recs'}
            

            for par1 in pars1.copy():
                for par2 in pars2.copy():
                    if par1==par2:
                        pars12.append(par1)
                        pars1.remove(par1)
                        pars2.remove(par2)

            for par2 in pars2.copy():
                for par3 in pars3.copy():
                    if par2==par3:
                        pars23.append(par2)
                        pars2.remove(par2)
                        pars3.remove(par3)

            for par1 in pars1.copy():
                for par3 in pars3.copy():
                    if par1==par3:
                        pars13.append(par1)
                        pars1.remove(par1)
                        pars3.remove(par3)

            for rec in recs:
                if rec not in graph_dict.keys():
                    graph_dict[rec] = {}
                    node_color_dict[rec] = rec_nodes
                set_of_recs.add(rec)
                for rec1 in recs:
                    if rec1 != rec:
                        if rec1 not in graph_dict[rec].keys():
                            graph_dict[rec][rec1] = []
                        if [color_rec, id] not in graph_dict[rec][rec1]:
                            graph_dict[rec][rec1].append([color_rec, id])
                
                for par1 in pars1:
                    if par1 != 'unknown':
                        set_of_parents.add(par1)
                        if par1 not in graph_dict[rec].keys():
                            graph_dict[rec][par1] = []
                        if [col1, id] not in graph_dict[rec][par1]:
                            graph_dict[rec][par1].append([col1, id])

                for par2 in pars2:
                    if par2 != 'unknown':
                        set_of_parents.add(par2)
                        if par2 not in graph_dict[rec].keys():
                            graph_dict[rec][par2] = []
                        if [col2, id] not in graph_dict[rec][par2]:
                            graph_dict[rec][par2].append([col2, id])

                for par3 in pars3:
                    if par3 != 'unknown':
                        set_of_parents.add(par3)
                        if par3 not in graph_dict[rec].keys():
                            graph_dict[rec][par3] = []
                        if [col3, id] not in graph_dict[rec][par3]:
                            graph_dict[rec][par3].append([col3, id])

                for par12 in pars12:
                    if par12 != 'unknown':
                        set_of_parents.add(par12)
                        if par12 not in graph_dict[rec].keys():
                            graph_dict[rec][par12] = []
                        if [col12, id] not in graph_dict[rec][par12]:
                            graph_dict[rec][par12].append([col12, id])
                #if id=='2':
                #    print('12',graph_dict['NOV_Cry1Da1_84.8'])

                for par23 in pars23:
                    if par23 != 'unknown':
                        set_of_parents.add(par23)
                        if par23 not in graph_dict[rec].keys():
                            graph_dict[rec][par23] = []
                        if [col23, id] not in graph_dict[rec][par23]:
                            graph_dict[rec][par23].append([col23, id])

                for par13 in pars13:
                    if par13 != 'unknown':
                        set_of_parents.add(par13)
                        if par13 not in graph_dict[rec].keys():
                            graph_dict[rec][par13] = []
                        if [col13, id] not in graph_dict[rec][par13]:
                            graph_dict[rec][par13].append([col13, id])


    global mix_nodes
    global par_nodes
    mix_nodes = '#89D5D2'+alpha_chanel
    par_nodes = '#9A73BE'+alpha_chanel


    global col_to_nodes_dict
    col_to_nodes_dict={mix_nodes:'mixed', par_nodes:'par', rec_nodes:'rec'}
            

    
    global mixed_set
    mixed_set=set()
    for par in set_of_parents:
        if par in set_of_recs:
            mixed_set.add(par)
            node_color_dict[par] = mix_nodes
        else:
            node_color_dict[par] = par_nodes


    G=make_multigraph(graph_dict, set_of_recs, set_of_parents)

    return(G, graph_dict, node_color_dict)

def draw_graph_nx(G, node_color_dict, out_name, size_tuple):
    new_G =G.copy()

    colors_nodes = list()
    for i in new_G.nodes():
        node_color = node_color_dict[i]
        colors_nodes.append(node_color)

    colors_edges = nx.get_edge_attributes(new_G, 'color').values()

    pos = nx.spring_layout(new_G, k=0.8, iterations=100, scale=0.25) #, scale=0.25

    fig = plt.figure(figsize=size_tuple)

    nx.draw(new_G, pos, node_color=colors_nodes, edge_color=colors_edges, with_labels=True, node_size = 2100, 
             font_size=4.5, fontname='times-bold', width=2.5, 
             style='filled', shape='circle', connectionstyle='arc3, rad = 0.1')
    plt.savefig(out_name)




def print_graph_karate(G, node_color_dict, layout, out_name):

    new_G = G.copy()
    
    karate_agr = nx.nx_agraph.to_agraph(new_G)
    karate_agr.node_attr.update(style='filled', shape='circle', gradientangle=90, width=3.5, fixedsize=True, penwidth=6, fontsize=19, fontname='times-bold', bgcolor="transparent") 
    karate_agr.edge_attr.update(penwidth=11.0, len=6.0, style='filled')

    karate_agr.edge_attr.update(dir='forward', penwidth=11.0, len=3.0, arrowhead='normal', arrowtail='normal')
    for i in karate_agr.nodes():
        n = karate_agr.get_node(i)
        if n not in node_color_dict.keys():
            node_color_dict[n] = 'yellow'
            print('NO')
        node_color = node_color_dict[n]
        n.attr['fillcolor'] = node_color

    karate_agr.draw(out_name, prog=layout)


def make_circos(tree, graph_dict):

    cry_order=get_order_from_tree(tree)
    cry_order_codes={cry_order[i]:i for i in range(len(cry_order))}

    data_for_circus_sectors=[]
    data_for_circus_sectors.append(['Cry','Sector','Color','two_rank','one_rank'])

    links=[]
    links.append(['Cry','Sector_from','Sector_to','Color'])


    for cry in cry_order:

        
        if transform_nov_names(cry) in mixed_set:

            data_for_circus_sectors.append([cry,cry_order_codes[cry],mix_nodes]+get_cry_rank(cry.split('_')[1]))

            for edge in graph_dict[transform_nov_names(cry)]:
                try:
                    links.append([cry, cry_order_codes[edge], cry_order_codes[cry],graph_dict[transform_nov_names(cry)][edge][0][0]])
                except:
                    links.append([cry, cry_order_codes[dict_of_NOV_names[edge]], cry_order_codes[cry],graph_dict[transform_nov_names(cry)][edge][0][0]])

        elif transform_nov_names(cry) in set_of_recs:
            for edge in graph_dict[transform_nov_names(cry)]:
                try:
                    links.append([cry, cry_order_codes[edge], cry_order_codes[cry],graph_dict[transform_nov_names(cry)][edge][0][0]])
                except:
                    links.append([cry, cry_order_codes[dict_of_NOV_names[edge]], cry_order_codes[cry],graph_dict[transform_nov_names(cry)][edge][0][0]])

            data_for_circus_sectors.append([cry,cry_order_codes[cry],rec_nodes]+get_cry_rank(cry.split('_')[1]))

        elif transform_nov_names(cry) in set_of_parents:
            data_for_circus_sectors.append([cry,cry_order_codes[cry],par_nodes]+get_cry_rank(cry.split('_')[1]))

        else:
            data_for_circus_sectors.append([cry,cry_order_codes[cry],'grey']+get_cry_rank(cry.split('_')[1]))

    write_csv(data_for_circus_sectors,'test_base_for_circus_plot.tsv')
    write_csv(links,'test_links_for_circus_plot.tsv')
    


def dict_vals_to_str(parse_dict):
    valuses_list=[]
    for key in parse_dict:
        valuses_list.append(str(key)+':'+str(parse_dict[key]))
    return('; '.join(valuses_list))

def res_str_to_list(Comp_id,res_str, res_type):
    ret_list = []
    for res_substr in res_str.split('; '):
        ret_list.append([Comp_id, res_substr.split(':')[0], res_substr.split(':')[1], res_type])
    return(ret_list)
    

def analyze_components(G, graph_dict, node_color_dict, filt_events_dict):
    i = 1
    print('Number of edges:', len(G.edges()))
    print('Number of nodes:', len(G.nodes()))
    comps = [G.subgraph(c).copy() for c in nx.weakly_connected_components(G)]

    wide_res_list = [['Component', 'ID_num', 'IDs', 'Event_types', 'Edge_num', 'Edges', 'Node_num','Nodes']]
    long_res_list = [['Component', 'Value', 'Num', 'Type']]

    for comp in comps:
        list_of_ids = []
        comp_copy = comp.copy()
        for node1 in comp_copy:
            for node2 in comp_copy:
                if node1 in graph_dict.keys() and node2 in graph_dict[node1].keys():
                    ids = [item[1] for item in graph_dict[node1][node2]]
                    for id in ids:
                        if id not in list_of_ids:
                            list_of_ids.append(id)

        edges = comp_copy.edges()
        nodes = comp_copy.nodes()

        egge_num = len(edges)
        node_num = len(nodes)
        ids_num = len(list_of_ids)
        edge_types = []
        
        for edge_tup in edges:
            if edge_tup[0] in graph_dict.keys() and edge_tup[1] in graph_dict[edge_tup[0]].keys():
                edge_cols = graph_dict[edge_tup[0]][edge_tup[1]]
            elif edge_tup[1] in graph_dict.keys() and edge_tup[0] in graph_dict[edge_tup[1]].keys():
                edge_cols = graph_dict[edge_tup[1]][edge_tup[0]]
            for edge in edge_cols:
                link_type = col_to_links_dict[edge[0]]
                edge_types.append(link_type)

        edge_types = Counter(edge_types)
        event_types = Counter([filt_events_dict[ID]['Type'] for ID in list_of_ids])
        node_types = Counter([col_to_nodes_dict[node_color_dict[node]] for node in nodes])

        wide_row = [i, ids_num, '; '.join(list_of_ids),dict_vals_to_str(event_types) ,egge_num,dict_vals_to_str(edge_types), node_num, dict_vals_to_str(node_types)]
        wide_res_list.append(wide_row )
        i += 1

    inds_dict = {3:'Event_type', 5:'Link_type', 7:'Node_type'}
    for res_row in wide_res_list[1:]:
        Comp_id = res_row[0]
        for ID in res_row[2].split('; '):
            long_res_list.append([Comp_id, ID, 1, 'Event_ID'])
        for ind in inds_dict:
            long_res_list.extend(res_str_to_list(Comp_id, res_row[ind], inds_dict[ind]))
        
    write_csv(wide_res_list,'Rec_components_stat_wide.csv')
    write_csv(long_res_list,'Rec_components_stat_long.csv')


def longest_simple_paths(graph, source, target):
    longest_paths = []
    longest_path_length = 0
    for path in nx.all_simple_paths(G, source=source, target=target):
        if len(path) > longest_path_length:
            longest_path_length = len(path)
            longest_paths.clear()
            longest_paths.append(path)
        elif len(path) == longest_path_length:
            longest_paths.append(path)
    return(longest_paths)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A script for generating graph corresponding to strains')
    parser.add_argument('-r', '--rec', dest='rec_tab', help='the path to the preparatory table for selection evaluation',
                        type=str)
    parser.add_argument('-t', '--tree', dest='tree_file', help='the path to the tree of the toxins',
                        type=str)
    args = parser.parse_args()

    G, graph_dict, node_color_dict  = prepare_graph_dict(args.rec_tab)
    filt_events_dict = make_events(args.rec_tab)

    S = [G.subgraph(c).copy() for c in nx.weakly_connected_components(G)]
    print(S[0].number_of_nodes(), S[0].nodes())

    make_circos(tree, graph_dict)

    #longest_paths = longest_simple_paths(S[0], source=0, target=3)
    #if longest_paths:
    #    print(f"Longest simple path contains {len(longest_paths[0])} nodes")
    #    print(longest_paths)

    #for node_pair in list(itertools.combinations(S[0].nodes(), 2)):
    #    node1, node2 = node_pair
    #    longest_paths = longest_simple_paths(S[0], source=node1, target=node2)
    #    if longest_paths:
    #        print(f"Longest simple path contains {len(longest_paths[0])} nodes")
    #        print(longest_paths)


    #Draw plots:
    #print_graph_karate(G, node_color_dict, 'dot', 'All_comp_dot.pdf')
    for i in range(1,len(S )+1):
        #Comp_neato_no_transp Comp_neato_90 Comp_dot_no_transp Comp_dot_90
        out_name = str(i)+'_Comp_dot_90.pdf'
        Comp = S[i-1].copy()
        print_graph_karate(Comp, node_color_dict, 'dot', out_name)
       

    analyze_components(G, graph_dict, node_color_dict, filt_events_dict)





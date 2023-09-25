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
from count_host_changes import get_parent_names
from collections import defaultdict
from annotate_clusters_consistency import write_csv, read_csv_to_list
from overall_toxicity_summary import make_tox_summary_dict, clean_tox_dict, make_strains_dictionary, reduce_dicts_to_third_rank, make_extended_stain_tox_dict, crop_cry_name
from get_breakpoints_from_alignments import  get_major_and_minor_parents, unpacking_list_of_lists
from cry_combinations_per_strain_and_assembly import get_jaccard_coeff, get_simpson_coeff

def prepare_dict_for_graph(tox_to_strain_dict, strain_to_tox_dict):
    
    graph_dict_with_tox = {}

    for strain in strain_to_tox_dict:
        if strain not in graph_dict_with_tox:
            graph_dict_with_tox[strain] = {}

        for strain_to_compare in strain_to_tox_dict:
            if strain != strain_to_compare:
                tox_intersections = list(set(strain_to_tox_dict[strain]).intersection(set(strain_to_tox_dict[strain_to_compare])))
                if len(tox_intersections) >0: 
                    graph_dict_with_tox[strain][strain_to_compare] = [len(tox_intersections), tox_intersections]


    for strain in strain_to_tox_dict:
        if len(graph_dict_with_tox[strain])<1:
            graph_dict_with_tox[strain][strain]= [len(strain_to_tox_dict[strain]),strain_to_tox_dict[strain]]

    return(graph_dict_with_tox)
    

def prepare_color_dict(graph_dict_with_tox, extended_tox_to_res_dict, spec_ord_dict):
    color_dict = {}

    hosts_to_cols_dict = {'Lepidoptera':'#00468BFF',
    'Diptera':'#ED0000FF',
    'Coleoptera':'#42B540FF',
    'Nematoda':'#0099B4FF',
    'Hemiptera':'#925E9FFF',
    'Human':'#FDAF91FF',
    'Hymenoptera':'#AD002AFF',
    'Orthoptera':'#ADB6B6FF',
    'grey':'#CFD6D8'}  #708090 
    host_orders_dict_for_nodes = defaultdict(dict)    

    for node in graph_dict_with_tox:
        host_orders_dict_for_nodes[node] = defaultdict(list)
        host_orders_dict_for_nodes[node]['Orders'] = []
        host_orders_dict_for_nodes[node]['Species'] = []

        if node not in extended_tox_to_res_dict['Species']:
            color_dict[node] = [hosts_to_cols_dict['grey']]
        else:
            colors = []
            hosts = []
            hosts_species = []

            for species in extended_tox_to_res_dict['Species'][node]:
                 
                if species != 'Unknown':
                    order = spec_ord_dict[species]
                    hosts_species.append(species)
                else:
                    order = extended_tox_to_res_dict['Orders'][node][0] 


                if order in hosts_to_cols_dict:
                    colors.append(hosts_to_cols_dict[order])
                    hosts.append(order)

            color_dict[node] = colors
            host_orders_dict_for_nodes[node]['Orders'] = sorted(hosts)
            host_orders_dict_for_nodes[node]['Species'] = sorted(list(set(hosts_species)))

           # if len(colors) > 4:
           #     print(node, hosts)


    return((color_dict, host_orders_dict_for_nodes))

def build_strains_graph(graph_dict):

    G = nx.Graph()
    for node in graph_dict:
        for neigh in graph_dict[node]:
            G.add_edge(node, neigh, weight=graph_dict[node][neigh][0])


    print('Number of edges:', G.number_of_edges())
    print('Number of nodes:', G.number_of_nodes())
    print('Number of connected components:', nx.number_connected_components(G))

    return(G)

def draw_graph(new_G, color_dict, out_name):

    #, iterations = 100
    pos = nx.spring_layout(new_G, scale=0.75) #scale=0.55, 1.75  spring_layout scale=4.95

    fig = plt.figure(figsize=(13, 13)) # 25 25 10 10
    ax = plt.axes([0, 0, 1, 1])
    ax.set_aspect('equal')
    weights_edges = [new_G[u][v]['weight']*0.25 for u, v in new_G.edges()] #0.3 0.7 0.9 
    nx.draw_networkx_edges(new_G, pos, ax=ax, width=weights_edges, alpha=0.4) #0.5 0.6

    plt.xlim(-1.3, 1.3) #-1.5 1.5
    plt.ylim(-1.3, 1.3)

    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform

    piesize = 0.05 #0.025  0.09
    p2 = piesize/2.0
    for n in new_G:
        xx, yy = trans(pos[n]) 
        xa, ya = trans2((xx, yy)) 	
        a = plt.axes([xa-p2, ya-p2, piesize, piesize])
        a.set_aspect('equal')
        colors = color_dict[n]

        num_of_parts = len(set(colors))
        fracs = []

        color_counter = Counter(colors)
        color_set = sorted(list(set(colors)))
        color_all_num = len(colors)

        for color in color_set:
            color_fraq = color_counter[color] / color_all_num

            fracs.append(100. * color_fraq)

        a.pie(fracs, colors=color_set, wedgeprops={"edgecolor":"k",'linewidth': 1, 'linestyle': 'solid', 'antialiased': True, 'alpha':0.75})

    plt.savefig(out_name)


def get_components_from_tox_list(tox_list, components_dict, tox_to_strain_dict):

    components_results = dict()

    
    for tox in tox_list:

        strains = list(set(tox_to_strain_dict[crop_cry_name(tox)]))
        
        components = list(set([components_dict[strain] for strain in strains]))
        if len(components) >0:
            components_results[tox] = '; '.join([str(el) for el in components])

    if len(components_results) >0:
        return(components_results)
    else:
        return(None)


def get_str_with_comps(comp_res, tox_list):
    if not(comp_res):
        return(None)

    else:
        tox_dat = '; '.join([tox+'('+comp_res[tox]+')' for tox in tox_list if tox in comp_res])
        comps_int = sorted(list(set([int(comp_res[tox]) for tox in tox_list if tox in comp_res])))
        comps = '; '.join([str(comp) for comp in comps_int])

        return((comps_int, comps, tox_dat))


def analyze_connected_components(G, color_dict, rec_df, tox_to_strain_dict, strain_to_tox_dict, graph_dict, host_orders_dict_for_nodes):

    components_dict = dict()

    comp_iter = 0
    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]


    components_to_nodes_dict = defaultdict(list)
    for s in S:
        node_list = s.nodes()
        components_to_nodes_dict[comp_iter] = node_list
        for node in node_list:
            components_dict[node]=comp_iter

        comp_iter+=1

    rec_components_results = [['ID', 'Domain', 'Rec_comps','Rec_tox', 'Min_comps','Min_tox', 'Min_jac','Maj_comps','Maj_tox','Maj_jac', 'Inter_par_jac']]
    rec_jaccard_long = [['ID', 'Domain', 'Coeff', 'Type']]

    comp_to_ev_dict = defaultdict(list)
    tox_types_dict = defaultdict(dict)
    tox_comps_simps = defaultdict(dict)
    rec_to_dom_dict = dict()

    for line in rec_df.iterrows():
       
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']

        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        
        rec_comps = get_components_from_tox_list(recs, components_dict, tox_to_strain_dict)
        min_comps = get_components_from_tox_list(minor_pars, components_dict, tox_to_strain_dict)
        maj_comps = get_components_from_tox_list(major_inter_pars, components_dict, tox_to_strain_dict)

        if (not min_comps and not maj_comps):
            continue

        if not rec_comps:
            continue

        rec_to_dom_dict[rec_id] = rec_type
        rec_comps_int, rec_comps_str, rec_tox_comps = get_str_with_comps(rec_comps, recs)

        for comp in rec_comps_int:
            if comp not in comp_to_ev_dict:
                comp_to_ev_dict[comp] = [rec_id]
            else:
                comp_to_ev_dict[comp].append(rec_id)
        
        for rec in recs:
            tox_types_dict[rec_id]['Recs'] = [crop_cry_name(tox) for tox in recs]

        if min_comps and maj_comps:

            tox_types_dict[rec_id]['Mins'] = [crop_cry_name(tox) for tox in minor_pars]
            tox_types_dict[rec_id]['Majs'] = [crop_cry_name(tox) for tox in major_inter_pars]

            min_comps_int, min_comps_str, min_tox_comps = get_str_with_comps(min_comps, minor_pars)
            maj_comps_int, maj_comps_str, maj_tox_comps = get_str_with_comps(maj_comps, major_inter_pars)

            min_simp = get_jaccard_coeff(rec_comps_int,min_comps_int)
            maj_simp = get_jaccard_coeff(rec_comps_int,maj_comps_int)
            inter_simp = get_jaccard_coeff(min_comps_int,maj_comps_int)
            all_sipm = get_jaccard_coeff(rec_comps_int, min_comps_int+maj_comps_int)

            rec_jaccard_long.append([rec_id, rec_type, min_simp, 'rec_min'])
            rec_jaccard_long.append([rec_id, rec_type, maj_simp, 'rec_maj'])
            rec_jaccard_long.append([rec_id, rec_type, all_sipm, 'rec_all'])
            rec_jaccard_long.append([rec_id, rec_type, min_simp, 'par_inter'])

            tox_comps_simps[rec_id]['Mins'] = min_simp
            tox_comps_simps[rec_id]['Maj'] = maj_simp

            rec_components_results.append([rec_id, rec_type, rec_comps_str, rec_tox_comps, min_comps_str, min_tox_comps, min_simp, maj_comps_str, maj_tox_comps, maj_simp, inter_simp])

            for comp in min_comps_int:
                if comp not in comp_to_ev_dict:
                    comp_to_ev_dict[comp] = [rec_id]
                else:
                    comp_to_ev_dict[comp].append(rec_id)


            for comp in maj_comps_int:
                if comp not in comp_to_ev_dict:
                    comp_to_ev_dict[comp] = [rec_id]
                else:
                    comp_to_ev_dict[comp].append(rec_id)


        elif min_comps:
            min_comps_int, min_comps_str, min_tox_comps = get_str_with_comps(min_comps, minor_pars)
            min_simp = get_jaccard_coeff(rec_comps_int,min_comps_int)
            rec_components_results.append([rec_id, rec_type, rec_comps_str, rec_tox_comps, min_comps_str, min_tox_comps, min_simp, '', '', '', ''])
            rec_jaccard_long.append([rec_id, rec_type, min_simp, 'rec_min'])


            tox_comps_simps[rec_id]['Mins'] = min_simp

            tox_types_dict[rec_id]['Mins'] = [crop_cry_name(tox) for tox in minor_pars]

            for comp in min_comps_int:
                if comp not in comp_to_ev_dict:
                    comp_to_ev_dict[comp] = [rec_id]
                else:
                    comp_to_ev_dict[comp].append(rec_id)

        elif maj_comps:
            maj_comps_int, maj_comps_str, maj_tox_comps = get_str_with_comps(maj_comps, major_inter_pars)
            maj_simp = get_jaccard_coeff(rec_comps_int,maj_comps_int)
            rec_components_results.append([rec_id, rec_type, rec_comps_str, rec_tox_comps,  '', '', '',maj_comps_str, maj_tox_comps, maj_simp,  ''])
            rec_jaccard_long.append([rec_id, rec_type, maj_simp, 'rec_maj'])


            tox_types_dict[rec_id]['Majs'] = [crop_cry_name(tox) for tox in major_inter_pars]

            tox_comps_simps[rec_id]['Maj'] = maj_simp

            for comp in maj_comps_int:
                if comp not in comp_to_ev_dict:
                    comp_to_ev_dict[comp] = [rec_id]
                else:
                    comp_to_ev_dict[comp].append(rec_id)


    print(components_to_nodes_dict[57], strain_to_tox_dict['FCC41'], comp_to_ev_dict[57], host_orders_dict_for_nodes['FCC41']['Orders'])
    print(components_to_nodes_dict[132], strain_to_tox_dict['Bun1-14'],  comp_to_ev_dict[132], host_orders_dict_for_nodes['Bun1-14']['Orders'])

    components_hosts_res = [['Component', 'Host_type', 'Num_sp', 'Comp_type', 'Num_strains','Coeff', 'Num_hosts', 'Num_events', 'Num_tox_types']]


    hosts_coeff_dict = {'Lepidoptera':10000,
    'Diptera':7000,
    'Coleoptera':4000,
    'Nematoda':2500,
    'Hemiptera':1200,
    'Human':700,
    'Hymenoptera':600,
    'Orthoptera':500,
    'Unknown':0}  #708090 


    for comp in  components_to_nodes_dict:

        if comp in comp_to_ev_dict:
            events = sorted(list(set(comp_to_ev_dict[comp])))
        else:
            events = []

        num_tox_types = len(sorted(list(set(unpacking_list_of_lists([strain_to_tox_dict[strain] for strain in components_to_nodes_dict[comp]])))))
         
        comp_hosts_ords = []
        for node in components_to_nodes_dict[comp]:
            comp_hosts_ords.extend(host_orders_dict_for_nodes[node]['Orders'])

        comp_hosts_ords_count = Counter(comp_hosts_ords)
        num_hosts = len(comp_hosts_ords_count)

        if len(comp_hosts_ords_count)==0:
            comp_hosts_ords_count['Unknown']=1

        coeff = sum([hosts_coeff_dict[host] for host in comp_hosts_ords_count])

        for host in comp_hosts_ords_count:
            components_hosts_res.append([comp, host, comp_hosts_ords_count[host], 'All', len(components_to_nodes_dict[comp]), coeff+len(components_to_nodes_dict[comp]), num_hosts,len(set(events)), num_tox_types])
        
        
        if comp in comp_to_ev_dict:
            events = sorted(list(set(comp_to_ev_dict[comp])))
            for host in comp_hosts_ords_count:
                components_hosts_res.append([comp, host, comp_hosts_ords_count[host], 'Rec', len(components_to_nodes_dict[comp]), coeff+len(components_to_nodes_dict[comp]), num_hosts,len(set(events)), num_tox_types])
        


    def get_num_hosts_for_comp(comp, components_to_nodes_dict, host_orders_dict_for_nodes):
        comp_hosts_ords = []
        for node in components_to_nodes_dict[int(comp)]:
            comp_hosts_ords.extend(host_orders_dict_for_nodes[node]['Orders'])

        comp_hosts_ords_count = Counter(comp_hosts_ords)

        return(list(comp_hosts_ords_count.keys()))

    def make_comp_summary_for_tox(event, tox_comps,components_to_nodes_dict,host_orders_dict_for_nodes, tox_type):
        summary_comparision = []
        for tox in tox_comps:
            comps = tox_comps[tox]
            comp_hosts = []
            summary_comparision.append([event,tox, tox_type, comps,  len(components_to_nodes_dict[int(comps)]), '; '.join(get_num_hosts_for_comp(comps, components_to_nodes_dict, host_orders_dict_for_nodes))]) 
        return(summary_comparision)


    hosts_ind_dict = {'Lepidoptera':3,
    'Diptera':4,
    'Coleoptera':5,
    'Nematoda':6,
    'Hemiptera':7,
    'Human':8,
    'Hymenoptera':9,
    'Orthoptera':10}



    def make_heatmap_componens(event, comp_list, components_to_nodes_dict,host_orders_dict_for_nodes, domain):

        comps_stat = []
        for comp in comp_list:
            comp = int(comp)
            comp_hosts_ords = []

            for node in components_to_nodes_dict[comp]:
                comp_hosts_ords.extend(host_orders_dict_for_nodes[node]['Orders'])
                comp_hosts_ords_count = Counter(comp_hosts_ords)

            if len(comp_hosts_ords_count)>0:

                type_hosts_comp = [0]*11
                type_hosts_comp[0:3] = [domain+'_'+str(event), domain, comp]

                for host in comp_hosts_ords_count:
                    type_hosts_comp[hosts_ind_dict[host]] = comp_hosts_ords_count[host]

                comps_stat.append(type_hosts_comp)

        return(comps_stat)

    def make_simpson_comparisions(event, rec_tox, par_tox, rec_comps, par_comps, components_to_nodes_dict,host_orders_dict_for_nodes, domain, par_type):

        rec_ords = []
        rec_species = []
 
        par_ords = []
        par_species = []

        for comp in rec_comps:
            comp = int(comp)

            for node in components_to_nodes_dict[comp]:
                rec_ords.extend(host_orders_dict_for_nodes[node]['Orders'])
                rec_species.extend(host_orders_dict_for_nodes[node]['Species'])


        for comp in par_comps:
            comp = int(comp)

            for node in components_to_nodes_dict[comp]:
                par_ords.extend(host_orders_dict_for_nodes[node]['Orders'])
                par_species.extend(host_orders_dict_for_nodes[node]['Species'])

        if len(rec_species)>0 and len(par_species)>0:
            simp_orders = get_simpson_coeff(rec_ords, par_ords)
            simp_species = get_simpson_coeff(rec_species, par_species)

            jac_orders = get_jaccard_coeff(rec_ords, par_ords)
            jac_species = get_jaccard_coeff(rec_species, par_species)

            ret_list = [[event, par_type, domain, 'Orders', simp_orders, 'Simpson'],[event, par_type, domain, 'Species', simp_species, 'Simpson'],[event, par_type, domain, 'Orders', jac_orders, 'Jaccard'],[event, par_type, domain, 'Species', jac_species, 'Jaccard']]
            return(ret_list)
        #

    type_components_per_events_summary = [['Event','Tox','Tox_type', 'Component','Num_nodes', 'Comp_hosts']]
    type_components_per_events_simpson = [['Event','Tox_type','Dom_type','Host_type', 'Coeff', 'Coeff_type']]
    hosts_per_component = [['Event','Domain', 'Component', 'Lepidoptera','Diptera','Coleoptera','Nematoda','Hemiptera', 'Human', 'Hymenoptera', 'Orthoptera']]

    for event in tox_types_dict:
        recs = tox_types_dict[event]['Recs']
        rec_domain = rec_to_dom_dict[event]
        if 'Mins' in tox_types_dict[event]:
            mins = tox_types_dict[event]['Mins']
        else:
            mins = []

        if 'Majs' in tox_types_dict[event]:
            majs = tox_types_dict[event]['Majs']
        else:
            majs = []

        rec_comps = get_components_from_tox_list(recs, components_dict, tox_to_strain_dict)
        min_comps = get_components_from_tox_list(mins, components_dict, tox_to_strain_dict)
        maj_comps = get_components_from_tox_list(majs, components_dict, tox_to_strain_dict)

        rec_summary = make_comp_summary_for_tox(event, rec_comps, components_to_nodes_dict,host_orders_dict_for_nodes, 'Rec')
        type_components_per_events_summary.extend(rec_summary)

        min_comps_flag=False
        if min_comps:
            min_summary = make_comp_summary_for_tox(event, min_comps, components_to_nodes_dict,host_orders_dict_for_nodes, 'Min')
            type_components_per_events_summary.extend(min_summary)
            min_comps_flag = True

        maj_comps_flag=False
        if maj_comps:
            maj_summary = make_comp_summary_for_tox(event, maj_comps, components_to_nodes_dict,host_orders_dict_for_nodes, 'Maj')
            type_components_per_events_summary.extend(maj_summary)
            maj_comps_flag = True
       
        if not min_comps:
           min_comps = {}

        if not maj_comps:
            maj_comps = {}

        all_comps = list(set(list(rec_comps.values())+list(min_comps.values())+list(maj_comps.values())))
        comps_for_heatmap = make_heatmap_componens(event, all_comps, components_to_nodes_dict,host_orders_dict_for_nodes,rec_domain)


        if len(min_comps)>0:
            min_res = make_simpson_comparisions(event, recs, mins, rec_comps.values(), min_comps.values(), components_to_nodes_dict,host_orders_dict_for_nodes, rec_domain, 'Min')
            if min_res:
                type_components_per_events_simpson.extend(min_res)

        if len(maj_comps)>0:
            maj_res = make_simpson_comparisions(event, recs, majs, rec_comps.values(), maj_comps.values(), components_to_nodes_dict,host_orders_dict_for_nodes, rec_domain, 'Maj')
            if maj_res:
                type_components_per_events_simpson.extend(maj_res)

        if len(rec_summary[0][5]) >0 and (min_comps_flag or maj_comps_flag):
            if len(min_summary[0][5]) >0 or len(maj_summary[0][5]) >0:
                for comp_stat in comps_for_heatmap:
                    hosts_per_component.append(comp_stat)

    write_csv(rec_components_results, 'Graph_components_per_event_rank3_wide.csv')
    write_csv(rec_jaccard_long, 'Graph_components_per_event_rank3_long.csv')
    write_csv(components_hosts_res, 'host_num_per_components.csv')
    write_csv(type_components_per_events_summary, 'Per_event_summary_for_components.csv')
    write_csv(hosts_per_component, 'Components_heatmap_only_annotated.csv')
    write_csv(type_components_per_events_simpson, 'Overlap_species_components.csv')

    #for node in components_to_nodes_dict[0]:
    #    if len(set(host_orders_dict_for_nodes[node]['Orders']))>3:
    #        #orders_set = list(set(host_orders_dict_for_nodes[node]['Orders']))
    #        #if len(graph_dict[node])>300:
    #        #    print(node, orders_set, len(graph_dict[node]))
    #         for linked_node in 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A script for generating graph corresponding to strains')
    parser.add_argument('-r', '--rec', dest='rec_tab', help='the path to the recombination results',
                        type=str)
    parser.add_argument('-s', '--str', dest='str_tab', help='the path to the table with strains host specificity information',
                        type=str)
    parser.add_argument('-p', '--prot', dest='prot_tab', help='the path to the table with proteins used for graph building',
                        type=str)
    parser.add_argument('-o', '--spec', dest='spec_ordres', help='the path to the table denoting orders of affected species',
                        type=str)


    args = parser.parse_args()

    tox_df = pd.read_csv(os.path.realpath(args.str_tab), sep='\t')
    print(tox_df)
    tox_summary_dict = clean_tox_dict(make_tox_summary_dict(tox_df))


    tox_to_strain_dict, strain_to_tox_dict =  reduce_dicts_to_third_rank(*make_strains_dictionary(args.prot_tab))
    tox_to_strain_dict_rank4, strain_to_tox_dict_rank4 =  make_strains_dictionary(args.prot_tab)

    extended_tox_to_res_dict_rank3 = make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank3')
    extended_tox_to_res_dict_rank4 = make_extended_stain_tox_dict(strain_to_tox_dict_rank4, tox_summary_dict, dict_type='Rank4')

    spec_ord_df = pd.read_csv(os.path.realpath(args.spec_ordres), sep='\t')
    spec_ord_dict = pd.Series(spec_ord_df.Host_order.values,index=spec_ord_df.Host_sp).to_dict()

    graph_dict = prepare_dict_for_graph(tox_to_strain_dict, strain_to_tox_dict)
    color_dict, host_orders_dict_for_nodes = prepare_color_dict(graph_dict, extended_tox_to_res_dict_rank3, spec_ord_dict)

    G = build_strains_graph(graph_dict)


    iter_num = 0
    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    #draw_graph(S[0], color_dict, '0_scale125_pie_0025.png')


    for s in S:
        if s.number_of_nodes() >= 2 and s.number_of_nodes()<50: #s.number_of_nodes() >3 and s.number_of_nodes()<50
            print(s.number_of_nodes(), iter_num)
            name = str(iter_num)+'.svg'
    #        if iter_num not in [4,9,10,11,34,46, 27,13, 95,7,87]:
            draw_graph(s, color_dict, name)
        iter_num += 1

    make_connected_components_summary(args.rec_tab)


    file_reader = read_csv_to_list(os.path.realpath(args.rec_tab), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    analyze_connected_components(G, color_dict, rec_df, tox_to_strain_dict, strain_to_tox_dict, graph_dict, host_orders_dict_for_nodes)



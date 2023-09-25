import os
import math
import copy
import click
import itertools
import numpy as np
from numpy import radians as rad
import pandas as pd
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, RegularPolygon
from collections import defaultdict, Counter
from count_host_changes import get_parent_names
from matplotlib.patches import FancyArrowPatch, Circle
from build_recombination_graph import prepare_graph_dict
from annotate_clusters_consistency import write_csv, read_csv_to_list
from get_breakpoints_from_alignments import  get_major_and_minor_parents
from build_strains_graph import prepare_dict_for_graph, prepare_color_dict, build_strains_graph, get_components_from_tox_list
from overall_toxicity_summary import clean_tox_dict, reduce_dicts_to_third_rank, make_strains_dictionary, make_extended_stain_tox_dict, make_tox_summary_dict, crop_cry_name

global colors_dict
colors_dict = {'domain1':'#a60b0b', 'domain2':'#2980b9', 'domain3':'#bdbd00', 'domain12':'#556B2F', 'domain23':'#DDA0DD', 'domain13':'#71C685'}

def get_all_rec_toxins(rec_df):
    raw_tox_list = []

    for line in rec_df.iterrows():
        recs = get_parent_names(line[1]['Rec'])
        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        raw_tox_list.extend(minor_pars+recs+major_inter_pars)

    raw_tox_list = list(set([el for el in raw_tox_list if el!='unknown']))
    return(raw_tox_list)



def make_components_dict(G):

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

    return(components_dict)

def clean_strain_graph_from_non_rec(G, rec_df, graph_dict, tox_to_strain_dict):

    rec_crys = get_all_rec_toxins(rec_df)
    components_dict = make_components_dict(G)

    fixed_graph_dict = defaultdict(dict)
    comps_with_recs = list(set(get_components_from_tox_list(rec_crys, components_dict, tox_to_strain_dict).values()))

    strains_with_rec_comps = list(set([strain for strain in components_dict if str(components_dict[strain]) in comps_with_recs]))

    print('Initital components: {0}; filtered components: {1}'.format(len([G.subgraph(c).copy() for c in nx.connected_components(G)]), len(comps_with_recs)))

    for strain in graph_dict:
        if strain in strains_with_rec_comps:
            fixed_graph_dict[strain] = graph_dict[strain]

    print('Initital strains: {0}; filtered strains: {1}'.format(len(components_dict), len(fixed_graph_dict)))

    fixed_G = build_strains_graph(fixed_graph_dict)
    return(fixed_graph_dict, fixed_G)



def reduce_species_graph_to_components(fixed_graph_dict, fixed_G, color_dict_fixed, host_orders_dict_for_nodes_fixed):

    components_with_colors_dict = defaultdict(list)

    initial_components_dict = make_components_dict(fixed_G)
    size_components_dict = defaultdict()
    
    for strain in initial_components_dict:
        component = initial_components_dict[strain]
        colors = color_dict_fixed[strain]
        components_with_colors_dict[component].extend(colors)
        if component not in size_components_dict:
            size_components_dict[component] = 1
        else:
            size_components_dict[component] += 1

    G_comp = nx.MultiGraph()
    for node in components_with_colors_dict:
        G_comp.add_node(node)

    return(G_comp, components_with_colors_dict,size_components_dict, initial_components_dict)


def make_per_event_link_summary(rec_df):
    ids_tox_pairs_dict = dict()

    for line in rec_df.iterrows():
        recs = get_parent_names(line[1]['Rec'])
        rec_id = line[1]['ID']
        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        domain = line[1]['Type']
        min_type = colors_dict[domain]
        if domain=='domain1':
            maj_type = colors_dict['domain23']
        elif domain=='domain2':
            maj_type = colors_dict['domain13']
        elif domain=='domain3':
            maj_type = colors_dict['domain12']
  
        for pair in itertools.product(minor_pars,recs):
            if 'unknown' in pair:
                continue
            else:
               key = (crop_cry_name(pair[0]), crop_cry_name(pair[1]))
               ids_tox_pairs_dict[key] = (rec_id, 'minor',min_type, min_type)

        for pair in itertools.product(major_inter_pars,recs):
            if 'unknown' in pair:
                continue
            else:
               key = (crop_cry_name(pair[0]), crop_cry_name(pair[1]))
               ids_tox_pairs_dict[key] = (rec_id, 'major', maj_type, min_type)

    return(ids_tox_pairs_dict)



def make_rec_edges_data(G_rec, rec_df):
    colors_edges = nx.get_edge_attributes(G_rec, 'color').values()
    edges = G_rec.edges()
    rec_crys = get_all_rec_toxins(rec_df)
    
    directed_edge_dict = defaultdict(list)
    rec_evet_summary = make_per_event_link_summary(rec_df)


    return(rec_evet_summary)




def combine_graph_with_strains(G_comp, directed_edge_dict, component_tox_dict, strain_to_tox_dict):

    G_comp_directed = G_comp.copy()
    tox_in_comps = defaultdict(list)
    tox_to_comp = dict()

    for strain in component_tox_dict:
        component = component_tox_dict[strain]
        tox_in_comps[component].extend(strain_to_tox_dict[strain])
  
    for component in tox_in_comps:
        for tox in list(set(tox_in_comps[component])):
            tox_to_comp[tox] = component


    for edge in directed_edge_dict:

        if edge[0] in  tox_to_comp and edge[1] in tox_to_comp:

            comp_in = tox_to_comp[edge[0]]
            comp_out = tox_to_comp[edge[1]]

            G_comp_directed.add_edge(comp_in,comp_out,color=directed_edge_dict[edge][2],weight=16)

    passed_tuples_tox = []
    passed_tuples_types = []

    G_comp_reduced_types = G_comp.copy()
    G_comp_reduced_events = G_comp.copy()

    components_link_orders_all = dict()
    components_link_orders_reduced = dict()

    for edge in directed_edge_dict:

        if edge[0] in  tox_to_comp and edge[1] in tox_to_comp:

            comp_in = tox_to_comp[edge[0]]
            comp_out = tox_to_comp[edge[1]]
            event_link = (comp_in, comp_out, directed_edge_dict[edge][0:2]) 
            event_only_link = (comp_in, comp_out, directed_edge_dict[edge][0]) 

            if event_link not in passed_tuples_tox:
                G_comp_reduced_types.add_edge(comp_in,comp_out,color=directed_edge_dict[edge][2],weight=16)

                key = tuple(sorted([comp_in,comp_out])+[directed_edge_dict[edge][2]])
                components_link_orders_all[key] =(comp_in, comp_out)

            if event_only_link not in passed_tuples_types:
                G_comp_reduced_events.add_edge(comp_in,comp_out,color=directed_edge_dict[edge][3],weight=16)

                key = tuple(sorted([comp_in,comp_out])+[directed_edge_dict[edge][3]])
                components_link_orders_reduced[key] =(comp_in, comp_out)


            passed_tuples_tox.append(event_link)
            passed_tuples_types.append(event_only_link)


    print('Number components: ', len([G_comp_directed.subgraph(c).copy() for c in nx.connected_components(G_comp_directed)]))
    print('Number of initial parents edges: ', len(G_comp_directed.edges()))
    print('Number of reduced parents edges: ', len(G_comp_reduced_types.edges()))
    print('Number of reduced events edges: ', len(G_comp_reduced_events.edges()))


    return(G_comp_directed, G_comp_reduced_types, G_comp_reduced_events, components_link_orders_all, components_link_orders_reduced)


def drawCirc(fig,ax,radius,centX,centY,angle_,theta2_,color_='black'):

    arc = Arc([centX,centY],radius,radius*0.75,angle=angle_,
          theta1=0,theta2=theta2_,capstyle='round',linestyle='-',lw=1.5,color=color_)
    ax.add_patch(arc)


    #endX=centX+(radius/2)*np.cos(rad(theta2_+angle_)) 
    #endY=centY+(radius/2)*np.sin(rad(theta2_+angle_))

    #ax.add_patch(                   
    #    RegularPolygon(
    #        (endX, endY),            
    #        3,                       
    #        radius/9,                
    #        rad(angle_+theta2_),     
    #        color=color_
    #    )
    #)


def draw_network(fig, G,pos,ax,components_link_orders, size_dict, sg=None):

    for n in G:
        radius = math.log(size_dict[n]+3.5, 10)/43
        c=Circle(pos[n],radius=radius,alpha=1) #0.05
        G.node[n]['patch']=c
        x,y=pos[n]

    seen={}
    for (u,v,d) in G.edges(data=True):
        #print(sorted(u,v))
        key = (*sorted([u,v]),d['color'])
        rad1 = math.log(size_dict[u]+5.5, 10)/33
        rad2 = math.log(size_dict[v]+5.5, 10)/33

        u,v = components_link_orders[key]

        n1=G.node[u]['patch']
        n2=G.node[v]['patch']

        rad=0.1
        if (u,v) in seen:
            rad=seen.get((u,v))
            rad=(rad+np.sign(rad)*0.1)*-1
        alpha=0.5
        color=d['color']


        if n1!=n2:
            e = FancyArrowPatch(n1.center, n2.center,patchA=n1,patchB=n2,
                            arrowstyle='-',
                            connectionstyle='arc3,rad=%s'%rad,
                            mutation_scale=10.0,
                            lw=1.5,
                            color=color, zorder =100) #100

            seen[(u,v)]=rad
            ax.add_patch(e)
            e.set_zorder(550) #0

        else:
            circ_rad = 0.23+rad1
            drawCirc(fig, ax,circ_rad,n1.center[0], n1.center[1]-circ_rad/2,90,-360,color_=color)
    try:
        return e
    except:
        pass




def draw_combined_graph(new_G, color_dict, size_dict, components_link_orders, out_name):


    fig = plt.figure(figsize=(7, 7)) # 25 25 10 10
    ax = plt.axes([0, 0, 1, 1])
    ax.set_aspect('equal')

    nodes_num = len(new_G.nodes())

    if (nodes_num <7 and nodes_num not in [3,4]):
        k = 0.09*math.sqrt(nodes_num)  
        scale = 0.09*math.sqrt(nodes_num) 

    elif nodes_num in [3,4] :
        k = 0.15
        scale = 0.15

    elif nodes_num ==10:
        k = 1.45
        scale = 1.35

    else:
        k=2.5
        scale = 2.05

    pos = nx.spring_layout(new_G, iterations = 700, k=k, scale=scale, seed=30)  #4.1 5.25 #v
    degree_dict = dict(new_G.degree)


    weights_edges = list(nx.get_edge_attributes(new_G, 'weight').values())
    weights_colors = list(nx.get_edge_attributes(new_G, 'color').values())

    if nodes_num == 1:
        plt.xlim(-0.22, 0.22) 
        plt.ylim(-0.32, 0.12)

    elif nodes_num > 1 and nodes_num <8:
        plt.xlim(-0.42, 0.42) 
        plt.ylim(-0.42, 0.42)

    else:
        plt.xlim(-2.17, 2.17) 
        plt.ylim(-2.17, 2.17)

    trans = ax.transData.transform
    trans2 = fig.transFigure.inverted().transform

    print(out_name, nodes_num, scale, k, plt.xlim())
    for n in new_G:
        piesize = math.log(size_dict[n]+3.5, 10)/30 #40
        p2 = piesize/2.0
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

        a.pie(fracs, colors=color_set, wedgeprops={"edgecolor":"k",'linewidth': 1, 'linestyle': 'solid', 'antialiased': True, 'alpha':1}) #0.85


    draw_network(fig, new_G,pos,ax, components_link_orders, size_dict)

    plt.savefig(out_name)
    plt.close(fig)


def draw_components_separately(in_G, color_dict, size_dict, components_link_orders, out_name):
    it =0 
    S = [in_G.subgraph(c).copy() for c in nx.connected_components(in_G)]

    for Sub_s in S:
        out_subname =  str(it)+ '_' + out_name
        if len(Sub_s.edges())>=1:
            draw_combined_graph(Sub_s, color_dict, size_dict, components_link_orders, out_subname)
            it+=1



def prepare_fixed_strain_graph_data(rec_df, tox_df, prot_tab, spec_ordres, G_rec):

    tox_summary_dict = clean_tox_dict(make_tox_summary_dict(tox_df))
    tox_to_strain_dict, strain_to_tox_dict =  reduce_dicts_to_third_rank(*make_strains_dictionary(prot_tab))
    tox_to_strain_dict_rank4, strain_to_tox_dict_rank4 =  make_strains_dictionary(prot_tab)

    extended_tox_to_res_dict_rank3 = make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank3')
    extended_tox_to_res_dict_rank4 = make_extended_stain_tox_dict(strain_to_tox_dict_rank4, tox_summary_dict, dict_type='Rank4')

    spec_ord_df = pd.read_csv(os.path.realpath(spec_ordres), sep='\t')
    spec_ord_dict = pd.Series(spec_ord_df.Host_order.values,index=spec_ord_df.Host_sp).to_dict()

    graph_dict = prepare_dict_for_graph(tox_to_strain_dict, strain_to_tox_dict)
    color_dict, host_orders_dict_for_nodes = prepare_color_dict(graph_dict, extended_tox_to_res_dict_rank3, spec_ord_dict)

    G = build_strains_graph(graph_dict)
    fixed_graph_dict, fixed_G = clean_strain_graph_from_non_rec(G, rec_df, graph_dict, tox_to_strain_dict)
    color_dict_fixed, host_orders_dict_for_nodes_fixed = prepare_color_dict(graph_dict, extended_tox_to_res_dict_rank3, spec_ord_dict)

    G_comp, components_with_colors_dict,size_components_dict, component_tox_dict = reduce_species_graph_to_components(fixed_graph_dict, fixed_G, color_dict_fixed, host_orders_dict_for_nodes_fixed)

    directed_edge_dict = make_rec_edges_data(G_rec, rec_df)
    G_combined_multi, G_reduced_tox, G_comp_reduced_events, components_link_orders_all, components_link_orders_reduced = combine_graph_with_strains(G_comp, directed_edge_dict, component_tox_dict, strain_to_tox_dict)

    #draw_components_separately(G_combined_multi, components_with_colors_dict,size_components_dict, components_link_orders_all, 'toxins.svg')
    draw_components_separately(G_reduced_tox, components_with_colors_dict,size_components_dict, components_link_orders_all , 'parents_reduced_parents.svg')
    draw_components_separately(G_comp_reduced_events, components_with_colors_dict,size_components_dict, components_link_orders_reduced , 'events_reduced_events.svg')




@click.command()
@click.option('--rec_file', '-r', help="recombination results", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--str_tab', '-s', help="the path to the table with strains host specificity information", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--prot_tab', '-p', help="the path to the table with proteins used for graph building", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--spec_ordres', '-o', help="the path to the table denoting orders of affected species", 
              type=click.Path(exists=True), metavar='<PATH>') 



def main(rec_file, str_tab, prot_tab, spec_ordres):

    tox_df = pd.read_csv(os.path.realpath(str_tab), sep='\t')
    G_rec, graph_dict_rec, node_color_dict  = prepare_graph_dict(rec_file, rec_data = 'filt')



    file_reader = read_csv_to_list(os.path.realpath(rec_file), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    prepare_fixed_strain_graph_data(rec_df, tox_df, prot_tab, spec_ordres, G_rec)
    
    

if __name__ == '__main__':
   main()

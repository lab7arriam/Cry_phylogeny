import os
import copy
import click
import os.path
import random
import math
import numpy as np
from scipy.stats import t
import pandas as pd
from ete3 import Tree
from itertools import product, combinations
from scipy.stats import fisher_exact
from collections import defaultdict
from annotate_clusters_consistency import write_csv, read_csv_to_list
from overall_toxicity_summary import make_tox_summary_dict,crop_cry_name, clean_tox_dict, make_strains_dictionary, get_parent_names, check_tox_for_hosts, make_extended_stain_tox_dict, get_hosts_list, extend_defaultdict_to_lists
from get_breakpoints_from_alignments import  get_major_and_minor_parents, unpacking_list_of_lists
from create_evol_data_on_fixed_recombinants import select_tree, gca, get_identity_heatmap


def get_common_and_non_common_hosts(hosts_of_recs_dict, hosts_of_parents_dict, common_dict, rec_type):
    hosts_of_parents = extend_defaultdict_to_lists(hosts_of_parents_dict)
    hosts_of_recs =  extend_defaultdict_to_lists(hosts_of_recs_dict)
    common_hosts = set(hosts_of_recs).intersection(set(hosts_of_parents))


    if len(hosts_of_parents)==0 or len(hosts_of_recs)==0:
        return(0)

    for item in hosts_of_parents:
        if item in common_hosts:
            common_dict[rec_type]['Par_common'] += 1
        else:
            common_dict[rec_type]['Par_uncommon'] += 1


    for item in hosts_of_recs:
        if item in common_hosts:
            common_dict[rec_type]['Rec_common'] += 1
        else:
            common_dict[rec_type]['Rec_uncommon'] += 1


    if set(hosts_of_recs)!=set(hosts_of_parents):
        return(1)


def mearge_defaultdicts_lists(dict_list1, dict_list2, mode ='list'):
    dict_list_mearged=defaultdict(list)
    key_set1 = set([key for key in dict_list1 if key not in dict_list2]) #keys from the first dictionary
    key_set2 = set([key for key in dict_list2 if key not in dict_list1]) #keys from the second dictionary
    key_set3 = set([key for key in dict_list1] + [key for key in dict_list2]) - key_set1 - key_set2 #shared keys
    
    for key1 in key_set1:
        dict_list_mearged[key1] = dict_list1[key1]

    for key2 in key_set2:
        dict_list_mearged[key2] = dict_list2[key2]

    for key3 in key_set3:
        if mode == 'list':
            dict_list_mearged[key3] = key_set1[key3] + dict_list2[key3]
        else:
            for subkey in dict_list1[key3].keys():
                dict_list_mearged[key3][subkey] = dict_list1[key3][subkey]+dict_list2[key3][subkey]

    return(dict_list_mearged)


def get_host_list_by_mode(tox_list, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type, Rank, mode):

    hosts = get_hosts_list(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type, Rank)
    str_hosts = get_hosts_list(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type, 'Strain')
    str_hosts_extr = get_hosts_list(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type, 'Strain_ext')

    if 'add' not in mode and 'ext' not in mode :
        return(hosts)

    elif 'add' in mode and 'ext' not in mode:
        if len(hosts)!=0:
            #return(mearge_defaultdict_lits(hosts,str_hosts)) try mearding all data from trains
            return(hosts)
        else:
            return(str_hosts)

    elif 'add' in mode and 'ext' in mode:
    
        if len(hosts)!=0:
            return(hosts)
        elif len(str_hosts)!=0:
            return(str_hosts)
        else:
            return(str_hosts_extr)

    else:
        return(str_hosts_extr)



def extend_dict_to_orders(host_res, spec_ord_dict):
    annot_res = []
    for host in  host_res:
        annot_res.append(host + '('+spec_ord_dict[host]+')')

    return(', '.join(sorted(annot_res)))


def parse_host_dict_old(host_dict, spec_ord_dict, host_type =  'Orders', Rank = 'Rank3'):

    res_list = []

    if len(host_dict)==0:
        return('Unknown')
     
    passed = set()

    for tox in host_dict:
        if Rank == 'Rank3':
            crop_tox = crop_cry_name(tox)

            if crop_tox in passed:
                continue
            if host_type == 'Species':

                res_list.append(crop_tox+": {" + extend_dict_to_orders(host_dict[tox], spec_ord_dict)+'}')

            else:
                res_list.append(crop_tox+": {" + ', '.join(sorted(host_dict[tox]))+'}')  

            passed.add(crop_tox)

        else:
            if host_type == 'Species':

                res_list.append(tox+": {" + extend_dict_to_orders(host_dict[tox], spec_ord_dict)+'}')

            else:
                res_list.append(tox+": {" + ', '.join(sorted(host_dict[tox]))+'}')  

    return('; '.join(res_list))


def parse_tox_for_heatmap(host_dict, spec_ord_dict, host_type =  'Orders', mode = 'heatmap'):

    res = []
    if len(host_dict)==0:
        return(['Unknown'])

    #print(host_dict)

    if mode == 'heatmap':

        for tox in host_dict:
            host_list = host_dict[tox]
            for host in host_list:
                if host_type=='Species':
                    res.append(spec_ord_dict[host])
                else:
                    res.append(host)

        return(res)

    else:
        res_dict = defaultdict(list)

        for tox in host_dict:
            host_list = host_dict[tox]
            for host in host_list:
                if host_type=='Species':
                    res_dict[tox].append(host + ' ('+spec_ord_dict[host]+')')
                else:
                    res.append(host)
        return(res_dict)


def make_tab_for_heatmap(hosts_of_recs, hosts_of_min_parents, hosts_of_maj_parents, spec_ord_dict, host_type =  'Orders'):
    hosts_lists = []
    for rec_host in parse_tox_for_heatmap(hosts_of_recs, spec_ord_dict, host_type ):
        hosts_lists.append(['Rec', rec_host])

    for min_host in parse_tox_for_heatmap(hosts_of_min_parents, spec_ord_dict, host_type ):
        hosts_lists.append(['Min', min_host])

    for maj_host in parse_tox_for_heatmap(hosts_of_maj_parents, spec_ord_dict, host_type ):
        hosts_lists.append(['Maj', maj_host])

    return(hosts_lists)

def make_dict_for_text(annot_dict, tox_list, tox_to_strain_dict):
    
    strains_ext_dict = {tox:tox_to_strain_dict[tox] for tox in tox_list if tox in tox_to_strain_dict}
    annot_dict_extended = defaultdict(list)


    if annot_dict==['Unknown']:
        return(annot_dict)
 
    for tox in tox_list:
        if tox in annot_dict:
            annot_dict_extended[tox]= sorted(annot_dict[tox])
        else:
            if tox not in tox_to_strain_dict:
                annot_dict_extended[tox] = ['Unknown']

            else:
                strain_hosts = []
                strains_with_annots = []
                for strain in tox_to_strain_dict[tox]:
                    if strain in annot_dict:
                        strains_with_annots.append(strain)
                        strain_hosts = strain_hosts + annot_dict[strain]

                strain_hosts = sorted(list(set(strain_hosts)))
                if len(strain_hosts) >0:
                    annot_dict_extended[tox + ' ('+  '; '.join(strains_with_annots)+')'] = strain_hosts
                  

                else:
                    annot_dict_extended[tox] = ['Unknown']


    if len(annot_dict_extended)!= len(tox_list):
        print(tox_list, annot_dict_extended)
    return(annot_dict_extended)  

def get_keys_num(ext_dict):
    if ext_dict==['Unknown']:
        return(0)
    else:
        return(sum([len(ext_dict[key]) for key in ext_dict]))

def extract_hosts_block_from_host_dict(host_dict, max_key_num):
    ret_list = [['',''] for el in range(max_key_num)]
    counter = 0

    if host_dict==['Unknown']:
        return(ret_list)

    else:
        num_keys = len(host_dict)
        host_dict_no_unknowns = ({(key):(host_dict[key] if host_dict[key]!=['Unknown'] else []) for key, value in host_dict.items()})
        host_dict_sorted = ({k: v for k, v in sorted(host_dict_no_unknowns.items(), key=lambda item: len(item[1]), reverse=True)})

        for tox in host_dict_sorted:
            passed_flag = 0

            if host_dict_sorted[tox] ==[]:
                ret_list[counter] = [tox, '']
                counter += 1

            for host in host_dict_sorted[tox]:
                if passed_flag == 0:
                    ret_list[counter] = [tox, host]
                else:
                    ret_list[counter] = ['',host]

                counter += 1
                passed_flag = 1

    return(ret_list)



def make_table_for_text(hosts_recs_species, hosts_min_species, hosts_maj_species, spec_ord_dict, recs, minor_pars, major_inter_pars, tox_to_strain_dict):
    if len(major_inter_pars) ==0:
        major_inter_pars = ['Unknown']

    recs_annot_hosts_dict = parse_tox_for_heatmap(hosts_recs_species, spec_ord_dict,'Species', 'text')
    mins_annot_hosts_dict = parse_tox_for_heatmap(hosts_min_species, spec_ord_dict,'Species', 'text')
    majs_annot_hosts_dict = parse_tox_for_heatmap(hosts_maj_species, spec_ord_dict,'Species', 'text')

    recs_extended = make_dict_for_text(recs_annot_hosts_dict, recs, tox_to_strain_dict)
    mins_extended = make_dict_for_text(mins_annot_hosts_dict, minor_pars, tox_to_strain_dict)
    majs_extended = make_dict_for_text(majs_annot_hosts_dict, major_inter_pars, tox_to_strain_dict)
   
    
    max_row_num = max([get_keys_num(recs_extended),get_keys_num(mins_extended),get_keys_num(majs_extended)])
    
    rec_block = extract_hosts_block_from_host_dict(recs_extended, max_row_num)
    min_block = extract_hosts_block_from_host_dict(mins_extended, max_row_num)
    maj_block = extract_hosts_block_from_host_dict(majs_extended, max_row_num)

    return([rec_block, min_block, maj_block])

def make_first_block(num_rows, rec_id, rec_type, mode, same_flag):
    ret_list = [['','', '', ''] for el in range(num_rows)]
    ret_list[0] = [rec_id, rec_type, mode, same_flag]
    return(ret_list)


def mearge_event_blocks(block_lists):
    mearged_list = [[] for el in range(len(block_lists[0]))]

    for i in range(len(mearged_list)):
        mearged_row = []
        for inner_ind in range(4):
            mearged_row = mearged_row + block_lists[inner_ind][i]
        mearged_list[i] = mearged_row

    return(mearged_list)

def parse_tox_df(rec_df, tox_summary_dict, tox_to_strain_dict, strain_to_tox_dict, spec_ord_dict, mode = 'Rank4', Rank = 'Rank4', host_type = 'Orders', par_type = 'Strict'):

    overall_tox_res_per_events_for_text = []
    overall_tox_res_per_events_for_heatmap = []
    mearded_hetamap_for_text = [['ID', 'Type', 'Mode', 'Same_flag', 'Recs','Rec_hosts', 'Mins','Min_hots', 'Majs', 'Maj_hosts']]
    heatmap_text_name = mode + '_text_heatmap.csv'

    if 'ext_3' in mode:
        extended_stain_tox_dict = make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank3')
    else:
        extended_stain_tox_dict = make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank4')

    common_dict = defaultdict(dict)

    for dom_type in ['domain1', 'domain2','domain3']:
        common_dict[dom_type] = {'Par_common':0,'Par_uncommon':0, 'Rec_common':0, 'Rec_uncommon':0}

    same_events = []
    non_same_events = [] 
    num_unknowns = []

    for line in rec_df.iterrows():
       
        all_hosts_of_event = set()
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']
        rec_unknown = line[1]['Unknown_flag']
        
        #if rec_id!='109':
        #    continue

        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        min_maj_strict_pars = minor_pars+major_inter_pars
        min_maj_all_pars = minor_pars+major_all_pars


        hosts_of_recs = get_host_list_by_mode(recs,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type,Rank, mode)
        hosts_of_min_parents = get_host_list_by_mode(minor_pars,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type,Rank, mode)
        hosts_of_maj_parents = get_host_list_by_mode(major_inter_pars,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type,Rank, mode)

        hosts_of_all_parents_strict = get_host_list_by_mode(min_maj_strict_pars,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type,Rank, mode)
        hosts_of_all_parents_all = get_host_list_by_mode(min_maj_all_pars,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type,Rank, mode)

        hosts_recs_species = get_host_list_by_mode(recs,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, 'Species',Rank, mode)
        hosts_min_species = get_host_list_by_mode(minor_pars,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, 'Species',Rank, mode)
        hosts_maj_species = get_host_list_by_mode(major_inter_pars,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, 'Species',Rank, mode)

        if par_type == 'Strict':

             #if rec_id=='52':
             #   print(Rank, mode, host_type, par_type, rec_id, rec_type, hosts_of_recs, hosts_of_all_parents_strict)

             res = get_common_and_non_common_hosts(hosts_of_recs, hosts_of_all_parents_strict, common_dict, rec_type)

             if len(hosts_of_recs)>0:
                 if len(hosts_of_all_parents_strict)>0:

                     heatmap_diff = make_tab_for_heatmap(hosts_recs_species, hosts_min_species, hosts_maj_species, spec_ord_dict, 'Species')

                     text_tox_res_list = make_table_for_text(hosts_recs_species, hosts_min_species, hosts_maj_species, spec_ord_dict, recs, minor_pars, major_inter_pars, tox_to_strain_dict)


                     if (res)==1:
                         first_text_block = make_first_block(len(text_tox_res_list[0]), rec_id, rec_type, mode, 'Different')
                         mearged_res = mearge_event_blocks([first_text_block]+text_tox_res_list)
                         mearded_hetamap_for_text.extend(mearged_res)


                         non_same_events.append(rec_id)
                         if rec_unknown!='no':
                             num_unknowns.append(rec_id)

                         
                         for row in heatmap_diff:
                             overall_tox_res_per_events_for_heatmap.append([rec_id, 'Different', rec_type]+ row+[mode, len(recs), len(minor_pars), len(major_inter_pars)])

                     else:
                         same_events.append(rec_id)
                         first_text_block = make_first_block(len(text_tox_res_list[0]), rec_id, rec_type, mode, 'Same')

                         mearged_res = mearge_event_blocks([first_text_block]+text_tox_res_list)
                         mearded_hetamap_for_text.extend(mearged_res)

                         if rec_unknown!='no':
                             num_unknowns.append(rec_id)

                         for row in heatmap_diff:
                             overall_tox_res_per_events_for_heatmap.append([rec_id, 'Same',  rec_type]+ row+[mode, len(recs), len(minor_pars), len(major_inter_pars)])


        elif par_type == 'All':
            #if rec_id=='2':
                #print(Rank, mode, host_type, par_type, rec_id, rec_type)
             res = get_common_and_non_common_hosts(hosts_of_recs, hosts_of_all_parents_all, common_dict, rec_type)

    print(mode,host_type , 'Same:', len(same_events), 'Different', len(non_same_events), 'Unknowns:', len(num_unknowns), 'Total:', len(same_events+non_same_events))
    dom_sum = sum([sum(list(common_dict['domain1'].values())), sum(list(common_dict['domain2'].values())), sum(list(common_dict['domain3'].values()))])

    write_csv(mearded_hetamap_for_text, heatmap_text_name)
    return((common_dict, overall_tox_res_per_events_for_heatmap))



def select_tree_for_length_distr(rec_type, trees) -> Tree:

    if rec_type.startswith('domain1'):
        return((trees[0],trees[1]))
    elif rec_type.startswith('domain2'):
        return [trees[1]]
    elif rec_type.startswith('domain3'):
        return ((trees[2], trees[1]))
    else:
        return None


def get_all_subtrees_from_tree(tree):
    node_content = tree.get_cached_content(store_attr='name')

    dist_list = []
    for node in node_content:
        if len(node)==368: 
            dist_list.append([len(node), len(node), tree.get_farthest_node()[1]])
        else:
            dist_list.append([len(node), len(node), node.dist])

    return(dist_list)
       

def analyze_tree_distance(rec_type, recs, minor_pars, major_pars, trees):
    dist_res = []

    adj_domain = ''
    if rec_type== 'domain1' or rec_type== 'domain3':
        adj_domain = 'domain2'
    
    if rec_type=='domain2' and (len(minor_pars)==0 or minor_pars==['unknown']):
        return(None)

    elif rec_type=='domain2':
        return(None)

    else:
        if not minor_pars==['unknown'] or len(minor_pars)==0:
            tree_transf, tree_adj = select_tree_for_length_distr(rec_type, trees)
            anc_transf = gca(tree_transf, recs+minor_pars)
            anc_dist =  anc_transf.dist

            if len(anc_transf)==368:
                anc_dist = anc_transf.get_farthest_node()[1]

            dist_res.append([rec_type, rec_type, 'Min',len(recs+minor_pars), len(anc_transf), anc_dist ])


        if not major_pars==['unknown'] or len(major_pars)==0:
            tree_transf, tree_adj = select_tree_for_length_distr(rec_type, trees)
            anc_transf = gca(tree_adj, recs+major_pars)
            anc_dist =  anc_transf.dist

            if len(anc_transf)==368:
                anc_dist = anc_transf.get_farthest_node()[1]

            dist_res.append([rec_type, adj_domain, 'Maj',len(recs+major_pars), len(anc_transf),  anc_dist])

    if len(dist_res)!=0:
        return(dist_res)


def make_distibution_of_tree_lengths(rec_df, tree_list:list):

    tree_dist_res = [['ID', 'Type', 'Tree_domain', 'Parent_type', 'Num_tox', 'Num_nodes', 'Dist']]
    ind=1
    for tree in tree_list:
        dist_res = get_all_subtrees_from_tree(tree)
        domain = 'domain' +str(ind)
        for dist_row in dist_res:
            tree_dist_res.append(['Ref', 'domain1', domain, 'Ref'] + dist_row)
            tree_dist_res.append(['Ref', 'domain2', domain, 'Ref'] + dist_row)
            tree_dist_res.append(['Ref', 'domain3', domain, 'Ref'] + dist_row)
        ind+=1
        
    set_pars_strict = set()
    set_pars_all = set()
    #n_passed,n_corrupted = 0,0

    for line in rec_df.iterrows():
       
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']

        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')


        set_pars_strict = set_pars_strict.union(set(recs+minor_pars+major_inter_pars))
        set_pars_all = set_pars_all.union(set(recs+minor_pars+major_all_pars))
        res = analyze_tree_distance(rec_type, recs,  minor_pars, major_inter_pars, tree_list)
        if res:
            for dom_res in res:
                tree_dist_res.append([rec_id]+dom_res)

    print(len(set_pars_strict), len(set_pars_all))

    #write_csv(tree_dist_res, 'Distribution_of_branch_lengths_no_domain2.csv')

def analyze_trees_lengths(rec_df, tree_path):
    trees = []
    for tf in sorted(os.listdir(os.path.realpath(tree_path))):
        tf_name = os.path.join(os.path.realpath(tree_path), tf)
        t = Tree(tf_name)
        trees.append(t)
    
    make_distibution_of_tree_lengths(rec_df, trees)


def get_mean_and_CL_from_list(x):
    x = np.array(x)
    m = x.mean() 
    s = x.std() 
    dof = len(x)-1 
    confidence = 0.95
    t_crit = np.abs(t.ppf((1-confidence)/2,dof))
    interval = s*t_crit/np.sqrt(len(x))

    return((round(m,3), round(interval,3)))


def get_domain_for_major(rec_type):
    if rec_type == 'domain3':
        dom_ind = 'domains12'  
    elif rec_type == 'domain2':
        dom_ind = 'domains13' 
    elif rec_type == 'domain1':
        dom_ind = 'domains23' 

    return(dom_ind)

def get_descriptive_data_for_par(rec_type, recs, pars, heatmap_dict, par_type = 'minor'):
    dom_ind = rec_type
    ret_num = 0
    ret_id_list = []

    if par_type=='major':
        dom_ind = get_domain_for_major(rec_type)

    if pars!=['unknown'] and len(pars) > 0:
        ret_num = len(pars)
        for tox_tuple in list(product(recs, pars)):
            key = '|'.join(sorted(list(tox_tuple)))
            ret_id_list.append(heatmap_dict[key][dom_ind])

        return((ret_num, ret_id_list))


def get_simpson_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/min(len(set(list1)),len(set(list2))))


def get_all_hosts_by_mode(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, mode):

    if mode =='Strains_add':
        host_list = get_hosts_list(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, 'Species', 'Rank3')
        if len(host_list) == 0:
            host_list = get_hosts_list(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, 'Species', 'Strain')
            return(list(set(unpacking_list_of_lists([host_list[tox] for tox in host_list]))))
        else:
            return(list(set(unpacking_list_of_lists([host_list[tox] for tox in host_list]))))

    else: 
        host_list = get_hosts_list(tox_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, 'Species', mode)

        return(list(set(unpacking_list_of_lists([host_list[tox] for tox in host_list]))))


def generate_coeff_for_tox_lists(rec_list, par_list, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict):
    results = []
    for mode in ['Rank4', 'Rank3', 'Strains_add']:
        hosts_rec = get_all_hosts_by_mode(rec_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, mode)
        host_par = get_all_hosts_by_mode(par_list,tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, mode)
        if len(hosts_rec)==0 or len(host_par)==0:
            continue
        else:
            simp = get_simpson_coeff(hosts_rec, host_par)
            results.append([mode, simp])
    if len(results)>0:
        return(results)
    else:
        return(None)

def summarize_nums_and_ids_in_events(rec_df, heatmap_dict, tox_summary_dict, tox_to_strain_dict, strain_to_tox_dict):
    
    ids_min = []
    ids_maj = []
    nums = [[],[],[]]


    rec_descriptive_stat = defaultdict(dict)

    extended_stain_tox_dict = make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank3')


    
    for mode in ['Rank3', 'Rank4', 'Strains_add']:
        rec_descriptive_stat[mode] = defaultdict(dict)
        for domain in ['domain1', 'domain2', 'domain3']:
            rec_descriptive_stat[mode][domain] = defaultdict(dict)
            for par_type in ['rec', 'min','maj']:
                rec_descriptive_stat[mode][domain][par_type] = defaultdict(list)
                rec_descriptive_stat[mode][domain][par_type]['id']=[]
                rec_descriptive_stat[mode][domain][par_type]['num']=[]


    for_simpson_calc = [['ID', 'Domain', 'Tox_type', 'Tox_num', 'Ident', 'Mode', 'Simp_coeff']]


    for line in rec_df.iterrows():
       
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']

        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        
        min_res = get_descriptive_data_for_par(rec_type, recs, minor_pars, heatmap_dict, 'minor')
        maj_res = get_descriptive_data_for_par(rec_type, recs, major_inter_pars, heatmap_dict, 'major')
        #all_res = get_descriptive_data_for_par(rec_type, recs, minor_pars+major_inter_pars, heatmap_dict, 'major')

        if min_res:
            res = generate_coeff_for_tox_lists(recs, minor_pars, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict)
            if res: 
                for row in res:
                    mode = row[0]
                    rec_descriptive_stat[mode][rec_type]['min']['num'].append(min_res[0])
                    rec_descriptive_stat[mode][rec_type]['min']['id'].extend(min_res[1])
                    for_simpson_calc.append([rec_id, rec_type, 'min',min_res[0] ,get_mean_and_CL_from_list(min_res[1])[0]] +row)         

        if maj_res:

            res = generate_coeff_for_tox_lists(recs, major_inter_pars, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict)

            if res: 
                for row in res:
                    mode = row[0]
                    rec_descriptive_stat[mode][rec_type]['maj']['num'].append(maj_res[0])
                    rec_descriptive_stat[mode][rec_type]['maj']['id'].extend(maj_res[1])

                    for_simpson_calc.append([rec_id, rec_type, 'maj',maj_res[0] ,get_mean_and_CL_from_list(maj_res[1])[0]] +row) 
                    

        all_res = generate_coeff_for_tox_lists(recs, major_inter_pars+minor_pars, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict)
        if all_res:
            for row in all_res:
                for_simpson_calc.append([rec_id, rec_type, 'all' ,'-','-'] +row)
        #if min_res or maj_res:
        #    rec_descriptive_stat[rec_type]['rec']['num'].append(len(recs))
            #if len(recs)>1:
            #    for tox_tuple in list(combinations(recs, 2)):
            #        key = '|'.join(sorted(list(tox_tuple)))
            #        rec_descriptive_stat[rec_type]['rec']['id'].append((heatmap_dict[key]['domain1']+heatmap_dict[key]['domain2']+heatmap_dict[key]['domain3'])/3)

            #else:
            #    rec_descriptive_stat[rec_type]['rec']['id'].append(100)


    write_csv(for_simpson_calc, 'Simpson_stat_for_events.csv')

    rec_descriptive_overall_stat = copy.deepcopy(rec_descriptive_stat)

    for mode in rec_descriptive_stat:
        for rec_type in rec_descriptive_stat[mode]:
            for par in rec_descriptive_stat[mode][rec_type]:
                if len(rec_descriptive_stat[mode][rec_type][par]['num'])==0:
                    if par in rec_descriptive_overall_stat[mode][rec_type]:
                        del(rec_descriptive_overall_stat[mode][rec_type][par])

    for mode in rec_descriptive_stat:
        for measure in ['id', 'num']:
            for rec_type in rec_descriptive_stat[mode]:
                if 'rec' in rec_descriptive_overall_stat[mode][rec_type]:
                    del rec_descriptive_overall_stat[mode][rec_type]['rec']

                if rec_type != 'rec': 
                    num_events_min = len(rec_descriptive_stat[mode][rec_type]['min']['num'])
                    if num_events_min>0:
                        pass
                        #print(mode, measure, rec_type,'min',num_events_min, str(get_mean_and_CL_from_list(rec_descriptive_stat[mode][rec_type]['min'][measure])[0])+' +/- '+str(get_mean_and_CL_from_list(rec_descriptive_stat[mode][rec_type]['min'][measure])[1]), sep ='\t')


                    num_events_maj = len(rec_descriptive_stat[mode][rec_type]['maj']['num'])
                    if num_events_maj>0:
                        pass
                        #print(mode, measure, rec_type,'maj',num_events_maj, str(get_mean_and_CL_from_list(rec_descriptive_stat[mode][rec_type]['maj'][measure])[0])+' +/- '+str(get_mean_and_CL_from_list(rec_descriptive_stat[mode][rec_type]['maj'][measure])[1]), sep ='\t')


                
    return(rec_descriptive_overall_stat)



def get_reduced_tox_list(tox_summary_dict, ref_list, tox_to_strain_dict, strain_to_tox_dict, extended_stain_tox_dict, mode):
    prots_to_use = []

    for cry in ref_list:
        hosts = get_all_hosts_by_mode([cry],tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, mode)
        if len(hosts)>0:
            prots_to_use.append(cry)

    return(prots_to_use)

def reduce_heatmap_by_list(prot_list, heatmap_dict):
    ret_heatmap = copy.deepcopy(heatmap_dict)

    for key in heatmap_dict:
        tox1, tox2 = key.split('|')
        if not (tox1 in prot_list and tox2 in prot_list):
            del ret_heatmap[key]

    return(ret_heatmap)


def generate_random_prot_pairs(heatmap_dict, domain, id_threshold, t, num):

    pairs_list = set()
    while len(pairs_list) < num:
       sample_pair = random.sample(heatmap_dict.keys(), 1)[0]
       if heatmap_dict[sample_pair][domain] <= id_threshold+t and heatmap_dict[sample_pair][domain] >= id_threshold-t:
           pairs_list.add(sample_pair)

    return(sorted(list(pairs_list)))


def calculate_num_of_pairs_for_dataset(heatmap_dict, domain, id_threshold, t):
    pairs_list = []
    for tox_pair in heatmap_dict:
       if heatmap_dict[tox_pair][domain] <= 100 and heatmap_dict[tox_pair][domain] >= id_threshold-t: # id_threshold+t
           pairs_list.append(tox_pair)


    num_prots = list(set(unpacking_list_of_lists([tox_pair.split('|') for tox_pair in pairs_list])))
    #print(pairs_list)

    return(len(pairs_list), len(num_prots))



def count_host_specificity(rec_table, tox_summary_dict, strains_dat, spec_ord_dict, heatmap_dict, ref_prots): #tree_path

    file_reader = read_csv_to_list(os.path.realpath(rec_table), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])
    tox_to_strain_dict, strain_to_tox_dict =  make_strains_dictionary(strains_dat)

    #tree stat (not valid)
    #analyze_trees_lengths(rec_df, tree_path)

    #modes = ['Rank4', 'Rank3', 'Strains_add',  'Strains_ext_4_add', 'Strains_ext_3_add', 'Strains_ext_4_all', 'Strains_ext_3_all'] #all comparision modes
    fisher_calcs_resulst = [['Host_type','Mode', 'Parent_type', 'Par_common', 'Par_uncommon', 'Rec_common', 'Rec_uncommon', 'domain','p_value']] #fisher results
    modes = [ 'Rank4', 'Rank3','Strains_add']  #'Rank4', 'Rank3',

    data_for_heatmap_res = [['ID', 'Change_flag', 'Domain', 'Tox_type', 'Host', 'Mode', 'Num_recs', 'Num_mins', 'Num_majs']] # data for hetmap with pie charts


    for mode in modes:
        if mode =='Rank4':
            Rank=mode
        else:
            Rank='Rank3'
        
        for host_type in ['Orders']: #'Orders','Species'
            for par_type in ['Strict']: #, 'All'
                print(mode)
                comm_dict, data_for_heatmap_tox = parse_tox_df(rec_df, tox_summary_dict, tox_to_strain_dict, strain_to_tox_dict, spec_ord_dict, mode, Rank, host_type, par_type)
                data_for_heatmap_res.extend(data_for_heatmap_tox)
                for domain in comm_dict:
                    res = fisher_exact([[comm_dict[domain]['Par_common'], comm_dict[domain]['Par_uncommon']], [comm_dict[domain]['Rec_common'], comm_dict[domain]['Rec_uncommon']]])
                    comm_res_list = [comm_dict[domain]['Par_common'], comm_dict[domain]['Par_uncommon'], comm_dict[domain]['Rec_common'], comm_dict[domain]['Rec_uncommon']]

                    fisher_calcs_resulst.append([host_type, mode, par_type] +comm_res_list +[domain, res[1]])
 
    #write_csv(fisher_calcs_resulst, 'Fisher_results_for_different_modes.csv')
    write_csv(data_for_heatmap_res, 'Hosts_for_heatmap_stat_per_species.csv')

    rec_descriptive_stats = summarize_nums_and_ids_in_events(rec_df, heatmap_dict, tox_summary_dict, tox_to_strain_dict, strain_to_tox_dict)
    extended_stain_tox_dict = make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank3')


    #reduced_tox_list = get_reduced_tox_list(tox_summary_dict, ref_prots , tox_to_strain_dict, strain_to_tox_dict, extended_stain_tox_dict, 'Strains_add')
    #reduced_heatmap = reduce_heatmap_by_list(reduced_tox_list, heatmap_dict)

    for mode in rec_descriptive_stats:
        reduced_tox_list = get_reduced_tox_list(tox_summary_dict, ref_prots , tox_to_strain_dict, strain_to_tox_dict, extended_stain_tox_dict, mode)
        reduced_heatmap = reduce_heatmap_by_list(reduced_tox_list, heatmap_dict)
        for domain in rec_descriptive_stats[mode]:
            for parent_type in rec_descriptive_stats[mode][domain]:
                num_events = len(rec_descriptive_stats[mode][domain][parent_type]['num'])
                mean_num, num_CI = get_mean_and_CL_from_list(rec_descriptive_stats[mode][domain][parent_type]['num'])
                mean_ident, ident_CI = get_mean_and_CL_from_list(rec_descriptive_stats[mode][domain][parent_type]['id'])

                if math.isnan(ident_CI):
                    ident_CI = mean_ident*0.05

                if parent_type == 'min':
                    num_pairs, num_prots = calculate_num_of_pairs_for_dataset(reduced_heatmap, domain, mean_ident, ident_CI)
                    print(mode, domain, parent_type, mean_ident, ident_CI, num_pairs, num_prots, sep = '\t')

                if parent_type == 'maj':
                    major_domain = get_domain_for_major(domain)
                    num_pairs, num_prots = calculate_num_of_pairs_for_dataset(reduced_heatmap, major_domain, mean_ident, ident_CI)

                    print(mode, domain, parent_type, mean_ident, ident_CI, num_pairs, num_prots, sep = '\t')
      
           
        #print(rec_descriptive_stats[mode]


    #for iter_num in range(50):

    #    tox_list = generate_random_prot_pairs(reduced_heatmap, 'domain3', 71.676, 2.44, 7)
    #    print(tox_list)

    #for tox pair in


@click.command()           
@click.option('--rec_tab', '-r', help="the path with recombination results", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--spec_tab', '-s', help="the path to the table denoting host specificity", 
              type=click.Path(exists=True), metavar='<PATH>')  
@click.option('--strains_tab', '-t', help="the path to the table denoting strains attributions", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--spec_ordres', '-o', help="the path to the table denoting orders of affected species", 
              type=click.Path(exists=True), metavar='<PATH>')
#@click.option('--tree_path', '-p', help="the path to renamed trees", 
#              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--ident_stat', '-i', help="the path to mean identity comparisions", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--ref_list', '-l', help="the path to reference proteins", 
              type=click.Path(exists=True), metavar='<PATH>')

def main(rec_tab, spec_tab, strains_tab, spec_ordres, ident_stat, ref_list): #tree_path

    tox_df = pd.read_csv(os.path.realpath(spec_tab), sep='\t')
    tox_summary_dict = clean_tox_dict(make_tox_summary_dict(tox_df))
    spec_ord_df = pd.read_csv(os.path.realpath(spec_ordres), sep='\t')
    spec_ord_dict = pd.Series(spec_ord_df.Host_order.values,index=spec_ord_df.Host_sp).to_dict()

    heatmap_dict=get_identity_heatmap(read_csv_to_list(ident_stat, headless=True,delim='\t'))

    for tox_pair in heatmap_dict:
        heatmap_dict[tox_pair]['domains12'] =  (heatmap_dict[tox_pair]['domain1']+heatmap_dict[tox_pair]['domain2'])/2
        heatmap_dict[tox_pair]['domains13'] = (heatmap_dict[tox_pair]['domain1']+heatmap_dict[tox_pair]['domain3'])/2
        heatmap_dict[tox_pair]['domains23'] = (heatmap_dict[tox_pair]['domain2']+heatmap_dict[tox_pair]['domain3'])/2

    ref_res = list(set([cry for cry in unpacking_list_of_lists(read_csv_to_list(os.path.realpath(ref_list), headless=False))]))
    count_host_specificity(rec_tab, tox_summary_dict, strains_tab, spec_ord_dict, heatmap_dict, ref_res)
    

if __name__ == '__main__':
   main()





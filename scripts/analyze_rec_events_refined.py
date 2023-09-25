import sys
import csv
from collections import defaultdict
import os
import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from annotate_clusters_consistency import write_csv, read_csv_to_list, create_dict_from_list
from make_data_for_heatmap_universal import make_sequence_dict, find_sequence_identity
import pandas as pd
from count_host_changes import get_parent_names
from itertools import combinations
from create_evol_data_on_fixed_recombinants import get_identity_heatmap, make_clusters_dict
import Levenshtein as lv


def return_domain_names(rtype: str) -> list:
    """
        Returns domain names for a recombination event either for domains exposed to recombination or the intact ones
    """
    if rtype.startswith('domain1'):
        return ['domain2', 'domain3']
    elif rtype.startswith('domain2'):
         return ['domain1','domain3']
    elif rtype.startswith('domain3'):
        return ['domain1','domain2']
    else:
        return None

def find_unknown_in_dataset(line, domain , heatmap_dict):
    ignore_set=list(set(line[1]['Rec'].split(';') +get_parent_names(line[1]['parent_1'])+ get_parent_names(line[1]['parent_2'])+get_parent_names(line[1]['parent_3'])))
    ids_old=list()
    new_par='-'

    for rec in line[1]['Rec'].split(';'):
        for par in get_parent_names(line[1]['Min_par']):
            ids_old.append(heatmap_dict['|'.join(sorted([rec, par]))][domain])
        id_thr=max(ids_old)

        for key_pair in heatmap_dict:
            if rec in key_pair:
                if heatmap_dict[key_pair][domain] > id_thr and key_pair.replace(rec,'').replace('|','') not in ignore_set:
                    id_thr=heatmap_dict[key_pair][domain]
                    new_par=key_pair.replace(rec,'').replace('|','')
    if line[1]['Type']==domain:
        ret_list=[line[1]['ID'],line[1][ 'Rec'], line[1][ 'parent_' +domain[-1]] , domain, 'Min', max(ids_old),id_thr,new_par]
    elif line[1]['Type'] != domain:
        ret_list=[line[1]['ID'],line[1][ 'Rec'], line[1][ 'parent_' +domain[-1]] , domain, 'Maj', max(ids_old),id_thr,new_par]
    return(ret_list)
                         


def get_unknowns(rdp_file, ref_recs: defaultdict(list), out_dir: str, id_tab):
    """
        Reads recombination results and returns sequences of unknown domains
    """
    file_reader = read_csv_to_list(os.path.realpath(rdp_file), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    un_list=list()
    doms_list=[[],[],[]]
    doms_stat=[]
    doms_stat.append(['ID', 'Rec', 'Domain', 'Dom_type'])
    unknowns_stat=[]
    unknowns_stat.append(['ID', 'Rec','Parent', 'Domain', 'Dom_type', 'Identity', 'New_identity','New_parent'])

    heatmap_dict=get_identity_heatmap(read_csv_to_list(id_tab, headless=True,delim='\t'))

    for line in rec_df.iterrows(): 
        if 'unknown' in line[1]['parent_1']:
            doms_list[0].extend(line[1]['Rec'].split(';'))
            if line[1]['Type']=='domain1':
                doms_stat.append([line[1]['ID'],line[1][ 'Rec'], 'domain1', 'Min'])
            else:
                doms_stat.append([line[1]['ID'],line[1][ 'Rec'], 'domain1', 'Maj'])

            unknowns_stat.append(find_unknown_in_dataset(line, 'domain1' , heatmap_dict))

        elif 'unknown' in line[1]['parent_2']:
            doms_list[1].extend(line[1]['Rec'].split(';'))
            if line[1]['Type']=='domain2':
                doms_stat.append([line[1]['ID'],line[1][ 'Rec'], 'domain2', 'Min'])
            else:
                doms_stat.append([line[1]['ID'],line[1][ 'Rec'], 'domain2', 'Maj'])
            unknowns_stat.append(find_unknown_in_dataset(line, 'domain2' , heatmap_dict))

        elif 'unknown' in line[1]['parent_3']:
            doms_list[2].extend(line[1]['Rec'].split(';'))
            if line[1]['Type']=='domain3':
                doms_stat.append([line[1]['ID'],line[1][ 'Rec'], 'domain3', 'Min'])
            else:
                doms_stat.append([line[1]['ID'],line[1][ 'Rec'], 'domain3', 'Maj'])
            unknowns_stat.append(find_unknown_in_dataset(line, 'domain3' , heatmap_dict))


    #write_csv(doms_stat,os.path.join(os.path.realpath(out_dir),'unknown_domains.tsv'))
    write_csv(unknowns_stat,os.path.join(os.path.realpath(out_dir),'unknown_stats.tsv'))

    #for i in range(3):
    #    rec_list=list()
    #    for tox in set(doms_list[i]):
    #        rec_list.append(SeqRecord(ref_recs[tox][i],id=tox,description=''))
    #    SeqIO.write(rec_list,os.path.join(os.path.realpath(out_dir),'unknown_domain'+str(i+1)+'.fasta'),"fasta")

def calculate_mismatch_rate(rdp_file : str, ref_recs: defaultdict(list)):
    """
        Calculates mismatch rate for all single recombination events 
    """
    file_reader = read_csv_to_list(os.path.realpath(rdp_file), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])
    ret_rows_single=list()
    ret_rows_all=list()
    ret_rows_10=list()

    ret_rows_single.append(['ID','Tox_type','Domain', 'lev_dist', 'align_dist'])
    ret_rows_all.append(['ID','Tox_type','Domain', 'lev_dist', 'align_dist'])
    ret_rows_10.append(['ID','Tox_type','Domain', 'lev_dist', 'align_dist'])

    for line in rec_df.iterrows(): 
        #single events
        if len(get_parent_names(line[1]['parent_1']))+len(get_parent_names(line[1]['parent_2']))+len(get_parent_names(line[1]['parent_3']))+len(get_parent_names(line[1]['Rec']))==4:
            if 'unknown' not in get_parent_names(line[1]['parent_1']):
                dom1_rec=ref_recs[get_parent_names(line[1]['Rec'])[0]][0]
                dom1_par=ref_recs[get_parent_names(line[1]['parent_1'])[0]][0]
                glob_matches = find_sequence_identity((dom1_rec, dom1_par), ret_mode='matches', gap_mode='pen')
                print(glob_matches)

                if line[1]['Type']=='domain1':
                    ret_rows_single.append([line[1]['ID'],'Min', 'domain1', lv.distance(dom1_rec, dom1_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_single.append([line[1]['ID'],'Maj', 'domain1', lv.distance(dom1_rec, dom1_par),glob_matches[1]-glob_matches[0]])

            if 'unknown' not in get_parent_names(line[1]['parent_2']):
                dom2_rec=ref_recs[get_parent_names(line[1]['Rec'])[0]][1]
                dom2_par=ref_recs[get_parent_names(line[1]['parent_2'])[0]][1]
                glob_matches = find_sequence_identity((dom2_rec, dom2_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain2':
                    ret_rows_single.append([line[1]['ID'],'Min', 'domain2', lv.distance(dom2_rec, dom2_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_single.append([line[1]['ID'],'Maj', 'domain2', lv.distance(dom2_rec, dom2_par),glob_matches[1]-glob_matches[0]])

            if 'unknown' not in get_parent_names(line[1]['parent_3']):
                dom3_rec=ref_recs[get_parent_names(line[1]['Rec'])[0]][2]
                dom3_par=ref_recs[get_parent_names(line[1]['parent_3'])[0]][2]
                glob_matches = find_sequence_identity((dom3_rec, dom3_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain3':
                    ret_rows_single.append([line[1]['ID'],'Min', 'domain3', lv.distance(dom3_rec, dom3_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_single.append([line[1]['ID'],'Maj', 'domain3', lv.distance(dom3_rec, dom3_par),glob_matches[1]-glob_matches[0]])

        #<=10 events
        if len(get_parent_names(line[1]['parent_1']))+len(get_parent_names(line[1]['parent_2']))+len(get_parent_names(line[1]['parent_3']))+len(get_parent_names(line[1]['Rec']))<=10:
            if 'unknown' not in get_parent_names(line[1]['parent_1']):
                dom1_rec=ref_recs[get_parent_names(line[1]['Rec'])[0]][0]
                dom1_par=ref_recs[get_parent_names(line[1]['parent_1'])[0]][0]
                glob_matches = find_sequence_identity((dom1_rec, dom1_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain1':
                    ret_rows_10.append([line[1]['ID'],'Min', 'domain1', lv.distance(dom1_rec, dom1_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_10.append([line[1]['ID'],'Maj', 'domain1', lv.distance(dom1_rec, dom1_par),glob_matches[1]-glob_matches[0]])

            if 'unknown' not in get_parent_names(line[1]['parent_2']):
                dom2_rec=ref_recs[get_parent_names(line[1]['Rec'])[0]][1]
                dom2_par=ref_recs[get_parent_names(line[1]['parent_2'])[0]][1]
                glob_matches = find_sequence_identity((dom2_rec, dom2_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain2':
                    ret_rows_10.append([line[1]['ID'],'Min', 'domain2', lv.distance(dom2_rec, dom2_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_10.append([line[1]['ID'],'Maj', 'domain2', lv.distance(dom2_rec, dom2_par),glob_matches[1]-glob_matches[0]])

            if 'unknown' not in get_parent_names(line[1]['parent_3']):
                dom3_rec=ref_recs[get_parent_names(line[1]['Rec'])[0]][2]
                dom3_par=ref_recs[get_parent_names(line[1]['parent_3'])[0]][2]
                glob_matches = find_sequence_identity((dom3_rec, dom3_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain3':
                    ret_rows_10.append([line[1]['ID'],'Min', 'domain3', lv.distance(dom3_rec, dom3_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_10.append([line[1]['ID'],'Maj', 'domain3', lv.distance(dom3_rec, dom3_par),glob_matches[1]-glob_matches[0]])

        #all events
        if  'unknown' not in get_parent_names(line[1]['parent_1']):
            for pair1 in list(combinations(get_parent_names(line[1]['parent_1'])+get_parent_names(line[1]['Rec']),2)):

                dom1_rec=ref_recs[pair1[0]][0]
                dom1_par=ref_recs[pair1[1]][0]
                glob_matches = find_sequence_identity((dom1_rec, dom1_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain1':
                    ret_rows_all.append([line[1]['ID'],'Min', 'domain1', lv.distance(dom1_rec, dom1_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_all.append([line[1]['ID'],'Maj', 'domain1', lv.distance(dom1_rec, dom1_par),glob_matches[1]-glob_matches[0]])

        if  'unknown' not in get_parent_names(line[1]['parent_2']):
            for pair1 in list(combinations(get_parent_names(line[1]['parent_2'])+get_parent_names(line[1]['Rec']),2)):

                dom2_rec=ref_recs[pair1[0]][1]
                dom2_par=ref_recs[pair1[1]][1]
                glob_matches = find_sequence_identity((dom2_rec, dom2_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain2':
                    ret_rows_all.append([line[1]['ID'],'Min', 'domain2', lv.distance(dom2_rec, dom2_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_all.append([line[1]['ID'],'Maj', 'domain2', lv.distance(dom2_rec, dom2_par),glob_matches[1]-glob_matches[0]])


        if  'unknown' not in get_parent_names(line[1]['parent_3']):
            for pair1 in list(combinations(get_parent_names(line[1]['parent_3'])+get_parent_names(line[1]['Rec']),2)):

                dom3_rec=ref_recs[pair1[0]][2]
                dom3_par=ref_recs[pair1[1]][2]
                glob_matches = find_sequence_identity((dom3_rec, dom3_par), ret_mode='matches', gap_mode='pen')

                if line[1]['Type']=='domain3':
                    ret_rows_all.append([line[1]['ID'],'Min', 'domain3', lv.distance(dom3_rec, dom3_par),glob_matches[1]-glob_matches[0]])
                else:
                    ret_rows_all.append([line[1]['ID'],'Maj', 'domain3', lv.distance(dom3_rec, dom3_par),glob_matches[1]-glob_matches[0]])

    write_csv(ret_rows_single,os.path.join(os.path.realpath(args.out_dir),'single_mismatch_rates.tsv'))
    write_csv(ret_rows_all,os.path.join(os.path.realpath(args.out_dir),'all_mismatch_rates.tsv'))
    write_csv(ret_rows_10,os.path.join(os.path.realpath(args.out_dir),'mismatch_rates_10.tsv'))

def count_rec_events(rdp_file : str):
    file_reader = read_csv_to_list(os.path.realpath(rdp_file), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])
    rec_stat=[]
    rec_stat.append(['ID','domain', 'Num_recs', 'Num_pars','Par_type'])
    rec_per_domain_stat=[]
    rec_per_domain_stat.append(['ID','domain', 'Num_recs', 'Num_pars', 'Num_tox'])

    for line in rec_df.iterrows(): 
        parent1=get_parent_names(line[1]['parent_1'])
        parent2=get_parent_names(line[1]['parent_2'])
        parent3=get_parent_names(line[1]['parent_3'])
        parent_set=list(set(parent1+parent2+parent3))
        parent12=list(set(parent1+parent2))
        parent13=list(set(parent1+parent3))
        parent23=list(set(parent2+parent3))

        recs=get_parent_names(line[1]['Rec'])
        domain=line[1]['Type']
        ID=line[1]['ID']

        for par_set, set_name in zip([parent_set,parent1,parent2,parent3,parent1+parent2,parent1+parent3, parent2+parent3],['All','Min1','Min2','Min3','Maj12','Maj13', 'Maj23']):
            try:
                par_set=list(set(par_set)).remove('unknown')

            except:
                par_set=list(set(par_set))

            if par_set==None:
                par_set=[]

            if domain=='domain1':
                if set_name=='Min2' or set_name=='Min3' or set_name=='Maj13' or set_name=='Maj12':
                    continue

            if domain=='domain2':
                if set_name=='Min1' or set_name=='Min3' or set_name=='Maj12' or set_name=='Maj23':
                    continue

            if domain=='domain3':
                if set_name=='Min1' or set_name=='Min2' or set_name=='Maj13' or set_name=='Maj23':
                    continue

            if len(par_set)!=0:
                rec_stat.append([ID,domain,len(recs), len(par_set), set_name.replace('1','').replace('2','').replace('3','').replace('12','').replace('13','').replace('23','')])

        if 'unknown' not in parent1 and len(parent1)>0 and parent1!=None:
            rec_per_domain_stat.append([ID, 'domain1', len(recs), len(parent1), len(recs) + len(parent1)])

        if 'unknown' not in parent2 and len(parent2)>0 and parent2!=None:
            rec_per_domain_stat.append([ID, 'domain2', len(recs), len(parent2), len(recs) + len(parent2)])

        if 'unknown' not in parent3 and len(parent3)>0 and parent3!=None:
            rec_per_domain_stat.append([ID, 'domain3', len(recs), len(parent3), len(recs) + len(parent3)])
            
          
    #write_csv(rec_stat,os.path.join(os.path.realpath(args.out_dir),'number_of_toxins_recombinants.tsv'))  
    write_csv(rec_per_domain_stat,os.path.join(os.path.realpath(args.out_dir),'per_domain_number_of_toxins_recombinants.tsv'))              
            
def find_ids(tox_list1, tox_list2, heatmap_dict, domain):
    mean_ids=[]
    for rec in tox_list1:
        for par in tox_list2:
            mean_ids.append(heatmap_dict['|'.join(sorted([rec, par]))][domain])
    return(sum(mean_ids)/len(mean_ids))


def find_mismatches(tox_list1, tox_list2, heatmap_dict, domain):
    mean_ids=[]
    for rec in tox_list1:
        for par in tox_list2:
            mean_ids.append(lv.distance(rec, par))
    return(sum(mean_ids)/len(mean_ids))
    

def mismathces_and_id_by_member(rdp_file : str,  id_tab, ref_recs: defaultdict(list)):
    file_reader = read_csv_to_list(os.path.realpath(rdp_file), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])
    rec_nums_list=[]
    rec_nums_list.append(['ID','Num_recs', 'Num_pars','Num_all', 'Domain', 'Curr_domain', 'Par_type', 'Mean_ident'])
    heatmap_dict=get_identity_heatmap(read_csv_to_list(id_tab, headless=True,delim='\t'))
    
    for line in rec_df.iterrows(): 
        parent1=get_parent_names(line[1]['parent_1'])
        parent2=get_parent_names(line[1]['parent_2'])
        parent3=get_parent_names(line[1]['parent_3'])
        parent12=list(set(parent1+parent2))
        parent13=list(set(parent1+parent3))
        parent23=list(set(parent2+parent3))

        for par_set, set_name in zip([parent1,parent2,parent3,parent1+parent2,parent1+parent3, parent2+parent3],['Min1','Min2','Min3','Maj12','Maj13', 'Maj23']):
            try:
                par_set=list(set(par_set)).remove('unknown')

            except:
                par_set=list(set(par_set))

            if par_set==None:
                par_set=[]


            recs=get_parent_names(line[1]['Rec'])
            domain=line[1]['Type']
            ID=line[1]['ID']

            if domain=='domain1':

                if len(par_set)>0 and set_name=='Min1':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain1')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain1', 'Min_par', ident])
                    #dom1_recs=ref_recs[get_parent_names(line[1]['Rec'])[0]][0]
                    #dom1_pars=ref_recs[get_parent_names(line[1]['parent_1'])[0]][0]
                    #print(ref_recs[recs[0]][0])

                if len(par_set)>0 and set_name=='Min2':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain2')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain2', 'Maj_par', ident])

                if len(par_set)>0 and set_name=='Min3':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain3')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain3', 'Maj_par', ident])

            if domain=='domain2':

                if len(par_set)>0 and set_name=='Min2':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain2')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain2', 'Min_par', ident])

                if len(par_set)>0 and set_name=='Min1':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain1')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain1', 'Maj_par', ident])

                if len(par_set)>0 and set_name=='Min3':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain3')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain3', 'Maj_par', ident])

            if domain=='domain3':

                if len(par_set)>0 and set_name=='Min3':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain3')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain3', 'Min_par', ident])

                if len(par_set)>0 and set_name=='Min2':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain2')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain2', 'Maj_par', ident])

                if len(par_set)>0 and set_name=='Min1':
                    ident=find_ids(recs, par_set, heatmap_dict, 'domain1')
                    rec_nums_list.append([ID,len(recs), len(par_set),len(recs)+len(par_set), domain, 'domain1', 'Maj_par', ident])

    write_csv(rec_nums_list,os.path.join(os.path.realpath(args.out_dir),'mean_ids_for_events.tsv'))   

def disclose_cluster(toxin: str, domain: str, clust_dict: defaultdict(dict)) -> list:
    """
        For a toxin and a domain specified  returns protein names included in the cluster
        Returns no more than four toxins from a cluster
    """   

    if toxin not in clust_dict:
        return([])
    
    else:
        if domain not in clust_dict[toxin]:
            return([])
        else:
            return(clust_dict[toxin][domain])

def calculate_protein_numbers_clustr(dom_seq_dict_ref, clust_dict, rec_tab):

    all_clust=[]
    all_clust.append(['ID', 'Tox_number', 'domain'])

    file_reader = read_csv_to_list(os.path.realpath(rec_tab), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    for line in rec_df.iterrows(): 
        recs=get_parent_names(line[1]['Rec'])
        domain=line[1]['Type']
        ID=line[1]['ID']
        dom1=get_parent_names(line[1]['parent_1'])+recs
        dom2=get_parent_names(line[1]['parent_2'])+recs
        dom3=get_parent_names(line[1]['parent_3'])+recs

        dom1_list=[]
        for tox in dom1:
            if tox!='' and tox!='unknown':
                dom1_list.append(tox)
                dom1_list.extend(disclose_cluster(tox,'domain1', clust_dict))
        all_clust.append([ID, len(dom1_list), 'domain1'])

        dom2_list=[]
        for tox in dom2:
            if tox!='' and tox!='unknown':
                dom2_list.append(tox)
                dom2_list.extend(disclose_cluster(tox,'domain2', clust_dict))
        all_clust.append([ID, len(dom2_list), 'domain2'])


        dom3_list=[]
        for tox in dom3:
            if tox!='' and tox!='unknown':
                dom3_list.append(tox)
                dom3_list.extend(disclose_cluster(tox,'domain3', clust_dict))
        all_clust.append([ID, len(dom3_list), 'domain3'])

    write_csv(all_clust,os.path.join(os.path.realpath(args.out_dir),'number_of_toxins_recombinants_clusters.tsv')) 


def calculate_mean_ids_clustr(dom_seq_dict_ref, clust_dict, rec_tab):

    all_ids=[]
    all_ids.append(['ID', 'Mean_id', 'domain', 'Num_tox'])

    file_reader = read_csv_to_list(os.path.realpath(rec_tab), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    for line in rec_df.iterrows(): 
        recs=get_parent_names(line[1]['Rec'])
        domain=line[1]['Type']
        ID=line[1]['ID']
        dom1=get_parent_names(line[1]['parent_1'])+recs
        dom2=get_parent_names(line[1]['parent_2'])+recs
        dom3=get_parent_names(line[1]['parent_3'])+recs

        dom1_list=[]
        dom1_ids=[]
        for tox in dom1:
            if tox!='' and tox!='unknown':
                dom1_list.append(tox)
                dom1_list.extend(disclose_cluster(tox,'domain1', clust_dict))

        if len(dom1_list)>=2:
            for tox_pair in list(set(list(combinations(dom1_list,2)))):
                seq_tuple=((dom_seq_dict_ref[tox_pair[0]][0],dom_seq_dict_ref[tox_pair[1]][0]))
                sim = find_sequence_identity(seq_tuple)[1]
                print('domain1', tox_pair)
                dom1_ids.append(sim)

        dom2_list=[]
        dom2_ids=[]
        for tox in dom2:
            if tox!='' and tox!='unknown':
                dom2_list.append(tox)
                dom2_list.extend(disclose_cluster(tox,'domain2', clust_dict))

        if len(dom2_list)>=2:
            for tox_pair in list(set(list(combinations(dom2_list,2)))):
                seq_tuple=((dom_seq_dict_ref[tox_pair[0]][1],dom_seq_dict_ref[tox_pair[1]][1]))
                sim = find_sequence_identity(seq_tuple)[1]
                print('domain2', tox_pair)
                dom2_ids.append(sim)

        dom3_list=[]
        dom3_ids=[]
        for tox in dom3:
            if tox!='' and tox!='unknown':
                dom3_list.append(tox)
                dom3_list.extend(disclose_cluster(tox,'domain3', clust_dict))

        if len(dom3_list)>=2:
            for tox_pair in list(set(list(combinations(dom3_list,2)))):
                seq_tuple=((dom_seq_dict_ref[tox_pair[0]][2],dom_seq_dict_ref[tox_pair[1]][2]))
                sim = find_sequence_identity(seq_tuple)[1]
                print('domain3', tox_pair)
                dom3_ids.append(sim)

        if len(dom1_ids)>=1:
            all_ids.append([ID, sum(dom1_ids)/len(dom1_ids), 'domain1', len(dom1_list)])
        if len(dom2_ids)>=1:
            all_ids.append([ID, sum(dom2_ids)/len(dom2_ids), 'domain2',len(dom2_list) ])
        if len(dom3_ids)>=1:
            all_ids.append([ID, sum(dom3_ids)/len(dom3_ids), 'domain3', len(dom3_list)])

    write_csv(all_ids,os.path.join(os.path.realpath(args.out_dir),'mean_ids_for_events_clust.tsv')) 





def make_ids_and_mismathces_distribution(rdp_file : str, ref_recs: defaultdict(list)):
    file_reader = read_csv_to_list(os.path.realpath(rdp_file), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])
    #ID, n_recs, n_min, n_maj1,n_maj2, n_all, domain, unknown_flag, group
    #domain toxin2 toxin2 ID_pair len1 len2  len_aln lev_dist align_dist 

    rec_id_and_mism_res=[['ID','n_recs', 'n_min', 'n_maj1','n_maj2', 'n_all', 
                         'domain', 'unknown_flag', 'group','curr_domain', 'toxin1', 'toxin2', 
                         'ID_pair', 'len1', 'len2' , 'len_aln', 'lev_dist', 'align_dist']]
    
    for row in rec_df.iterrows(): 


        recs=get_parent_names(row[1]['Rec'])
        domain=row[1]['Type']
        ID=row[1]['ID']
        dom_comp_order=dict()
        unknown_flag='no'
          
        if row[1]['Type']=='domain1':
            min=get_parent_names(row[1]['parent_1'])
            min = [m.strip() for m in min]

            maj1=get_parent_names(row[1]['parent_2'])
            maj1 = [m.strip() for m in maj1]

            maj2=get_parent_names(row[1]['parent_3'])
            maj2 = [m.strip() for m in maj2]

            dom_comp_order['domain1']='min'
            dom_comp_order['domain2']='maj1'
            dom_comp_order['domain3']='maj2'


        if row[1]['Type']=='domain2':
            min=get_parent_names(row[1]['parent_2'])
            min = [m.strip() for m in min]

            maj1=get_parent_names(row[1]['parent_1'])
            maj1 = [m.strip() for m in maj1]

            maj2=get_parent_names(row[1]['parent_3'])
            maj2 = [m.strip() for m in maj2]

            dom_comp_order['domain1']='maj1'
            dom_comp_order['domain2']='min'
            dom_comp_order['domain3']='maj2'

        if row[1]['Type']=='domain3':
            min=get_parent_names(row[1]['parent_3'])
            min = [m.strip() for m in min]

            maj1=get_parent_names(row[1]['parent_1'])
            maj1 = [m.strip() for m in maj1]

            maj2=get_parent_names(row[1]['parent_2'])
            maj2 = [m.strip() for m in maj2]

            dom_comp_order['domain1']='maj1'
            dom_comp_order['domain2']='maj2'
            dom_comp_order['domain3']='min'

        if 'unknown' in min:
            unknown_flag=unknown_flag.replace('no', '') + 'minor'

        if 'unknown' in maj1:
            unknown_flag=unknown_flag.replace('no', '') + 'major1'

        if 'unknown' in maj2:
            unknown_flag=unknown_flag.replace('no', '') + 'major2'

        n_all=len(set(recs).union(set(maj1)).union(set(maj2)).union(set(min)))
        n_recs=len(set(recs))
        n_min=len(set(min))
        n_maj1=len(set(maj1))
        n_maj2=len(set(maj2))

        #                dom1_rec=ref_recs[pair1[0]][0]
        #        dom1_par=ref_recs[pair1[1]][0]


        def get_ids_for_group(tox_list, ref_recs):
            #domain toxin2 toxin2 ID_pair len1 len2  len_aln lev_dist align_dist 
            ret_list=[]

            if len(tox_list)==1:
                for domain in ['domain1', 'domain2', 'domain3']:
                    domain_index=int(domain.replace('domain',''))-1
                    ret_list.append([domain, tox_list[0], tox_list[0], 100, 
                                    len(ref_recs[tox_list[0]][domain_index]), len(ref_recs[tox_list[0]][domain_index]),len(ref_recs[tox_list[0]][domain_index]), 0, 0])
            else:
                for tox_pair in list(set(list(combinations(tox_list,2)))):
                    for domain in ['domain1', 'domain2', 'domain3']:
                        domain_index=int(domain.replace('domain',''))-1
                        seq1 = ref_recs[tox_pair[0]][domain_index]
                        seq2 = ref_recs[tox_pair[1]][domain_index]
                        len1 = len(seq1)
                        len2 = len(seq2)

                        seq_tuple=((seq1,seq2))
                        sim = find_sequence_identity(seq_tuple)[0]
                        
                        aln_len = find_sequence_identity(seq_tuple)[2]
                        glob_matches = find_sequence_identity(seq_tuple, ret_mode='matches', gap_mode='pen')
                        
                        ret_list.append([domain, tox_pair[0], tox_pair[1], sim, len1, len2, aln_len, glob_matches[1]-glob_matches[0], lv.distance(seq1, seq2)])
            return(ret_list)
                        

        groups_dict=defaultdict(list)
        for group in ['rec_rec', 'rec_min', 'rec_maj1', 'rec_maj2', 'min_min', 'min_maj1', 'min_maj2', 'maj1_maj2', 'maj1_maj1', 'maj2_maj2']:
            group_list=[]
            if 'rec' in group:
              group_list.extend(recs)  
            if 'min' in group:
              group_list.extend(min) 
            if 'maj1' in group:
              group_list.extend(maj1) 
            if 'maj2' in group:
              group_list.extend(maj2)
            
            groups_dict[group]=list(set(group_list))
        
        print(ID)
        for group in groups_dict:
            if 'unknown' in groups_dict[group]:
                pass
            else:
                group_comp_results = get_ids_for_group(groups_dict[group], ref_recs)
                for comp_row in group_comp_results:
                    rec_id_and_mism_res.append([ID, n_recs, n_min, n_maj1,n_maj2, n_all, domain, unknown_flag, group]+ comp_row)

    write_csv(rec_id_and_mism_res,os.path.join(os.path.realpath(args.out_dir),'Mearged_stat_id_and_mism.tsv')) 



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reads recombination results and gets unknown domains, mean mismatch rate, etc.')
    parser.add_argument('-o', '--out', dest='out_dir', help='the path to the output directory',
                        type=str)
    parser.add_argument('-r', '--rec', dest='rec_tab', help='the path to the raw RDP table',
                        type=str)
    parser.add_argument('-p', '--parsed_t', dest='parsed_trees', help='the path to the parsed tree tables',
                        type=str)
    parser.add_argument('-rn', '--ref', dest='ref_nucl', help='the path to the reference nucleotide sequences',
                        type=str)
    parser.add_argument('-i', '--id', dest='id_tab', help='the path to the table with identity',
                        type=str)
    parser.add_argument('-c', '--cd', dest='clust_dat', help='the path to the table with clusters',
                        type=str)
    parser.add_argument('-an', '--all', dest='all_nucl', help='the path to the mearged nucleotide sequences',
                        type=str)

    args = parser.parse_args()

    dom_seq_dict_ref=make_sequence_dict(args.ref_nucl)
    dom_seq_dict_all=make_sequence_dict(args.all_nucl)
    #get_unknowns(args.rec_tab,dom_seq_dict, args.out_dir,  args.id_tab)
    #calculate_mismatch_rate(args.rec_tab,dom_seq_dict_all)
    count_rec_events(args.rec_tab)
    clust_dict = make_clusters_dict(read_csv_to_list(args.clust_dat, headless=True,delim='\t'))
    mismathces_and_id_by_member(args.rec_tab, args.id_tab, dom_seq_dict_ref)
    calculate_protein_numbers_clustr(dom_seq_dict_ref, clust_dict, args.rec_tab)
    calculate_mean_ids_clustr(dom_seq_dict_all, clust_dict, args.rec_tab)
    
    make_ids_and_mismathces_distribution(args.rec_tab,dom_seq_dict_all)


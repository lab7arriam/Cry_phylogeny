import os
import click
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter

from RDP_filter_combined import clean_name
from annotate_clusters_consistency import write_csv, read_csv_to_list
from make_data_for_heatmap_universal import make_sequence_dict, find_sequence_identity
from get_breakpoints_from_alignments import fasta_file_to_dict, get_major_and_minor_parents, get_start_stop_easy


def get_parent_names(proc_str):
    try:
        rec_str = proc_str.split('(')[1].replace(')','').split(';')
    except:
        rec_str = proc_str.replace(')','').split(';')

    return([r.strip() for r in rec_str])



def get_init_RDP_dict(rdp_file : list) -> defaultdict(dict):
    """
    Reads RDP-generated table and filters events marked with swing dash
    Returns a list with parsed rows for further annotaion
    """
    event_dict=defaultdict(dict)
    tilda_names=set()
    envent_dict_filt = defaultdict(dict)
    ind_dict = {8:'recs',10:'major_par',9:'minor_par'}
    for row in rdp_file:
        if len(row)>1:

            if len(row[0])>0:

                key=row[0].replace(' ','')
 
                if '~' in row[1]:
                    tilda_names.add(key)

                if len(row)==21:

                    event_dict[key]['recs']=[clean_name(row[8])]
                    event_dict[key]['coords']=sorted([int(clean_name(x)) for x in  [row[5],row[4]]])
                    event_dict[key]['align_coords']=sorted([int(clean_name(x)) for x in  [row[2],row[3]]])
                    event_dict[key]['major_par']=[clean_name(row[10])]
                    event_dict[key]['minor_par']=[clean_name(row[9])]

                elif len(row)==11:
                    for ind in ind_dict:
                        if len(clean_name(row[ind]))>0:
                            event_dict[key][ind_dict[ind]].append(clean_name(row[ind]))


    for event in event_dict:
        if event not in tilda_names: 

            recs = ';'.join([rec for rec in event_dict[event]['recs'] if len(rec)>0])
            min_p = ';'.join([rec for rec in event_dict[event]['minor_par'] if len(rec)>0])
            maj_p = ';'.join([rec for rec in event_dict[event]['major_par'] if len(rec)>0])

            key_rec_min = recs +'|'+min_p+'|'+maj_p
            if key_rec_min in envent_dict_filt:
                key_rec_min+='_103_event'
            envent_dict_filt[key_rec_min] = event_dict[event]

    return(envent_dict_filt)


def make_dict_with_filt_events(rec_filt):
    filt_res = read_csv_to_list(os.path.realpath(rec_filt), headless=False)
    filt_df = pd.DataFrame(filt_res[1:], columns=filt_res[0])

    filt_events_dict = defaultdict(dict)
    for line in filt_df.iterrows(): 
        new_id = line[1]['Rec'] +'|'+line[1]['Min_par']+'|'+line[1]['Maj_par']

        if new_id in filt_events_dict:
            new_id+='_103_event'

        for index in ['ID', 'Rec','Min_par','Maj_par','Unknown_flag','Type','start', 'stop', 'd1_start', 'd1_stop', 'd2_start', 'd2_stop', 'd3_start', 'd3_stop','parent_1','parent_2','parent_3']:
            filt_events_dict[new_id][index]= line[1][index]

        filt_events_dict[new_id]['intersect_maj'] = get_major_and_minor_parents(line, mode='intersect')[2]
        filt_events_dict[new_id]['new_all_maj'] = get_major_and_minor_parents(line, mode='intersect')[1]
        filt_events_dict[new_id]['new_min_list'] = get_major_and_minor_parents(line, mode='intersect')[0]

    return(filt_events_dict)

def make_defaultdict_with_breaks_from_df(df_tab):

    df = pd.read_csv(os.path.realpath(df_tab), sep='\t')
    ret_dict = defaultdict(dict)

    for index, row in df.iterrows():
        ret_dict.setdefault(row['ID'], []).append(row)

    return(ret_dict)

def make_dict_with_RDP_annot_table(rdp_calc_tab):
    rdp_ret_dict = defaultdict(dict)
    rdp_df = pd.read_csv(os.path.realpath(rdp_calc_tab), sep='\t')

    for index, rdp_row in rdp_df.iterrows():
        new_id = rdp_row['Rec'] +'|'+rdp_row['Min_par']+'|'+rdp_row['Maj_par']
        if new_id in rdp_ret_dict:
            new_id+='_103_event'
        for inner_key in ['Rec', 'Min_par', 'Maj_par', 'Unknown_flag', 'Type']:
            rdp_ret_dict[new_id][inner_key] = rdp_row[inner_key]
        rdp_ret_dict[new_id]['coords'] = [rdp_row['start'], rdp_row['stop']]
    return(rdp_ret_dict)

def check_sublist_of_list(sub_list, full_list):

    sub_list=sorted(sub_list)
    full_list=sorted(full_list)
    els=[]

    for el in sub_list:
        if el in full_list:
            els.append(el)

    return(els==sub_list)



def compare_maj_lists(filt_events_dict):
     major_parents_list=[['ID','Rec','Unknown_flag','Type', 'Min_par','Num_min_par', 'Min_par_new','Num_min_par_new','Maj_par','Num_maj_par', 'Maj_par_new' , 'Num_maj_par_new','Par_all_init', 'Par_new', 'Rec_num','parent_1','parent_2','parent_3']]
     major_par_for_heatmap = [['ID','Tox_type', 'Num_tox','Domain', 'Unknown_type','tile_fac']]

     for key in filt_events_dict:

        mins_init=get_parent_names(filt_events_dict[key]['Min_par'])
        mins_new=filt_events_dict[key]['new_min_list']

        majs_true_init = get_parent_names(filt_events_dict[key]['Maj_par'])
        majs_new = filt_events_dict[key]['intersect_maj']
        majs_all_new = filt_events_dict[key]['new_all_maj']

        par_1 = len([par for par in filt_events_dict[key]['parent_1'].split(';') if par!='unknown'])
        par_2 = len([par for par in filt_events_dict[key]['parent_2'].split(';') if par!='unknown'])
        par_3 = len([par for par in filt_events_dict[key]['parent_3'].split(';') if par!='unknown'])
        num_recs = len(filt_events_dict[key]['Rec'].split(';'))

        if filt_events_dict[key]['Unknown_flag']!='minor':
            init_min_num=len(mins_init)

            if not check_sublist_of_list(mins_init,mins_new):
                print('Minors failed')

        else:
            init_min_num=0

        if filt_events_dict[key]['Unknown_flag']!='major':
            init_maj_num=len(majs_true_init)
            if not (check_sublist_of_list(majs_true_init, majs_new) or check_sublist_of_list(majs_new, majs_true_init)):
                print('Majors failed')
        else:
            init_maj_num=0

        if majs_all_new==['unknown']:
            majs_all_new= []

        if majs_new==['unknown']:
            majs_new= []

        if mins_new==['unknown']:
            mins_new= []

        new_min_num = len(mins_new)
        new_maj_num = len(majs_new)

        majors_compare_row = []

        for sub_key in ['ID','Rec','Unknown_flag','Type', 'Min_par']:
            majors_compare_row.append(filt_events_dict[key][sub_key])
        majors_compare_row.extend([init_min_num, ';'.join(mins_new), new_min_num, filt_events_dict[key]['Maj_par'], init_maj_num, ';'.join(majs_new), new_maj_num, init_min_num+init_maj_num, len(majs_all_new+mins_new)])

        for sub_key in ['parent_1','parent_2','parent_3']:
            majors_compare_row.append(filt_events_dict[key][sub_key])
        #print(majors_compare_row)
        major_parents_list.append(majors_compare_row)

        major_par_for_heatmap.append([filt_events_dict[key]['ID'], 'Recs', num_recs, filt_events_dict[key]['Type'], filt_events_dict[key]['Unknown_flag'], 'fac']) #rec
        major_par_for_heatmap.append([filt_events_dict[key]['ID'], 'Mins', new_min_num, filt_events_dict[key]['Type'], filt_events_dict[key]['Unknown_flag'], 'fac']) #minor
        major_par_for_heatmap.append([filt_events_dict[key]['ID'], 'Majs', new_maj_num, filt_events_dict[key]['Type'], filt_events_dict[key]['Unknown_flag'], 'fac'])#major
        major_par_for_heatmap.append([filt_events_dict[key]['ID'], 'Dom1', par_1, filt_events_dict[key]['Type'], filt_events_dict[key]['Unknown_flag'], 'fac'])#par1
        major_par_for_heatmap.append([filt_events_dict[key]['ID'], 'Dom2', par_2, filt_events_dict[key]['Type'], filt_events_dict[key]['Unknown_flag'], 'fac'])#par2
        major_par_for_heatmap.append([filt_events_dict[key]['ID'], 'Dom3', par_3, filt_events_dict[key]['Type'], filt_events_dict[key]['Unknown_flag'], 'fac'])#par3

     maj_min_compare_reshaped = [['ID', 'Num1', 'Num2', 'Type']]
     #4--> , 'Min_par', 'Num_min_par', 'Min_par_new','Num_min_par_new','Maj_par','Num_maj_par', 'Maj_par_new' , 'Num_maj_par_new','Par_all_init', 'Par_new'
     for row in major_parents_list[1:]:
         maj_min_compare_reshaped.append([row[0],row[5], row[7],'Min'])
         maj_min_compare_reshaped.append([row[0],row[9], row[11],'Maj'])
         maj_min_compare_reshaped.append([row[0],row[12], row[13],'All'])
     write_csv(major_parents_list, 'Comparisions_between_parents_filtered_results.csv')
     write_csv(major_par_for_heatmap, 'Num_rec_filtered_heatmap.csv')
     write_csv(maj_min_compare_reshaped, 'Comparisions_between_parents_filtered_reshaped.csv')



def transform_alignment_coords(key, RDP_init_dict, filt_events_dict):
    RDP_start_rec, RDP_stop_rec = RDP_init_dict[key]['coords']
    RDP_start_aln, RDP_stop_aln = RDP_init_dict[key]['align_coords']
    filt_start, filt_stop = int(filt_events_dict[key]['start']), int(filt_events_dict[key]['stop'])

    if RDP_start_rec==filt_start and RDP_stop_rec==filt_stop:
        return(RDP_start_aln, RDP_stop_aln)

    else:
        if RDP_start_rec==filt_stop and filt_start==1:
            RDP_stop_aln=RDP_start_aln
            RDP_start_aln=1
        else:
            RDP_start_aln=RDP_stop_aln
            RDP_stop_aln=5112
        return(RDP_start_aln, RDP_stop_aln)

def full_coords_from_RDP_alignments(tox, coords, nucl_aln, full_nucl):
    start, stop = coords
    tox_aln_seq = str(nucl_aln[tox].seq).replace('N','A')
    tox_full_seq = str(full_nucl[tox].seq).replace('N','A')
    
    extr_break_seq = tox_aln_seq[start:stop].replace('-','').replace('N','A')
    start_ind = tox_full_seq.find(extr_break_seq)
    stop_ind = start_ind+ len(extr_break_seq)

    if tox_full_seq[start_ind:stop_ind]!=extr_break_seq:
        print('{} failed to be a substring'.format(tox))
        print(tox_aln_seq)
        print(tox_full_seq)
        print(extr_break_seq)


    return((start_ind, stop_ind))
    

def get_true_coordinates_from_RDP_filt(RDP_init_dict, nucl_aln, full_nucl, filt_events_dict):
    RDP_breaks_dict_extr = defaultdict(dict)
    for key in filt_events_dict:
        RDP_breaks_dict_extr[key]['recs'] = {rec:() for rec in RDP_init_dict[key]['recs']}
        RDP_breaks_dict_extr[key]['minor_par'] = {min_par:() for min_par in RDP_init_dict[key]['minor_par']}
        RDP_breaks_dict_extr[key]['major_par'] = {maj_par:() for maj_par in RDP_init_dict[key]['major_par']}
        RDP_breaks_dict_extr[key]['all_tox'] = dict()

        rdp_coords = transform_alignment_coords(key, RDP_init_dict, filt_events_dict)

        for tox_group in RDP_breaks_dict_extr[key]:
            for tox in  RDP_breaks_dict_extr[key][tox_group]:
                adjust_coords = full_coords_from_RDP_alignments(tox,rdp_coords, nucl_aln, full_nucl)
                RDP_breaks_dict_extr[key][tox_group][tox] = adjust_coords

        all_parents = list(set(get_parent_names(filt_events_dict[key]['parent_1']) +get_parent_names(filt_events_dict[key]['parent_2']) + get_parent_names(filt_events_dict[key]['parent_3'])))
        if 'unknown' in all_parents:
            all_parents.remove('unknown')

        for parent in all_parents:
            adjust_coords = full_coords_from_RDP_alignments(parent,rdp_coords, nucl_aln, full_nucl)
            RDP_breaks_dict_extr[key]['all_tox'][parent] = adjust_coords

    return(RDP_breaks_dict_extr)

def compare_breakpoints_alignments(RDP_coords_full_adj, break_dict, filt_events_dict):
    break_compare_res = [['ID','Tox', 'Type', 'Domain','Curr_Domain', 'Aln_start', 'Aln_stop','RDP_start', 'RDP_stop', 'Dif_len', 'Dif_start', 'Dif_stop']]

    for key in filt_events_dict:
        ID = int(filt_events_dict[key]['ID'])
        recs = filt_events_dict[key]['Rec'].split(';')
        for rec in recs:
            recs_ajd_coords = RDP_coords_full_adj[key]['recs'][rec]
            recs_init_coords = break_dict[ID][0]['r_start_full'], break_dict[ID][0]['r_stop_full']
            domain = break_dict[ID][0]['Domain']
            dif_len = abs((recs_ajd_coords[1]-recs_ajd_coords[0])-(recs_init_coords[1]-recs_init_coords[0]))
            dif_start = abs(recs_ajd_coords[0]-recs_init_coords[0])
            dif_stop = abs(recs_ajd_coords[1]-recs_init_coords[1])
            
            break_compare_res.append([ID, 'Rec', rec, domain, domain, recs_init_coords[0], recs_init_coords[1], recs_ajd_coords[0], recs_ajd_coords[1],dif_len,dif_start,dif_stop])

        for break_row in break_dict[ID]:
            tox =  break_row['Parent']
            tox_type =  break_row['Parent_type']
            domain =  break_row['Domain']
            curr_domain =  break_row['Curr_Domain']
            recs_ajd_coords = RDP_coords_full_adj[key]['all_tox'][tox]
            recs_init_coords = break_row['p_start_full'], break_row['p_stop_full']
            dif_len = abs((recs_ajd_coords[1]-recs_ajd_coords[0])-(recs_init_coords[1]-recs_init_coords[0]))
            dif_start = abs(recs_ajd_coords[0]-recs_init_coords[0])
            dif_stop = abs(recs_ajd_coords[1]-recs_init_coords[1])

            break_compare_res.append([ID, tox_type, tox, domain, curr_domain, recs_init_coords[0], recs_init_coords[1], recs_ajd_coords[0], recs_ajd_coords[1],dif_len,dif_start,dif_stop])

    write_csv(break_compare_res, 'Comparisions_between_breakpoints.csv')

def make_symbols_dict(seqs_lists: str) -> defaultdict(dict):

    Counters = list(map(lambda x:Counter(x), zip(*seqs_lists)))
    per_symbol_dict = {i:dict(Counters[i]) for i in range(len(Counters))}

    return(per_symbol_dict)

def calculate_similarity_score(coord_dict: dict):
    num_sybm = sum(coord_dict[symb] for symb in coord_dict)
    try:
        max_num = max([coord_dict[symb] for symb in coord_dict if symb!='-']) 
    except:
        max_num=0
    if max_num == num_sybm-1:
        return(0)
    else:
       return(max_num/num_sybm)

def assign_type_to_coords(coord, doms_coords_aln, break_point):
    reg_types = []

    if coord <= break_point+30 and coord >= break_point-30 and coord!=break_point:
        reg_type='break_flank'
    elif coord ==  break_point:
        reg_type='breakpoint'
    elif coord > break_point+30:
        reg_type='100_nucl_right'
    elif coord < break_point-30:
        reg_type='100_nucl_left'
    else:
        print('Wrong flank', coord)

    reg_types.append(reg_type)
    
    if coord >= doms_coords_aln[0][0] and coord<=doms_coords_aln[0][1]:
        region = 'domain1'
    elif coord >= doms_coords_aln[1][0] and coord<=doms_coords_aln[1][1]:
        region = 'domain2'
    elif coord >= doms_coords_aln[2][0] and coord<=doms_coords_aln[2][1]:
        region = 'domain3'

    elif coord > doms_coords_aln[0][1] and coord < doms_coords_aln[1][0]:
        region = 'linker12'

    elif coord > doms_coords_aln[1][1] and coord < doms_coords_aln[2][0]:
        region = 'linker23'

    elif coord < doms_coords_aln[0][0]:
        region = 'N-term'

    elif coord > doms_coords_aln[2][1]:
        region = 'C-term'
    else:
        print('Wrong domain flank')

    reg_types.append(region)

    return(reg_types)


    

def make_id_distribution_for_break_points(break_point, parent_list, aln_rec, mafft_res_dict, doms_coords_aln, ID, flank_num, break_len, aln_dir, dirs):
    id_distr_results=[['ID', 'Coord','Identity', 'Flank_num', 'Reg_type', 'Reg_domain','Flank_class']]
    before = break_point-break_len
    if before<0:
        before = 0
                
    after = break_point+break_len
    if after > len(aln_rec):
        after = len(aln_rec)

    parent_seqs = []

    break_substrings = defaultdict(dict)
    parent_only_records = []
    parent_only_recs_dict=dict()

    for par_tox in parent_list:
        parent_seqs.append(str(mafft_res_dict[par_tox].seq)[before:after])
        paretn_only_rec = SeqRecord(Seq(str(mafft_res_dict[par_tox].seq)[before:after].replace('-','')),
                               id=mafft_res_dict[par_tox].id,
                               description='')
        parent_only_records.append(paretn_only_rec)
        parent_only_recs_dict[par_tox] = Seq(str(mafft_res_dict[par_tox].seq)[before:after])
        extr_seq = str(mafft_res_dict[par_tox].seq)[before:after]

        break_check_rec = SeqRecord(Seq(extr_seq),
                               id=mafft_res_dict[par_tox].id,
                               description='')

        break_check_rec.seq = Seq(str(break_check_rec.seq)[0:131].replace('-','')) #[break_point-95:break_point+95]
        break_substrings[par_tox] = break_check_rec

    break_coord_res_init_aln = dict()

    #(35, 225) - for 95
    #130 - middle

    for par_tox in parent_list:
        before_break_seq = str(break_substrings[par_tox].seq)
        parent_seq = str(parent_only_recs_dict[par_tox])
        parent_no_gaps_seq = parent_seq.replace('-','')

        break_coords = get_start_stop_easy(before_break_seq,parent_no_gaps_seq , parent_seq)
        break_coord_res_init_aln[par_tox]=break_coords

    #for par_tox in break_coord_res_init_aln:

        #break_coords = break_coord_res_init_aln[par_tox]
        #break_point = (break_coords[1]+break_coords[0])/2
        #start, stop = break_coords
        #if break_point!=130.0:
        #    print(parent_only_recs_dict[par_tox])
        #    print(break_coords)
        #    if start<35:
        #        start=35
        #    if stop<225:
        #        start=225
        #    print(str(parent_only_recs_dict[par_tox][break_coords[0]:break_coords[1]]), break_point, str(parent_only_recs_dict[par_tox][start:stop]))

    aln_name = os.path.join(os.path.realpath(aln_dir),str(ID)+'_par_seqs.fasta')
    new_aln_name = aln_name.replace('_par_seqs.fasta','_par_aligned.fasta')
    SeqIO.write(parent_only_records, aln_name, "fasta")

    mafft_call = subprocess.call("mafft --quiet --maxiterate 1000 --localpair --op 7 --ep 3.5  --thread 8  {0} > {1} ".format(aln_name,new_aln_name ), shell = True)
    data_clean_call = subprocess.call("mv {0} {1} ".format(aln_name, dirs[3]), shell = True)

    mafft_res_dict_pars = fasta_file_to_dict(new_aln_name)

    coords_compare = []
    breaks_for_ID = dict()
    extr_par_seqs = []

    for par_tox in mafft_res_dict_pars:
        old_coords = break_coord_res_init_aln[par_tox]
        new_coords =  get_start_stop_easy(str(break_substrings[par_tox].seq),str(mafft_res_dict_pars[par_tox].seq).replace('-','') , str(mafft_res_dict_pars[par_tox].seq))

        coords_compare.append([ID, par_tox, old_coords[0],old_coords[1],new_coords[0],new_coords[1],new_coords==old_coords])

        old_break_point = old_coords[1]
        new_break_point = new_coords[1]

        old_break = parent_only_recs_dict[par_tox][old_break_point-10:old_break_point+11]
        new_break = mafft_res_dict_pars[par_tox].seq[new_break_point-10:new_break_point+11]

        old_break_nuc = parent_only_recs_dict[par_tox][old_break_point-3:old_break_point+4]
        new_break_nuc = mafft_res_dict_pars[par_tox].seq[new_break_point-3:new_break_point+4] #new_break_point-3:new_break_point+4]

        breaks_for_ID[par_tox]=new_coords

        extr_par_seqs.append(str(mafft_res_dict_pars[par_tox].seq[new_break_point-100:new_break_point+100]))
        #if True:
        #    print(ID)
        #    print('SUBSTRING', break_substrings[par_tox].seq)
        #    print('OLD_SUBSTRING', parent_only_recs_dict[par_tox][old_coords[0]:old_coords[1]])
        #    print('NEW_SUBSTRING', mafft_res_dict_pars[par_tox][new_coords[0]:new_coords[1]].seq)    
        #    print(old_coords, new_coords, old_break, new_break, old_break_nuc,new_break_nuc,old_break_point, new_break_point, len(parent_only_recs_dict[par_tox]), len(mafft_res_dict_pars[par_tox].seq))


    if len(set(breaks_for_ID.values()))>1:
        return((0,0))

        for par_tox in breaks_for_ID:
            old_coords = break_coord_res_init_aln[par_tox]
            old_seq = parent_only_recs_dict[par_tox][old_coords[0]:old_coords[1]]
            
            new_coords = breaks_for_ID[par_tox]
            new_seq = mafft_res_dict_pars[par_tox].seq[new_coords[0]:new_coords[1]]

            print('>' +str(ID) +' '+par_tox +' OLD')
            print(old_seq)

            print('>' +str(ID) +' '+par_tox+' NEW')
            print(new_seq)

    per_symbol_dict = make_symbols_dict(parent_seqs)
    per_sim_score = {}

    new_coord_dict = {}
    new_ind = 0

    for old_coord in range(before,after):
        new_coord_dict[new_ind] = old_coord
        new_ind+=1

    for coord in per_symbol_dict:
        sim_score = calculate_similarity_score(per_symbol_dict[coord])
        per_sim_score[coord] = sim_score


    per_symbol_dict_pars_only = make_symbols_dict(extr_par_seqs)
    per_sim_score_pars_only = {}

    for coord in per_symbol_dict_pars_only:
        sim_score = calculate_similarity_score(per_symbol_dict_pars_only[coord])
        per_sim_score_pars_only[coord] = sim_score


    #for new_coord in per_sim_score:
    #    old_coord = new_coord_dict[new_coord]
    #    coord_types = assign_type_to_coords(old_coord, doms_coords_aln, break_point)
    #    id_distr_results.append([ID, new_coord, per_sim_score[new_coord], flank_num] + coord_types)


    for new_coord in per_sim_score_pars_only:
        old_coord = new_coord_dict[new_coord+30]
        #print(new_coord, per_sim_score_pars_only[new_coord ])
        #print(new_coord+30, per_sim_score[new_coord+30])

        coord_types = assign_type_to_coords(old_coord, doms_coords_aln, break_point)
        id_distr_results.append([ID, new_coord, per_sim_score_pars_only[new_coord], flank_num] + coord_types)

    regs_classes = set([dist_res[5] for dist_res in id_distr_results[1:]])

    if flank_num=='Break1':
        if regs_classes=={'domain2', 'linker23', 'domain3'} or regs_classes=={'domain3'} or regs_classes=={'linker23', 'domain3'}:
            for iter_count in range(1,len(id_distr_results)):
                id_distr_results[iter_count].append('Zone3')
        elif regs_classes=={'domain1'} or regs_classes=={'domain1', 'N-term'}:
            for iter_count in range(1,len(id_distr_results)):
                id_distr_results[iter_count].append('Zone1')

        elif regs_classes=={'domain2'} or regs_classes=={'domain2', 'domain1', 'linker12'} or regs_classes=={'domain2', 'linker12'} :
            for iter_count in range(1,len(id_distr_results)):
                id_distr_results[iter_count].append('Zone2')

    elif flank_num=='Break2':
        if regs_classes=={'C-term', 'domain3'} or regs_classes=={'domain3'}:
            for iter_count in range(1,len(id_distr_results)):
                id_distr_results[iter_count].append('Zone4')


        elif regs_classes=={'domain1'} or regs_classes=={'domain1', 'linker12'} or regs_classes=={'domain2', 'domain1', 'linker12'} or regs_classes=={'domain2'}:
            for iter_count in range(1,len(id_distr_results)):
                id_distr_results[iter_count].append('Zone2')

        elif regs_classes=={'domain2'} or regs_classes=={'domain2', 'linker23', 'domain3'} :
            for iter_count in range(1,len(id_distr_results)):
                id_distr_results[iter_count].append('Zone3')

    return(id_distr_results, coords_compare)

def get_flanking_regions_coords(RDP_coords_full_adj, filt_events_dict, full_nucl_fasta, break_dict, aln_dir):

    dirs = []
    for index in ['raw_fastas', 'msa_all', 'ident_stat', 'par_fastas','msa_pars', 'ident_pars_only']:
        dir_name = os.path.join(os.path.realpath(aln_dir),index)
        dir_call = subprocess.call("mkdir {}".format(dir_name), shell = True)
        dirs.append(dir_name)

    id_distr_merged=[['ID', 'Coord','Identity', 'Flank_num', 'Reg_type', 'Reg_domain','Flank_class']]
    coords_compare = [['ID', 'Tox', 'Old_start', 'Old_stop','New_start', 'New_stop', 'Equal_flag']]
    for key in filt_events_dict:
        recs = filt_events_dict[key]['Rec'].split(';')
        ref_rec = recs[0]
        mins = filt_events_dict[key]['new_min_list']
        majs = filt_events_dict[key]['intersect_maj']
        ID = int(filt_events_dict[key]['ID'])

        if 'unknown' not in mins and 'unknown' not in majs and len(majs)!=0:
            seqs = []
            for tox in recs+mins+majs:
                record = full_nucl_fasta[tox]
                record.seq = Seq(str(record.seq).upper())
                seqs.append(full_nucl_fasta[tox])

            
            aln_name = os.path.join(os.path.realpath(aln_dir),str(ID)+'_all_seqs.fasta')
            new_aln_name = aln_name.replace('_all_seqs.fasta','_all_aligned.fasta')
            SeqIO.write(seqs, aln_name, "fasta")

            mafft_call = subprocess.call("mafft --quiet --maxiterate 1000 --localpair --op 7 --ep 3.5  --thread 8  {0} > {1} ".format(aln_name,new_aln_name ), shell = True)
            data_clean_call = subprocess.call("mv {0} {1} ".format(aln_name, dirs[0]), shell = True)

            mafft_res_dict = fasta_file_to_dict(new_aln_name)
            RDP_breaks = RDP_coords_full_adj[key]['recs'][ref_rec]

            aln_rec = str(mafft_res_dict[ref_rec].seq)
            rec_full = str(full_nucl_fasta[ref_rec].seq)
            rec_extr = rec_full[RDP_breaks[0]:RDP_breaks[1]]

            flank1_name = os.path.join(os.path.realpath(aln_dir),str(ID)+'_flank1.csv')
            flank2_name = os.path.join(os.path.realpath(aln_dir),str(ID)+'_flank2.csv')

            doms_coords_full = []
            for dom_ind in ['d1','d2','d3']:
                doms_coords_full.append((break_dict[ID][0][dom_ind+'_r_start_full'], break_dict[ID][0][dom_ind+'_r_stop_full']))

            doms_coords_aln=[]
            for dom_inds in doms_coords_full:
                dom_extr = rec_full[dom_inds[0]:dom_inds[1]]
                doms_coords_aln.append(get_start_stop_easy(dom_extr, rec_full, aln_rec))
                
            mafft_start, mafft_stop = get_start_stop_easy(rec_extr, rec_full, aln_rec)


            if abs(mafft_start-doms_coords_aln[0][0])>15:
               break_point1 = mafft_start
            else:
               break_point1 = None

            if abs(mafft_stop-doms_coords_aln[2][1])>15:
               break_point2 = mafft_stop
            else:
               break_point2 = None

            if break_point1:
                flank1_res, coords_compare_row = make_id_distribution_for_break_points(break_point1, mins+majs, aln_rec, mafft_res_dict, doms_coords_aln, ID, 'Break1', 130, aln_dir, dirs)
                if flank1_res!=0:
                
                    coords_compare.extend(coords_compare_row)
                    #write_csv(flank1_res, flank1_name)
                    #data_clean_call = subprocess.call("mv {0} {1} ".format(flank1_name, dirs[2]), shell = True)
                    id_distr_merged.extend(flank1_res[1:])

            if break_point2:
                flank2_res, coords_compare_row = make_id_distribution_for_break_points(break_point2, mins+majs, aln_rec, mafft_res_dict, doms_coords_aln, ID, 'Break2', 130, aln_dir, dirs)
                if flank2_res!=0:
                    coords_compare.extend(coords_compare_row)
                    #write_csv(flank2_res, flank2_name)
                    #data_clean_call = subprocess.call("mv {0} {1} ".format(flank2_name, dirs[2]), shell = True)
                    id_distr_merged.extend(flank2_res[1:])

            data_clean_call = subprocess.call("mv {0} {1} ".format(new_aln_name, dirs[1]), shell = True)

    write_csv(id_distr_merged, os.path.join(os.path.realpath(aln_dir),'Parents_only_Distribution_of_flank_identities_merged.csv'))
    write_csv(coords_compare, 'Parents_only_Breakpoints_coords_par_compare.csv')

def generate_for_for_num_comparisions(unknown_flag, num_recs, num_mins, num_majs, dom_type, filt_flag):
    if unknown_flag=='no':
        if 'partial' not in dom_type:
            compare_results = [num_recs, num_mins, num_majs,dom_type, unknown_flag, filt_flag]
        else:
            compare_results = [num_recs, num_mins, num_majs,dom_type.replace('_partial',''), unknown_flag, 'raw_data']

    elif unknown_flag=='minor':
        num_mins=0
        if 'partial' not in dom_type:
             compare_results = [num_recs, num_mins, num_majs,dom_type, unknown_flag, filt_flag]
        else:
             compare_results = [num_recs, num_mins, num_majs,dom_type.replace('_partial',''), unknown_flag, 'raw_data']

    elif unknown_flag=='major':
        num_majs=0
        if 'partial' not in dom_type:
            compare_results = [num_recs, num_mins, num_majs,dom_type, unknown_flag, filt_flag]
        else:
            compare_results = [num_recs, num_mins, num_majs,dom_type.replace('_partial',''), unknown_flag, 'raw_data']
    return(compare_results)

def compare_num_of_participants_in_events(RDP_annot_dict, filt_events_dict):
    num_compare_res = [['Num_recs','Num_mins', 'Num_majs', 'Domain', 'Unknown_type', 'Dataset']]

    ind_dict={0:'Recs', 1:'Mins', 2:'Majs'}
    events_heatmap_data = [['ID', 'Num_tox', 'Tox_type', 'Domain', 'Unknown_type', 'Dataset']]
    #Get events IDs
    ind_events = 0
    for key in RDP_annot_dict:
        RDP_annot_dict[key]['ID'] = str(ind_events)
        ind_events+=1

    for key in RDP_annot_dict:

        unknown_flag = RDP_annot_dict[key]['Unknown_flag']
        num_recs = len(RDP_annot_dict[key]['Rec'].split(';'))
        num_mins = len(RDP_annot_dict[key]['Min_par'].split(';'))
        num_majs = len(RDP_annot_dict[key]['Maj_par'].split(';'))
        dom_type = RDP_annot_dict[key]['Type']
        filt_flag = 'filt_partials'
        add_row = generate_for_for_num_comparisions(unknown_flag, num_recs, num_mins, num_majs, dom_type, filt_flag)
        num_compare_res.append(add_row)

        for ind in range(3):
            events_heatmap_data.append([RDP_annot_dict[key]['ID'], add_row[ind], ind_dict[ind]] + add_row[3:])

        if add_row[5]==filt_flag:
            new_add_row=add_row[0:5]+['raw_data']
            num_compare_res.append(new_add_row)
            for ind in range(3):
                events_heatmap_data.append([RDP_annot_dict[key]['ID'], new_add_row[ind], ind_dict[ind]] + new_add_row[3:])
        

    filt_recs_per_parents = ['ID','']
    for key in filt_events_dict:
        unknown_flag = filt_events_dict[key]['Unknown_flag']
        num_recs = len(filt_events_dict[key]['Rec'].split(';'))
        num_mins = len([el for el in filt_events_dict[key]['new_min_list'] if el!='unknown'])
        num_majs = len([el for el in filt_events_dict[key]['intersect_maj'] if el!='unknown'])
        dom_type = filt_events_dict[key]['Type']
        filt_flag = 'filt_tree'
        #if num_mins == 0 and num_majs==0:
        #    print(filt_events_dict[key])
        #    continue
        add_row = generate_for_for_num_comparisions(unknown_flag, num_recs, num_mins, num_majs, dom_type, filt_flag)
        num_compare_res.append(add_row)

        for ind in range(3):
            events_heatmap_data.append([RDP_annot_dict[key]['ID'], add_row[ind], ind_dict[ind]] + add_row[3:])

    compare_res_reshaped = [['Num','Tox_type', 'Domain', 'Unknown_type', 'Dataset']]
    for row in num_compare_res[1:]:
        compare_res_reshaped.append([row[0],'Rec']+row[3:])
        compare_res_reshaped.append([row[1],'Min']+row[3:])
        compare_res_reshaped.append([row[2],'Maj']+row[3:])


    write_csv(num_compare_res, 'Comparisions_between_number_of_recs_per_filtration.csv')
    write_csv(compare_res_reshaped, 'Comparisions_between_number_of_recs_per_filtration_reshaped.csv')
    write_csv(events_heatmap_data, 'Num_recs_per_event_for_heatmap.csv')

def transform_RDP_coords_from_annot(RDP_init_dict, RDP_annot_dict, key):

    RDP_start_rec, RDP_stop_rec = RDP_init_dict[key]['coords']
    RDP_start_aln, RDP_stop_aln = RDP_init_dict[key]['align_coords']
    annot_start, annot_stop = int(RDP_annot_dict[key]['coords'][0]), int(RDP_annot_dict[key]['coords'][1])

    if RDP_start_rec==annot_start and RDP_stop_rec==annot_stop:
        return(RDP_start_aln, RDP_stop_aln)

    else:
        if RDP_start_rec==annot_stop and annot_start==1:
            RDP_stop_aln=RDP_start_aln
            RDP_start_aln=1
        else:
            RDP_start_aln=RDP_stop_aln
            RDP_stop_aln=5112
        return(RDP_start_aln, RDP_stop_aln)

def return_domain_names(domain) -> list:
    if domain.startswith('domain1'):
        return ['domain2', 'domain3']
    elif domain.startswith('domain2'):
        return ['domain1','domain3']
    elif domain.startswith('domain3'):
        return ['domain1','domain2']
    else:
        return None


def compare_identities_with_unknowns(RDP_init_dict, RDP_annot_dict, filt_events_dict, align_fasta, nucl_dict):
    ind_events = 0
    id_results = [['ID', 'Domain', 'Curr_domain', 'Rec', 'Par', 'Par_type','Dom_type','Aln_type', 'Dataset','Ident', 'Unknown_flag']]

    for key in RDP_annot_dict:
        RDP_annot_dict[key]['ID'] = str(ind_events)
        RDP_init_dict[key]['ID'] = str(ind_events)
        ind_events+=1
        new_aln_start, new_aln_stop = transform_RDP_coords_from_annot(RDP_init_dict, RDP_annot_dict, key)
        RDP_init_dict[key]['new_alns_coords'] = [new_aln_start, new_aln_stop]

    for key in filt_events_dict:
        ID = filt_events_dict[key]['ID']
        recs = filt_events_dict[key]['Rec'].split(';')
        mins = filt_events_dict[key]['new_min_list']
        majs = [maj for maj in filt_events_dict[key]['intersect_maj'] if maj!='unknown']
        unknown = RDP_annot_dict[key]['Unknown_flag']

        domain = filt_events_dict[key]['Type']
        dom_ind = int(domain.replace('_partial','').replace('domain', ''))-1
        parent_domains = return_domain_names(domain)
        parent_domains_inds = [int(domain.replace('_partial','').replace('domain', ''))-1 for domain in parent_domains] 



        #if unknown =='no':
        #    continue

        if len(mins)==0:
            mins = ['unknown']

        if len(majs)==0:
            majs = ['unknown']

        if majs == ['unknown'] and mins == ['unknown']:
            continue

        if majs == ['unknown']:
            unknown = 'major'

        print(ID, unknown, majs, mins)

        adj_aln_start, adj_aln_stop = map(int, RDP_init_dict[key]['new_alns_coords'])

        print('STARTING MINOR TEST')
        for tox_tuple in list(product(recs, mins)):
            if unknown=='minor':
                id_results.append([ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Extr', 'filt_tree', 0, unknown])
                id_results.append([ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Full_dom', 'filt_tree', 0, unknown])
                print('MINOR PASSED')
                for dom_iter in [0,1]:
                    curr_dom = parent_domains[dom_iter]
                    curr_dom_ind = parent_domains_inds[dom_iter]

                    id_results.append([ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Extr', 'filt_tree', 0, unknown])
                    id_results.append([ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Full_dom', 'filt_tree', 0, unknown])

                break
 
            print(tox_tuple[0])
            tox1_dom = str(nucl_dict[tox_tuple[0]][dom_ind])
            tox2_dom = str(nucl_dict[tox_tuple[1]][dom_ind])

            rec_extr, min_extr = (str(align_fasta[tox_tuple[0]].seq)[adj_aln_start:adj_aln_stop+1].replace('-',''), str(align_fasta[tox_tuple[1]].seq)[adj_aln_start:adj_aln_stop+1].replace('-',''))

            ident_extr = find_sequence_identity(( rec_extr, min_extr))[1]
            ident_dom = find_sequence_identity(( tox1_dom, tox2_dom))[1]

            id_res_row_extr = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Extr', 'filt_tree', ident_extr]
            id_res_row = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Full_dom', 'filt_tree', ident_dom]
            id_results.append(id_res_row_extr)
            id_results.append(id_res_row)


            for dom_iter in [0,1]:
                curr_dom = parent_domains[dom_iter]
                curr_dom_ind = parent_domains_inds[dom_iter]
                tox1_dom = nucl_dict[tox_tuple[0]][curr_dom_ind]
                tox2_dom = nucl_dict[tox_tuple[1]][curr_dom_ind]

                ident_dom = find_sequence_identity(( tox1_dom, tox2_dom))[1]
                id_res_row_extr = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Extr', 'filt_tree', ident_dom]
                id_res_row = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Full_dom', 'filt_tree', ident_dom]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

        for tox_tuple in list(product(recs, majs)):
            print('MAJOR OUTER TEST')

            if unknown=='major':
                print('MAJORS PASSED')
                id_results.append([ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Extr', 'filt_tree', 0, unknown])
                id_results.append([ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Full_dom', 'filt_tree', 0, unknown])

                for dom_iter in [0,1]:
                    curr_dom = parent_domains[dom_iter]
                    curr_dom_ind = parent_domains_inds[dom_iter]
                    id_results.append([ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Extr', 'filt_tree', 0, unknown])
                    id_results.append([ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Full_dom', 'filt_tree', 0, unknown])

                break

            tox1_dom = str(nucl_dict[tox_tuple[0]][dom_ind])
            tox2_dom = str(nucl_dict[tox_tuple[1]][dom_ind])
            ident_par = find_sequence_identity(( tox1_dom, tox2_dom))[1]

            id_res_row_extr = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Extr', 'filt_tree', ident_par]
            id_res_row = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Full_dom', 'filt_tree', ident_par]
            id_results.append(id_res_row_extr)
            id_results.append(id_res_row)

            for dom_iter in [0,1]:
                curr_dom = parent_domains[dom_iter]
                curr_dom_ind = parent_domains_inds[dom_iter]
                tox1_dom = nucl_dict[tox_tuple[0]][curr_dom_ind]
                tox2_dom = nucl_dict[tox_tuple[1]][curr_dom_ind]
                ident_par = find_sequence_identity(( tox1_dom, tox2_dom))[1]

                id_res_row_extr = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Extr', 'filt_tree', ident_par]
                id_res_row = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Full_dom', 'filt_tree', ident_par]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

    write_csv(id_results, 'Identity_per_dataset_with_unknowns.csv')


def compare_identities_for_datasets(RDP_init_dict, RDP_annot_dict, filt_events_dict, align_fasta, nucl_dict):
    ind_events = 0
    id_results = [['ID', 'Domain', 'Curr_domain', 'Rec', 'Par', 'Par_type','Dom_type','Aln_type', 'Dataset','Ident']]
    mismatches_res = [['ID', 'Domain', 'Curr_domain', 'Rec', 'Par', 'Par_type','Dom_type','Num_mism', 'Max_length','Aln_length']]

    for key in RDP_annot_dict:
        RDP_annot_dict[key]['ID'] = str(ind_events)
        RDP_init_dict[key]['ID'] = str(ind_events)
        ind_events+=1
        new_aln_start, new_aln_stop = transform_RDP_coords_from_annot(RDP_init_dict, RDP_annot_dict, key)
        RDP_init_dict[key]['new_alns_coords'] = [new_aln_start, new_aln_stop]

    print('Filtered events')

    for key in filt_events_dict:
        ID = filt_events_dict[key]['ID']
        recs = filt_events_dict[key]['Rec'].split(';')
        mins = filt_events_dict[key]['new_min_list']
        majs = [maj for maj in filt_events_dict[key]['intersect_maj'] if maj!='unknown']
        unknown = RDP_annot_dict[key]['Unknown_flag']

        domain = filt_events_dict[key]['Type']
        dom_ind = int(domain.replace('_partial','').replace('domain', ''))-1
        parent_domains = return_domain_names(domain)
        parent_domains_inds = [int(domain.replace('_partial','').replace('domain', ''))-1 for domain in parent_domains] 

        print(ID)
        if unknown!='no' or len(mins)==0 or len(majs)==0:
            continue
        
        adj_aln_start, adj_aln_stop = map(int, RDP_init_dict[key]['new_alns_coords'])

        for tox_tuple in list(product(recs, mins)):
            tox1_dom = str(nucl_dict[tox_tuple[0]][dom_ind])
            tox2_dom = str(nucl_dict[tox_tuple[1]][dom_ind])

            rec_extr, min_extr = (str(align_fasta[tox_tuple[0]].seq)[adj_aln_start:adj_aln_stop+1].replace('-',''), str(align_fasta[tox_tuple[1]].seq)[adj_aln_start:adj_aln_stop+1].replace('-',''))

            ident_extr = find_sequence_identity(( rec_extr, min_extr))[1]
            ident_dom = find_sequence_identity(( tox1_dom, tox2_dom))[1]

            id_res_row_extr = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Extr', 'filt_tree', ident_extr]
            id_res_row = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Full_dom', 'filt_tree', ident_dom]
            id_results.append(id_res_row_extr)
            id_results.append(id_res_row)

            mism_res = find_sequence_identity((tox1_dom, tox2_dom), ret_mode='mathes_with_aln', gap_mode='pen')
            mism_row = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred',mism_res[1]-mism_res[0],mism_res[1],mism_res[2]]
            mismatches_res.append(mism_row)

            for dom_iter in [0,1]:
                curr_dom = parent_domains[dom_iter]
                curr_dom_ind = parent_domains_inds[dom_iter]
                tox1_dom = nucl_dict[tox_tuple[0]][curr_dom_ind]
                tox2_dom = nucl_dict[tox_tuple[1]][curr_dom_ind]

                ident_dom = find_sequence_identity(( tox1_dom, tox2_dom))[1]
                id_res_row_extr = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Extr', 'filt_tree', ident_dom]
                id_res_row = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Full_dom', 'filt_tree', ident_dom]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

                mism_res = find_sequence_identity((tox1_dom, tox2_dom), ret_mode='mathes_with_aln', gap_mode='pen')
                mism_row = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred',mism_res[1]-mism_res[0],mism_res[1],mism_res[2]]
                mismatches_res.append(mism_row)


        for tox_tuple in list(product(recs, majs)):
            tox1_dom = str(nucl_dict[tox_tuple[0]][dom_ind])
            tox2_dom = str(nucl_dict[tox_tuple[1]][dom_ind])
            ident_par = find_sequence_identity(( tox1_dom, tox2_dom))[1]

            id_res_row_extr = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Extr', 'filt_tree', ident_par]
            id_res_row = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Full_dom', 'filt_tree', ident_par]
            id_results.append(id_res_row_extr)
            id_results.append(id_res_row)

            mism_res = find_sequence_identity((tox1_dom, tox2_dom), ret_mode='mathes_with_aln', gap_mode='pen')
            mism_row = [ID, domain, domain, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred',mism_res[1]-mism_res[0],mism_res[1],mism_res[2]]
            mismatches_res.append(mism_row)

            for dom_iter in [0,1]:
                curr_dom = parent_domains[dom_iter]
                curr_dom_ind = parent_domains_inds[dom_iter]
                tox1_dom = nucl_dict[tox_tuple[0]][curr_dom_ind]
                tox2_dom = nucl_dict[tox_tuple[1]][curr_dom_ind]
                ident_par = find_sequence_identity(( tox1_dom, tox2_dom))[1]

                id_res_row_extr = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Extr', 'filt_tree', ident_par]
                id_res_row = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Full_dom', 'filt_tree', ident_par]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

                mism_res = find_sequence_identity((tox1_dom, tox2_dom), ret_mode='mathes_with_aln', gap_mode='pen')
                mism_row = [ID, domain, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred',mism_res[1]-mism_res[0],mism_res[1],mism_res[2]]
                mismatches_res.append(mism_row)



    print('All events')
    for key in RDP_annot_dict:
        ID = RDP_annot_dict[key]['ID']
        recs = RDP_annot_dict[key]['Rec'].split(';')
        mins = RDP_annot_dict[key]['Min_par'].split(';')
        majs = RDP_annot_dict[key]['Maj_par'].split(';')
        domain = RDP_annot_dict[key]['Type']
        domain_clear = RDP_annot_dict[key]['Type'].replace('_partial','')
        unknown = RDP_annot_dict[key]['Unknown_flag']
        dom_ind = int(domain.replace('_partial','').replace('domain', ''))-1
        parent_domains = return_domain_names(domain)
        parent_domains_inds = [int(domain.replace('_partial','').replace('domain', ''))-1 for domain in parent_domains] 

        if unknown!='no':
            continue

        print(ID)

        adj_aln_start, adj_aln_stop = map(int, RDP_init_dict[key]['new_alns_coords'])

        for tox_tuple in list(product(recs, mins)):
            tox1_dom = str(nucl_dict[tox_tuple[0]][dom_ind])
            tox2_dom = str(nucl_dict[tox_tuple[1]][dom_ind])

            rec_extr, min_extr = (str(align_fasta[tox_tuple[0]].seq)[adj_aln_start:adj_aln_stop+1].replace('-',''), str(align_fasta[tox_tuple[1]].seq)[adj_aln_start:adj_aln_stop+1].replace('-',''))

            ident_extr = find_sequence_identity(( rec_extr, min_extr))[1]
            ident_dom = find_sequence_identity(( tox1_dom, tox2_dom))[1]

            if 'partial' in domain:
                id_res_row_extr = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Extr', 'Raw', ident_extr]
                id_res_row = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Full_dom', 'Raw', ident_dom]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)
            else:
                id_res_row_extr = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Extr', 'filt_partials', ident_extr]
                id_res_row = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Full_dom', 'filt_partials', ident_dom]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

                id_res_row_extr = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Extr', 'Raw', ident_extr]
                id_res_row = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Minor', 'Transferred','Full_dom', 'Raw', ident_dom]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)


            for dom_iter in [0,1]:
                curr_dom = parent_domains[dom_iter]
                curr_dom_ind = parent_domains_inds[dom_iter]
                tox1_dom = nucl_dict[tox_tuple[0]][curr_dom_ind]
                tox2_dom = nucl_dict[tox_tuple[1]][curr_dom_ind]

                ident_dom = find_sequence_identity(( tox1_dom, tox2_dom))[1]

                if 'partial' in domain:
                    id_res_row_extr = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Extr', 'Raw', ident_dom]
                    id_res_row = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Full_dom', 'Raw', ident_dom]
                    id_results.append(id_res_row_extr)
                    id_results.append(id_res_row)
                else:
                    id_res_row_extr = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Extr', 'filt_partials', ident_dom]
                    id_res_row = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Full_dom', 'filt_partials', ident_dom]
                    id_results.append(id_res_row_extr)
                    id_results.append(id_res_row)

                    id_res_row_extr = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Extr', 'Raw', ident_dom]
                    id_res_row = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Minor', 'Non_transferred','Full_dom', 'Raw', ident_dom]
                    id_results.append(id_res_row_extr)
                    id_results.append(id_res_row)

        for tox_tuple in list(product(recs, majs)):

            tox1_dom = str(nucl_dict[tox_tuple[0]][dom_ind])
            tox2_dom = str(nucl_dict[tox_tuple[1]][dom_ind])
            ident_par = find_sequence_identity(( tox1_dom, tox2_dom))[1]
            #print(tox_tuple)


            if 'partial' in domain:
                id_res_row_extr = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Extr', 'Raw', ident_par]
                id_res_row = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Full_dom', 'Raw', ident_par]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)
            else:
                id_res_row_extr = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Extr', 'filt_partials', ident_par]
                id_res_row = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Full_dom', 'filt_partials', ident_par]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

                id_res_row_extr = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Extr', 'Raw', ident_par]
                id_res_row = [ID, domain_clear, domain_clear, tox_tuple[0], tox_tuple[1], 'Major', 'Non_transferred','Full_dom', 'Raw', ident_par]
                id_results.append(id_res_row_extr)
                id_results.append(id_res_row)

            for dom_iter in [0,1]:
                curr_dom = parent_domains[dom_iter]
                curr_dom_ind = parent_domains_inds[dom_iter]
                tox1_dom = nucl_dict[tox_tuple[0]][curr_dom_ind]
                tox2_dom = nucl_dict[tox_tuple[1]][curr_dom_ind]
                ident_par = find_sequence_identity(( tox1_dom, tox2_dom))[1]

                if 'partial' in domain:
                    id_res_row_extr = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Extr', 'Raw', ident_par]
                    id_res_row = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Full_dom', 'Raw', ident_par]
                    id_results.append(id_res_row_extr)
                    id_results.append(id_res_row)
                else:
                    id_res_row_extr = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Extr', 'filt_partials', ident_par]
                    id_res_row = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Full_dom', 'filt_partials', ident_par]
                    id_results.append(id_res_row_extr)
                    id_results.append(id_res_row)

                    id_res_row_extr = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Extr', 'Raw', ident_par]
                    id_res_row = [ID, domain_clear, curr_dom, tox_tuple[0], tox_tuple[1], 'Major', 'Transferred','Full_dom', 'Raw', ident_par]
                    id_results.append(id_res_row_extr)
                    id_results.append(id_res_row)


    write_csv(id_results, 'Identity_per_dataset.csv')
    write_csv(mismatches_res, 'Mismathces_with_len.csv')
    

        


@click.command() 
@click.option('--nucl_aln', '-na', help="the file with processed nucleotide alignment", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--full_nucl', '-fn', help="the file with full hucleotide sequences", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--rec_tab', '-rt', help="the path to the raw RDP table", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--rec_filt', '-rf', help="the path to the filtered RDP table", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--break_tab', '-bt', help="the path to the table with predicted breakpoint coordinates", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--rdp_calc_tab', '-rc', help="the path to the table with filtered RDP events without tree congruence", 
              type=str, metavar='<STR>',default=None)
@click.option('--aln_dir', '-ad', help="the path to the directory with alignments", 
              type=str, metavar='<STR>',default=None)
@click.option('--nuc_dir', '-nd', help="the path to the domain sequences", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(nucl_aln, full_nucl, rec_tab, rec_filt, break_tab, rdp_calc_tab, aln_dir, nuc_dir):

    #File with nucleotide alignment
    align_fasta = fasta_file_to_dict(nucl_aln)
 
    #File eith full-sequence  fasta 
    full_nucl_fasta = fasta_file_to_dict(full_nucl)
    #5112 - alignment length

    #Inital raw RDP file
    RDP_init_dict = get_init_RDP_dict(read_csv_to_list(rec_tab, headless=True,delim=','))

    #Dictionary with filtered events by tree congruence
    filt_events_dict = make_dict_with_filt_events(rec_filt)

    #Breakpoints calculated from alignments
    break_dict = make_defaultdict_with_breaks_from_df(break_tab)

    #Compare sets of major and minor parents before and after alignment
    compare_maj_lists(filt_events_dict)

    #Get coordinates for RDP-obtained breakpoints in full sequences
    RDP_coords_full_adj = get_true_coordinates_from_RDP_filt(RDP_init_dict, align_fasta, full_nucl_fasta, filt_events_dict)

    #Compare breakpoint coordinates for RDP-emanated and alignment-depdendent approaches
    compare_breakpoints_alignments(RDP_coords_full_adj, break_dict, filt_events_dict)

    #Compare flanking regions from parents
    get_flanking_regions_coords(RDP_coords_full_adj, filt_events_dict, full_nucl_fasta, break_dict, aln_dir)
    
    #Dictionary with annotated RDP events
    RDP_annot_dict = make_dict_with_RDP_annot_table(rdp_calc_tab)

    #Compare the number of participants in recombination events
    compare_num_of_participants_in_events(RDP_annot_dict, filt_events_dict)

    #Read sequences of the domains
    nucl_dict = make_sequence_dict(nuc_dir)

    compare_identities_for_datasets(RDP_init_dict, RDP_annot_dict, filt_events_dict, align_fasta, nucl_dict)
    compare_identities_with_unknowns(RDP_init_dict, RDP_annot_dict, filt_events_dict, align_fasta, nucl_dict)

if __name__ == '__main__':
   main()


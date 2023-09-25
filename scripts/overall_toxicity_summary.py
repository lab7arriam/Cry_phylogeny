import os
import os.path
import click
import pandas as pd
import copy
from collections import defaultdict
from annotate_clusters_consistency import write_csv, read_csv_to_list
from get_breakpoints_from_alignments import  get_major_and_minor_parents, unpacking_list_of_lists


def get_parent_names(proc_str):
    try:
        rec_str = proc_str.split('(')[1].replace(')','').split(';')
    except:
        rec_str = proc_str.replace(')','').split(';')

    return([r.strip() for r in rec_str])



def make_strains_dictionary(strains_dat):
    
    tox_to_strain_dict = defaultdict(list)
    strain_to_tox_dict = defaultdict(list)

    strains_df = pd.read_csv(os.path.realpath(strains_dat), sep='\t')


    for index, strain_row in strains_df.iterrows():
        strain_to_tox_dict[strain_row['Strain']].append(strain_row['Cry'])
        tox_to_strain_dict[strain_row['Cry']].append(strain_row['Strain'])


    return((tox_to_strain_dict, strain_to_tox_dict))



def make_extended_stain_tox_dict(strain_to_tox_dict, tox_summary_dict, dict_type='Rank4'):

    strains_tox_extended = copy.deepcopy(tox_summary_dict['Strains'])

    for strain in strain_to_tox_dict:
        for tox in strain_to_tox_dict[strain]:
            crop_tox = crop_cry_name(tox)
            for key_index in ['Orders', 'Species']:

                if tox in tox_summary_dict['Rank4'][key_index]:
                    strains_tox_extended[key_index][strain].extend(tox_summary_dict['Rank4'][key_index][tox])

                if dict_type=='Rank3':
                    if crop_tox in tox_summary_dict['Rank3'][key_index]:
                        strains_tox_extended[key_index][strain].extend(tox_summary_dict['Rank3'][key_index][crop_tox])


    for key_index in strains_tox_extended:
        for strain in strains_tox_extended[key_index]:
            strains_tox_extended[key_index][strain]=list(set(strains_tox_extended[key_index][strain]))

    return(strains_tox_extended)



def crop_cry_name(cry_name):
    if 'NOV' in cry_name:
        id_cry = float(cry_name.split('_')[2])
        if id_cry >= 95:
            return('BT_'+ crop_cry_name(cry_name.split('_')[1]))
        else:
            return(cry_name)

    else:

        end_ind = len(cry_name)

        for symb in list(cry_name)[::-1]:
            if symb.isnumeric():
                end_ind-=1
            else:
                break

        return(cry_name[:end_ind])


def clean_tox_dict(tox_summary_dict):
    for Tox_type in tox_summary_dict:
        for Host_type in tox_summary_dict[Tox_type]:
            for key_index in tox_summary_dict[Tox_type][Host_type]:
                host_set = list(set(tox_summary_dict[Tox_type][Host_type][key_index]))
                for Host in host_set:
                    tox_summary_dict[Tox_type][Host_type][key_index] = host_set
    return(tox_summary_dict)



def update_defaultdict(update_dict,key, append_element):
    if key in update_dict:
        update_dict[key].append(append_element)
    else:
        update_dict[key]=[append_element]


def make_tox_summary_dict(tox_df):

    tox_summary_dict = defaultdict(dict)
    for dat_type in ['Rank4','Rank3','Strains']:
        tox_summary_dict[dat_type]={'Orders':defaultdict(list),'Species':defaultdict(list)}

    for index, tox_row in tox_df.iterrows():
 
        key = tox_row['Strain/Toxin']
        tox_type = tox_row['Rank'] 
        host_order = tox_row['Host_order'] 
        host_species = tox_row['Host_sp']

        if tox_type=='Strain':
            update_defaultdict(tox_summary_dict['Strains']['Orders'], key, host_order)
            update_defaultdict(tox_summary_dict['Strains']['Species'], key, host_species)

        elif tox_type=='4':
            update_defaultdict(tox_summary_dict['Rank4']['Orders'], key, host_order)
            update_defaultdict(tox_summary_dict['Rank4']['Species'], key, host_species)
            crop_key = crop_cry_name(key)

            if 'BT_' in crop_key:

                update_defaultdict(tox_summary_dict['Rank3']['Orders'], crop_key, host_order)
                update_defaultdict(tox_summary_dict['Rank3']['Species'], crop_key, host_species)

            else:

                update_defaultdict(tox_summary_dict['Rank3']['Orders'], crop_key, host_order)
                update_defaultdict(tox_summary_dict['Rank3']['Species'], crop_key, host_species)

        else:
            update_defaultdict(tox_summary_dict['Rank3']['Orders'], key, host_order)
            update_defaultdict(tox_summary_dict['Rank3']['Species'], key, host_species)

    return(tox_summary_dict)





def check_tox_for_hosts(search_term, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type='Orders', search_type = 'Rank4'):
    
    resulting_hosts = []

    if search_type == 'Rank4':
        if search_term in tox_summary_dict['Rank4'][host_type]:
            resulting_hosts.extend(tox_summary_dict['Rank4'][host_type][search_term])

    elif search_type == 'Rank3':
        crop_tox = crop_cry_name(search_term)
        if crop_tox in tox_summary_dict['Rank3'][host_type]:
            resulting_hosts.extend(tox_summary_dict['Rank3'][host_type][crop_tox])

    elif search_type == 'Strain':
        if search_term in tox_summary_dict['Strains'][host_type]:
            resulting_hosts.extend(tox_summary_dict['Strains'][host_type][search_term])


    elif search_type == 'Strain_ext':
        if search_term in extended_stain_tox_dict[host_type]:
            resulting_hosts.extend(extended_stain_tox_dict[host_type][search_term])

    return([el for el in resulting_hosts if el!='Unknown'])


def get_hosts_list(tox_list, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type='Orders', search_type = 'Rank4'):
    tox_list = [el for el in tox_list if el!='unknown']
    strains_list = []
    host_results = defaultdict(list)

    for search_term in tox_list:
        if search_type in ['Rank3','Rank4']:
            tox_spec = check_tox_for_hosts(search_term, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type, search_type)
            if len(tox_spec) > 0:
                host_results[search_term] = tox_spec

    for tox in tox_list:
        if tox in tox_to_strain_dict:
            strains_list.extend(tox_to_strain_dict[tox])
        
    strains_list = list(set(strains_list))
    
    if search_type in ['Strain','Strain_ext']:
        for search_term in strains_list:
            tox_spec = check_tox_for_hosts(search_term, tox_summary_dict, tox_to_strain_dict, extended_stain_tox_dict, host_type, search_type)
            if len(tox_spec) > 0:
                host_results[search_term] = tox_spec

    return(host_results)



def summarize_toxicity(tox_dat: str):
    tox_df = pd.read_csv(os.path.realpath(tox_dat), sep='\t')
    

    tox_summary_dict = make_tox_summary_dict(tox_df)

    tox_summary_rows = [['Strain/Toxin', 'Tox_type','Host_type', 'Host']]
    for Tox_type in tox_summary_dict:
        print(Tox_type, len(tox_summary_dict[Tox_type]['Orders']))
        for Host_type in tox_summary_dict[Tox_type]:
            for key_index in tox_summary_dict[Tox_type][Host_type]:
                host_set = list(set(tox_summary_dict[Tox_type][Host_type][key_index]))
                for Host in host_set:
                    tox_summary_rows.append([key_index, Tox_type, Host_type, Host])

    write_csv(tox_summary_rows, 'Overall_toxicity_stat.csv')



def make_asmbls_dictionary(asmbl_stat, plasm_stat, mode = 'asmbl'):
    asmbl_csv=read_csv_to_list(asmbl_stat, headless=False)
    asmbl_df=pd.DataFrame(asmbl_csv[0:], columns=['Prot_acc', 'Old_ident', 'Old_cry', 'Assembly', 'Contig', 'Level', 'New_name', 'New_ident', 'New_flag'])

    tox_to_asmbls_dict = defaultdict(list)
    asmbls_to_tox_dict = defaultdict(list)

    acc_to_plasm_dict = dict()


    plasm_csv = read_csv_to_list(plasm_stat, headless=False)
    for row in plasm_csv:
        acc_to_plasm_dict[row[1]]=row[2]


    for index, asmbl_row in asmbl_df.iterrows():
        contig_acc = asmbl_row['Contig']
        if mode == 'asmbl':
            asmbls_to_tox_dict[asmbl_row['Assembly']].append(asmbl_row['New_name'])
            tox_to_asmbls_dict[asmbl_row['New_name']].append(asmbl_row['Assembly'])

        elif mode == 'plasm':
            if contig_acc in acc_to_plasm_dict:
                if acc_to_plasm_dict[contig_acc]=='Plasmid':
                    asmbls_to_tox_dict[contig_acc].append(asmbl_row['New_name'])
                    tox_to_asmbls_dict[asmbl_row['New_name']].append(contig_acc)


        elif mode == 'chrom':
            if contig_acc in acc_to_plasm_dict:
                if acc_to_plasm_dict[contig_acc]=='Chromosome':
                    asmbls_to_tox_dict[contig_acc].append(asmbl_row['New_name'])
                    tox_to_asmbls_dict[asmbl_row['New_name']].append(contig_acc)


        elif mode == 'contig':
            asmbls_to_tox_dict[contig_acc].append(asmbl_row['New_name'])
            tox_to_asmbls_dict[asmbl_row['New_name']].append(contig_acc)


    return((reduce_dicts_to_third_rank(tox_to_asmbls_dict, asmbls_to_tox_dict)))
    #return((tox_to_asmbls_dict, asmbls_to_tox_dict))



def reduce_dicts_to_third_rank(tox_to_res_dict, res_to_tox_dict):
    tox_to_res_dict_reduced = defaultdict(list)
    res_to_tox_dict_reduced = defaultdict(list)

    for tox in tox_to_res_dict:
        crop_tox = crop_cry_name(tox)
        tox_to_res_dict_reduced[crop_tox].extend(tox_to_res_dict[tox])

    
    for res in res_to_tox_dict:
        crop_tox_list = list(set([crop_cry_name(cry) for cry in res_to_tox_dict[res]]))
        res_to_tox_dict_reduced[res] = crop_tox_list


    for tox in res_to_tox_dict_reduced:
        res_to_tox_dict_reduced[tox] = list(set(res_to_tox_dict_reduced[tox]))

    return((tox_to_res_dict_reduced, res_to_tox_dict_reduced))




def get_tox_types_for_hosts(rec_id, recs, minor_pars,major_inter_pars,major_all_pars, min_maj_strict_pars, min_maj_all_pars):
    rec_list = []
    par_strict_list = []
    par_all_list = []

    for tox in recs:
        crop_tox = crop_cry_name(tox)
        rec_list.append(crop_tox)

    for tox in min_maj_strict_pars:
        crop_tox = crop_cry_name(tox)
        par_strict_list.append(crop_tox)

    for tox in min_maj_all_pars:
        crop_tox = crop_cry_name(tox)
        par_all_list.append(crop_tox)

    return([rec_list, par_strict_list, par_all_list])


def get_tox_sets_from_rec_events(rec_tab, ref_list):

    file_reader = read_csv_to_list(os.path.realpath(rec_tab), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    rec_list_rank3 = []
    par_strict_list_rank3 = []
    par_all_list_rank3 = []


    rec_list_rank4 = []
    par_strict_list_rank4 = []
    par_all_list_rank4 = []

    for line in rec_df.iterrows():
       
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']
        
        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        min_maj_strict_pars = minor_pars+major_inter_pars
        min_maj_all_pars = minor_pars+major_all_pars

        res = get_tox_types_for_hosts(rec_id, recs, minor_pars, major_inter_pars, major_all_pars, min_maj_strict_pars, min_maj_all_pars)

        rec_list_rank3.extend(res[0])
        par_strict_list_rank3.extend(res[1])
        par_all_list_rank3.extend(res[2])
        

        rec_list_rank4.extend(recs)
        par_strict_list_rank4.extend(min_maj_strict_pars)
        par_all_list_rank4.extend(min_maj_all_pars)


    tox_sets_rank3 = [list(set([el for el in rec_list_rank3 if el!='unknown'])), list(set([el for el in par_strict_list_rank3 if el!='unknown'])), list(set([el for el in par_all_list_rank3 if el!='unknown']))]
    tox_sets_rank4 = [list(set([el for el in rec_list_rank4 if el!='unknown'])), list(set([el for el in par_strict_list_rank4 if el!='unknown'])), list(set([el for el in par_all_list_rank4 if el!='unknown']))]


    ref_res_rank3 = list(set([crop_cry_name(cry) for cry in unpacking_list_of_lists(read_csv_to_list(os.path.realpath(ref_list), headless=False))]))
    ref_no_rec_rank3 = []

    for ref_tox in ref_res_rank3:
        if ref_tox not in unpacking_list_of_lists(tox_sets_rank3):
            ref_no_rec_rank3.append(ref_tox)

    ref_res_rank4 = list(set([cry for cry in unpacking_list_of_lists(read_csv_to_list(os.path.realpath(ref_list), headless=False))]))
    ref_no_rec_rank4 = []

    for ref_tox in ref_res_rank4:
        if ref_tox not in unpacking_list_of_lists(tox_sets_rank4):
            ref_no_rec_rank4.append(ref_tox)

    
    tox_sets_rank_3 = [ref_no_rec_rank3]+tox_sets_rank3
    tox_sets_rank_4 = [ref_no_rec_rank4]+tox_sets_rank4

    tox_sets_dict_rank_3 = {'ref_no_rec': tox_sets_rank_3[0], 'rec': tox_sets_rank_3[1], 'par_strict':tox_sets_rank_3[2], 'par_all':tox_sets_rank_3[3]}
    tox_sets_dict_rank_4 = {'ref_no_rec': tox_sets_rank_4[0], 'rec': tox_sets_rank_4[1], 'par_strict':tox_sets_rank_4[2], 'par_all':tox_sets_rank_4[3]}

    return((tox_sets_dict_rank_3, tox_sets_dict_rank_4))



def extend_defaultdict_to_lists(list_dict):
    res = []
    for key in list_dict:
        res.extend(list_dict[key])
    return(res)


def compare_tox_num_hosts_per_tox_types(tox_summary_dict, rec_tab, tox_to_res_dict, res_to_tox_dict, ref_list):

    tox_sets_dict_rank_3, tox_sets_dict_rank_4 = get_tox_sets_from_rec_events(rec_tab, ref_list)

    extended_tox_to_res_dict_rank3 = make_extended_stain_tox_dict(res_to_tox_dict, tox_summary_dict, dict_type='Rank3')
    extended_tox_to_res_dict_rank4 = make_extended_stain_tox_dict(res_to_tox_dict, tox_summary_dict, dict_type='Rank4')

    
    num_hosts_for_toxin_groups = [['Tox_type', 'Tox', 'Host_type', 'Group', 'Num_hosts']]

   
   #per toxins:
    for host_type in ['Species', 'Orders']:
        for mode in ['Rank3', 'Rank4', 'Strain']: # 'Strains_add'
            for tox_set in tox_sets_dict_rank_3:   
                for tox in tox_sets_dict_rank_3[tox_set]:
                    if mode=='Strain':
                        hosts_rank =  get_hosts_list([tox],tox_summary_dict, tox_to_res_dict, extended_tox_to_res_dict_rank4, host_type, 'Rank3')
                        hosts_strains =  get_hosts_list([tox],tox_summary_dict, tox_to_res_dict, extended_tox_to_res_dict_rank4, host_type, mode)

                        if len(hosts_rank)>0:
                            hosts = set(extend_defaultdict_to_lists(hosts_rank))
                        else:
                            hosts = set(extend_defaultdict_to_lists(hosts_strains))

                    else:
                        hosts =  set(extend_defaultdict_to_lists(get_hosts_list([tox],tox_summary_dict, tox_to_res_dict, extended_tox_to_res_dict_rank4, host_type, mode)))

                    num_hosts_for_toxin_groups.append([tox_set, tox, host_type, mode, len(hosts)])

    write_csv(num_hosts_for_toxin_groups, 'num_hosts_per_toxin_group.csv')


def compare_num_hosts_for_datasets_with_tox_types(tox_summary_dict, rec_tab, tox_to_res_dict, res_to_tox_dict, ref_list, mode = 'Strains'):

    tox_sets_dict_rank_3, tox_sets_dict_rank_4 = get_tox_sets_from_rec_events(rec_tab, ref_list)

    num_hosts_for_dataset_groups = []

    for host_type in ['Species', 'Orders']: #, 'Orders'
        for tox_set in tox_sets_dict_rank_3:
            for group, tox_sets_dict in zip(['Rank4', 'Rank3'],[tox_sets_dict_rank_4, tox_sets_dict_rank_3]):

                res_set = list(set(unpacking_list_of_lists([tox_to_res_dict[tox] for tox in tox_sets_dict[tox_set]])))
            
                for res in res_set:
                    host_str_all = tox_summary_dict['Strains'][host_type][res].copy()
                    host_str_trim =  tox_summary_dict['Strains'][host_type][res].copy()

                    for tox in res_to_tox_dict[res]:
                        if len(tox_summary_dict['Rank3']['Species'][tox]) >= 15:
                            host_str_all.extend(tox_summary_dict['Rank3'][host_type][tox])

                        else:
                            #print(tox, 'trim')
                            host_str_all.extend(tox_summary_dict['Rank3'][host_type][tox])
                            host_str_trim.extend(tox_summary_dict['Rank3'][host_type][tox])


                    num_hosts_for_dataset_groups.append([tox_set, res, host_type,group, mode, len(set(host_str_trim)),'Trimmed' ])
                    num_hosts_for_dataset_groups.append([tox_set, res, host_type,group, mode, len(set(host_str_all)),'Untrimmed' ])
                        

    return(num_hosts_for_dataset_groups)


@click.command()           
@click.option('--tox_dat', '-t', help="The file with toxicity spectrum", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--rec_tab', '-r', help="the path with recombination results", 
              type=click.Path(exists=True), metavar='<PATH>')    
@click.option('--strains_stat', '-s', help="the path to the table denoting strains attributions", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--ref_list', '-f', help="the path to the list of reference proteins", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--asmbl_stat', '-a', help="the path to cry proteins disctribution on assemblies", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--plasm_stat', '-p', help="the path to attributions of nucleotide acessions in the assemblies", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(tox_dat, rec_tab, strains_stat, ref_list, asmbl_stat, plasm_stat): 
    summarize_toxicity(tox_dat)
    tox_to_strain_dict, strain_to_tox_dict =  reduce_dicts_to_third_rank(*make_strains_dictionary(strains_stat))
    tox_to_asmbls_dict, asmbls_to_tox_dict =  make_asmbls_dictionary(asmbl_stat, plasm_stat, mode = 'asmbl')

    tox_df = pd.read_csv(os.path.realpath(tox_dat), sep='\t')
    tox_summary_dict = clean_tox_dict(make_tox_summary_dict(tox_df))

    compare_tox_num_hosts_per_tox_types(tox_summary_dict, rec_tab, tox_to_strain_dict, strain_to_tox_dict, ref_list)

    num_hosts_for_dataset_groups = [['Tox_type', 'Accession', 'Host_type', 'Group', 'Mode', 'Num_hosts', 'Trim_flag']]
    num_hosts_for_dataset_groups.extend(compare_num_hosts_for_datasets_with_tox_types(tox_summary_dict, rec_tab, tox_to_asmbls_dict, asmbls_to_tox_dict, ref_list, 'Assemblies'))
    num_hosts_for_dataset_groups.extend(compare_num_hosts_for_datasets_with_tox_types(tox_summary_dict, rec_tab, tox_to_strain_dict, strain_to_tox_dict, ref_list, 'Strains'))

    write_csv(num_hosts_for_dataset_groups, 'num_hosts_per_dataset_in_toxin_group.csv')

if __name__ == '__main__':
   main()

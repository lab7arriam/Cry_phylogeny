import os
import click
import copy
import pandas as pd
from collections import defaultdict
from create_evol_data_on_fixed_recombinants import get_parent_names
from annotate_clusters_consistency import write_csv, read_csv_to_list, create_dict_from_list


def get_recombination_tox_list(rec_file: str) -> list:
    rec_list = []
    rec_res = read_csv_to_list(os.path.realpath(rec_file), headless=False)
    rec_df = pd.DataFrame(rec_res[1:], columns=rec_res[0])

    for line in rec_df.iterrows(): 
        for index in ['parent_1','parent_2','parent_3']:
            rec_list.extend([tox for tox in get_parent_names(line[1][index]) if tox!='unknown'])
        rec_list.extend([tox for tox in get_parent_names(line[1]['Rec'])])

    return(list(set(rec_list)))

def extend_rec_with_clusters(rec_list: str, clust_file: str) -> list:
    rec_extend = rec_list.copy()
    clust_res = read_csv_to_list(os.path.realpath(clust_file), headless=True)
    for clust_line in clust_res:
        if clust_line[0] in rec_list:
            rec_extend.extend(clust_line[7].split(';'))

    rec_extend=list(set(rec_extend))
    return(rec_extend)
    

def parse_MGE_bed(MGE_bed: str, extr_type = 'global') -> list:
    bed_res = read_csv_to_list(os.path.realpath(MGE_bed), headless=False)
    MGE_list = []
    MGE_dict = defaultdict(list)

    for line in bed_res:
        MGE_list.append(line[3].replace('ID=',''))
        if extr_type == 'cry':
            MGE_dict[line[3].replace('ID=','')].append(line[7])

    if extr_type == 'global':
        return(MGE_list)
    else:
        return(MGE_dict)


def get_MGEs_from_dir(mge_subdir: str, MGE_type: str, extr_type = 'global') -> defaultdict(dict):

    MGE_asmbl_dict = defaultdict(dict)

    for mge_bed in os.listdir(mge_subdir):
        mge_bed_path = os.path.join(os.path.realpath(mge_subdir), mge_bed)
        asmbl_name = '_'.join([mge_bed.split('_')[0], mge_bed.split('_')[1]])

        if extr_type == 'global':
            MGE_asmbl_dict[asmbl_name] = dict()
            MGE_asmbl_dict[asmbl_name][MGE_type] = parse_MGE_bed(mge_bed_path)

        else:
            MGE_asmbl_dict[asmbl_name] = dict()
            MGE_asmbl_dict[asmbl_name][MGE_type] = defaultdict(list)
            MGE_crys = parse_MGE_bed(mge_bed_path, extr_type = extr_type)

            for mge in MGE_crys:
                MGE_asmbl_dict[asmbl_name][MGE_type][mge] = MGE_crys[mge]

    return(MGE_asmbl_dict)


def filter_asmbl_from_mge_dict(asmbl_dict: dict, MGE_dict: defaultdict(dict), MGE_type: str) -> defaultdict(dict):
    MGE_dict_redone = defaultdict(dict)
    for asmbl in asmbl_dict:
        MGE_dict_redone[asmbl] = dict()
        MGE_dict_redone[asmbl][MGE_type] = defaultdict(list)
        if asmbl in MGE_dict:
            MGE_dict_redone[asmbl][MGE_type] = MGE_dict[asmbl][MGE_type]

    return(MGE_dict_redone)


def summarize_MGEs(mge_dir: str, asmbl_dict: dict, extr_type = 'global') -> list:

    if extr_type == 'global':
        parse_subdir = 'source'
    else:
        parse_subdir = 'cry_inter'

    MGE_types_dict = {'GIs':'GIs', 'IS':'IS', 'phages':'phrophages'}
    MGE_res_list = []

    for mge_type_key in MGE_types_dict:
        mge_subdir = os.path.join(os.path.realpath(mge_dir), mge_type_key+'/'+parse_subdir)
        mge_res_dict = get_MGEs_from_dir(mge_subdir, MGE_types_dict[mge_type_key], extr_type = extr_type)

        if  extr_type == 'global': 
            mge_res_filt = filter_asmbl_from_mge_dict(asmbl_dict, mge_res_dict, MGE_types_dict[mge_type_key])
            MGE_res_list.append(mge_res_filt)
        else:
            if len(mge_res_dict) > 0:
                MGE_res_list.append(mge_res_dict)

    return(MGE_res_list)

def make_table_with_MGEs_nums(MGE_res_list: list, asmbl_dict: dict):

    mge_num_res = [['Accession', 'MGE_type', 'MGE_num', 'Dataset']]
    mge_num_res_non_redundant = [['Accession', 'MGE_type', 'MGE_num', 'Dataset']]

    for asmbl in asmbl_dict:
        asmbl_type = asmbl_dict[asmbl][1]
        for mge_dict in MGE_res_list:
            mge_type = ''.join(list(mge_dict[asmbl].keys()))
            mge_num = len(mge_dict[asmbl][mge_type])
            mge_num_res.append([asmbl, mge_type, mge_num, asmbl_type])
            mge_num_res_non_redundant.append([asmbl, mge_type, mge_num, asmbl_type])

            if asmbl_type == 'rec':
                mge_num_res.append([asmbl, mge_type, mge_num, 'cry'])
                mge_num_res.append([asmbl, mge_type, mge_num, 'all'])

            if asmbl_type == 'cry':
                mge_num_res.append([asmbl, mge_type, mge_num, 'all'])

    write_csv(mge_num_res,'MGES_num_summary.csv')
    write_csv(mge_num_res_non_redundant,'MGES_num_summary_non_redundant.csv')

def get_cry_list_asmbl(cry_dir: str) -> defaultdict(list):
    cry_asmbl_dict = defaultdict(list)

    for cry_bed in os.listdir(cry_dir):
        cry_bed_path = os.path.join(os.path.realpath(cry_dir), cry_bed)
        asmbl_name = '_'.join([cry_bed.split('_')[0], cry_bed.split('_')[1]])
        cry_asmbl_dict[asmbl_name] = parse_MGE_bed(cry_bed_path)

    return(cry_asmbl_dict)

def get_crys_from_mge(mge_dict_cry: defaultdict(dict)):
    cry_list = []
    for mge in mge_dict_cry:
        for cry in mge_dict_cry[mge]:
            cry_list.append(cry)
    return(cry_list)

def get_num_recs_from_mge(rec_list: list, mge_dict_cry: defaultdict(dict)):
    rec_crys = []
    for mge in mge_dict_cry:
        for cry in mge_dict_cry[mge]:
           if cry in rec_list:
               rec_crys.append(cry)
    return(rec_crys)
     
def summarize_crys_in_asmbls(MGE_cry_list: list, MGE_glob_list: list, asmbl_dict: dict, rec_list: list, cry_asmbl_dict: list):

    MGE_crys_stat = [['Accession', 'Dataset', 'MGE_type', 'MGE_num_glob', 'MGE_num_cry','Cry_num', 'Rec_num', 'Cry_mge', 'Rec_mge']]
    total_stat = [['Accession', 'Cry_num', 'Rec_num', 'Cry_GIs', 'Rec_GIs', 'Cry_IS', 'Rec_IS']]

    total_crys = 0
    total_recs = 0


    GIs_crys = 0
    GIs_recs = 0

    IS_crys = 0
    IS_recs = 0

    for asmbl in cry_asmbl_dict:
        asmbl_type = asmbl_dict[asmbl][1]
        cry_list = cry_asmbl_dict[asmbl]
        cry_num = len(cry_list)
        rec_num = len([cry for cry in cry_list if cry in rec_list])
        total_crys += cry_num
        total_recs += rec_num

        for ind in [0,1]:
            mge_dict_glob = MGE_glob_list[ind]
            mge_type_glob = ''.join(list(mge_dict_glob[asmbl].keys()))
            mge_num_glob = len(mge_dict_glob[asmbl][mge_type_glob])

            mge_dict_cry = MGE_cry_list[ind]

            if asmbl in mge_dict_cry:
                mge_num_cry = len(mge_dict_cry[asmbl][mge_type_glob])
                crys_in_mges = len(get_crys_from_mge(mge_dict_cry[asmbl][mge_type_glob]))
                recs_in_mges = len(get_num_recs_from_mge(rec_list, mge_dict_cry[asmbl][mge_type_glob]))

                if mge_type_glob=='GIs':
                    GIs_crys += crys_in_mges
                    GIs_recs += recs_in_mges

                elif mge_type_glob=='IS':
                    IS_crys += crys_in_mges
                    GIs_recs += recs_in_mges

            else:
                mge_num_cry = 0
                crys_in_mges = 0
                recs_in_mges = 0

            MGE_crys_stat.append([asmbl, asmbl_type, mge_type_glob, mge_num_glob, mge_num_cry, cry_num, rec_num, crys_in_mges, recs_in_mges])

    write_csv(MGE_crys_stat,'MGES_cry_summary.csv')
    print(total_crys, GIs_crys, IS_crys, total_recs, GIs_recs, IS_recs)


@click.command()
@click.option('--rec_file', '-r', help="recombination results", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--asmbl_tab', '-a', help="list of assemblies", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--mge_dir', '-m', help="directory with MGEs predictions", 
              type=click.Path(exists=True), metavar='<PATH>') 
@click.option('--cry_dir', '-c', help="directory with cry beds", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--clust_res', '-cl', help="clusterization results", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(rec_file, asmbl_tab, mge_dir, cry_dir, clust_res):
    rec_list = get_recombination_tox_list(rec_file)
    #rec_list =  extend_rec_with_clusters(rec_list, clust_res)

    asmbl_dict = create_dict_from_list(read_csv_to_list(os.path.realpath(asmbl_tab), headless=False), key_ind = 0 )

    mge_nums = summarize_MGEs(mge_dir, asmbl_dict)
    mge_crys = summarize_MGEs(mge_dir, asmbl_dict, extr_type = 'cry')
   
    make_table_with_MGEs_nums(mge_nums, asmbl_dict)

    cry_asmbl_dict = get_cry_list_asmbl(cry_dir)
    summarize_crys_in_asmbls(mge_crys, mge_nums, asmbl_dict, rec_list, cry_asmbl_dict)
    
if __name__ == '__main__':
   main()

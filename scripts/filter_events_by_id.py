import argparse
import csv
from collections import defaultdict
from annotate_clusters_consistency import write_csv, read_csv_to_list, create_dict_from_list
import pandas as pd
import os.path

def choose_ids(row: pd.Series, sel_type)-> tuple:

    if row['Type'].startswith('domain3'):
        if sel_type=='recs':
            #return(row['rec_3'],row['r_min_3'],row['r_maj_3'])
            ###results are identical
            return(row['rec_a'],row['r_min_a'],row['r_maj_a'])

        elif sel_type=='parents':
            return(row['r_min_3'],row['r_min_1'],row['r_min_2'])

    if row['Type'].startswith('domain2'):
        if sel_type=='recs':
            #return(row['rec_2'],row['r_min_2'],row['r_maj_2'])
            return(row['rec_a'],row['r_min_a'],row['r_maj_a'])

        elif sel_type=='parents':
            return(row['r_min_2'],row['r_min_1'],row['r_min_3'])

    if row['Type'].startswith('domain1'):
        if sel_type=='recs':
            #return(row['rec_1'],row['r_min_1'],row['r_maj_1'])
            return(row['rec_a'],row['r_min_a'],row['r_maj_a'])

        elif sel_type=='parents':
            return(row['r_min_1'],row['r_min_2'],row['r_min_2'])


def compare_domains_results(id1,id2,id3,mode) -> bool:
   if mode=='forward':
       return(float(id1) >= sum([float(id2), float(id3)])/2)

   elif mode=='reverse':
       return(float(id1) < sum([float(id2), float(id3)])/2)


def filter_by_recomb_id(row: pd.Series, mode, sel_type) -> list:

    result = row.copy()

    if compare_domains_results(*choose_ids(result,sel_type=sel_type),mode):
        return(list(result))


def iterate_over_df_rows(recomb_df: pd.DataFrame, mode, sel_type) -> list:

    ret_list=[]
    for index, row in recomb_df.iterrows():
        filt = filter_by_recomb_id(row, mode=mode, sel_type=sel_type)
        if filt:
             ret_list.append(filt)
    return(ret_list)
    

def parse_recombination_events(recomb_list, out_dir):

    col_names=recomb_list[0]
    recomb_df=pd.DataFrame(recomb_list[1:], columns = recomb_list[0])

    print('recs filtered',iterate_over_df_rows(recomb_df,mode='reverse', sel_type='recs') )
    print('pars filtered',iterate_over_df_rows(recomb_df,mode='reverse', sel_type='parents') )

    #write_csv(iterate_over_df_rows(recomb_df,mode='reverse', sel_type='recs'),os.path.join(os.path.realpath(out_dir),'rec_id_not_passed.csv'),*col_names)
    write_csv(iterate_over_df_rows(recomb_df,mode='reverse', sel_type='parents'),os.path.join(os.path.realpath(out_dir),'parents_id_not_passed.csv'),*col_names)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyzes RDP-detected recombination events and filers them based on sequence identity')
    parser.add_argument('-r', '--rec_f', dest='rec_file', help='path to the recombination table file',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_dir', help='the path to the output directory',
                        type=str)
    args = parser.parse_args()
    recomb_df=parse_recombination_events(read_csv_to_list(args.rec_file, headless=False,delim='\t'), args.out_dir)


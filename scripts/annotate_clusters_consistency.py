#/usr/bin/python3.7

import argparse
import csv
from collections import defaultdict
import os.path
from Bio import SeqIO, Entrez
from Bio import pairwise2 as pw2
from itertools import combinations, chain
import re
import time


#rename patent names
rename_dict={"NOV_Cry2Ah2_92.7_AAQ52362.1":"NOV_Cry2Ah2_92.7_AX098631.1", "NOV_Cry40Aa1_67.9_ACC34441.1":"NOV_Cry40Aa1_67.9_CQ868314.1", "NOV_Cry8Bb1_98.4_ABJ38812.1":"NOV_Cry8Bb1_98.4_CS131028.1", "NOV_Cry2Ak1_76.8_ABU37568.1":"NOV_Cry2Ak1_76.8_AX513524.1", "NOV_Cry1Ga1_77.9_AAQ52381.1":"NOV_Cry1Ga1_77.9_AX098669.1","NOV_Cry1Ab18_92.4_AAE49246.1":"NOV_Cry1Ab18_92.4_AX383797.1", "BT_Cry1Jc1":"NOV_Cry1Jc1_99.1_AAS93799.1", "NOV_Cry24Aa1_98.3_WP_086403584.1":"NOV_Cry24Aa1_98.3_WP_086403584.1" }

def find_sequence_identity(seq_tuple: tuple) -> float:
    """
      Performes global pairwise alignmet of two sequences
      Returns blast-like identity between sequences
    """

    global_align = pw2.align.globalxx(seq_tuple[0], seq_tuple[1])
    seq_length = min(len(seq_tuple[0]), len(seq_tuple[1]))

    matches = global_align[0][2]
    percent_match_loc = (matches / seq_length) * 100

    al_length = global_align[0][4]
    percent_match_glob = (matches / al_length) * 100

    return([percent_match_loc, percent_match_glob])

def write_csv(row_list: list,out_name: str,*header_strings : str):
    """
       A universal function for writing lists to csv-files
       If input is list of lists uses writerows function else iteratively writes 
       If strings for header are specified, writes header to the output, saves headless table otherwise
    """
    with open(out_name,'w',newline='',encoding='utf-8') as result_file:
        wr = csv.writer(result_file, delimiter='\t')
        if header_strings:
            wr.writerow([name for name in header_strings])
        if type(row_list[0]) is list:
            wr.writerows(row_list)
        else:
            for row in row_list:
                wr.writerow([row])

def create_dict_from_list(parse_list: list, key_ind: int, *val_inds: int, multi=False) -> defaultdict(list):
    """
        Creates list dict from parsing list
        Key refers to the asserted index (key_ind), value - by default takes all the indexes (including the key index)
    """
    parse_dict=defaultdict(list)
    for string in parse_list:
        if not val_inds:
            if not multi:
                parse_dict[string[key_ind]]=list(string)
            else:
                if string[key_ind] not in parse_dict:
                    parse_dict[string[key_ind]]=[list(string)]
                else:
                    parse_dict[string[key_ind]].append(list(string))
        else:
            if not multi:
                parse_dict[string[key_ind]]=[string[i] for i in range(len(string)) if i in val_inds]
            else:
                if string[key_ind] not in parse_dict:
                    parse_dict[string[key_ind]]=[[string[i] for i in range(len(string)) if i in val_inds]]
                else:
                    parse_dict[string[key_ind]].append([string[i] for i in range(len(string)) if i in val_inds])
    return(parse_dict)


def read_csv_to_list(in_file: str, headless=True, delim='\t') -> list:
    """
        Reads csv file and returns list without header by default 
        If headless argument is false, parses the whole file
    """
    ret_list=list()
    with open(in_file,'r',encoding='utf-8') as csv_file:
        my_reader = csv.reader(csv_file, delimiter=delim) 
        if headless:
            next(my_reader)
        for row in my_reader:
            ret_list.append(list(row))
    return(ret_list)

def clean_bt_name(toxin: str) -> str:
    """
        Assings name to the raw fasta headings with 3-d Cry toxins
        Adds 'NOV' to the name in case it is a novel toxin and 'BT' otherwise

    """
    new_id = toxin.replace('(+)', '').replace('(-)', '').replace(')_', '_').replace('(', '_')
    if '_' in new_id:
        new_id='NOV_'+new_id
    else:
        new_id='BT_'+new_id
    return(new_id)

def parse_cd_hit_clusters(in_file: str, delim='\t') -> list:
    """
        Reads .cluster CD-HIT file
        For each cluster returns the number of sequences and the length distribution
    """

    ret_row=[]
    clust_dict=defaultdict(dict)

    with open(in_file,'r') as csv_file:
        my_reader3 = csv.reader(csv_file, delimiter=delim)
        for row in my_reader3:
            if row[0][0]==">":
                clust=row[0].split(' ')[1]
                clust_dict[clust]['length']=list()
                clust_dict[clust]['name']=''
                clust_dict[clust]['names']=list()
            else:
                if '*' in row[1].split():
                    cry_name=row[1].split()[1].replace('>', '').replace('...', '')
                    clust_dict[clust]['name']=cry_name
                    clust_dict[clust]['names'].append(cry_name)
                else:
                    cry_name=row[1].split()[1].replace('>', '').replace('...', '')
                    clust_dict[clust]['names'].append(cry_name)
                    clust_dict[clust]['length'].append(row[1].split()[0].replace('aa,',''))

    for key in clust_dict:

        try:
            length_range = max([int(x) for x in clust_dict[key]['length']]) - min([int(x) for x in clust_dict[key]['length']])

        except(ValueError):
            length_range = 0

        ret_row.append([clean_bt_name(clust_dict[key]['name']),clust_dict[key]['name'],';'.join(list(set(clust_dict[key]['length']))),
                       length_range,len(clust_dict[key]['names']),len(set([clean_bt_name(x) for x in list(set(clust_dict[key]['names']))])),';'.join(list(set(clust_dict[key]['names']))),';'.join([clean_bt_name(x) for x in list(set(clust_dict[key]['names']))])])

    return(ret_row)

def filter_rec_clusters(clusters: list, recs: list) -> list:
    """
        Parses file with CD-HT clusters and returns toxins pertaining to the recombinant events
    """

    rec_nest=[rec[1].split(';') for rec in recs if len(rec[1].split(';'))<=3]
    rec_list=list(set([obj for recs in rec_nest for obj in recs]))
    ret_list=[]

    for row in clusters:
        if row[0] in rec_list:
            ret_list.append(row)
    return(ret_list)

def replace_patent_names(parsed_clusters: list) -> defaultdict(list):
    """
        Reads the dictionary of cluster names and replaces new names related to patents
    """

    refs_dict=defaultdict(list)

    for row in parsed_clusters:
        if row[0]!='BT_Cry41Ca1':
            for key, value in zip(rename_dict.keys(),rename_dict.values()):
                row[0]=row[0].replace(key,value)
                row[7]=row[7].replace(key,value)
            refs_dict[row[0]]=row

    return(refs_dict)


def get_intra_cluster_identity (refs_dict: defaultdict(list), recs_dict: dict, mode='nucl') -> list:
    """
        Gets a dictionary with sequences and a dictionary with clusters
        Returns a list of clusters' annotation with the number of unique sequences
    """
    ret_rows=list()

    if mode=='nucl':
        parse_ind=7
    else:
        parse_ind=9

    for clust in refs_dict:
        print(clust)
        #calculate non-duplicate sequences compared to reference
        filt_inter=set()
        if len(refs_dict[clust][parse_ind].split(';'))>1:
            ref_seq=recs_dict[refs_dict[clust][0]]
            non_dups=list()
            for rec in refs_dict[clust][parse_ind].split(';'):
                if rec!=refs_dict[clust][0]:
                    try:
                        if mode=='nucl':
                            id_thr=find_sequence_identity((ref_seq,recs_dict[rec]))[0]
                        else:
                            id_thr=find_sequence_identity((ref_seq.translate(),recs_dict[rec].translate()))[0]
                        if id_thr<100.0:
                            non_dups.append(rec) 
                    except:
                        pass

            filt_inter=set()
            
            if len(non_dups)>1:
                #calculate non-duplicate sequences in pre-filtered set
                comb_list=list(combinations(non_dups,2))
                non_passed=set()
                filt_inter.add(comb_list[0][0])
                for pair in list(combinations(non_dups,2)):
                    if mode=='nucl':
                        th2=find_sequence_identity((recs_dict[pair[0]],recs_dict[pair[1]]))[0]
                    else:
                        th2=find_sequence_identity((recs_dict[pair[0]].translate(),recs_dict[pair[1]].translate()))[0]
                    if pair[0] not in non_passed:
                        if th2==100.0:
                            non_passed.add(pair[1])
                        if pair[0] not in non_passed and th2!=100.0:
                            filt_inter.add(pair[0])
        

        ret_rows.append(refs_dict[clust] + [len(filt_inter)+1, ';'.join([clust]+list(filt_inter))])

    return(ret_rows)

def get_seq_dict(ref_nucl: str, clust_nucl: str) -> dict:
    """
        Reads the paths to the reference and clusters sequences and returns a dictionary with fasta-records
    """
    recs_dict=dict()
    for record in SeqIO.parse(os.path.relpath(ref_nucl),"fasta"):
        recs_dict[record.id]=record.seq.upper()

    for record in SeqIO.parse(os.path.relpath(clust_nucl),"fasta"):
        recs_dict[record.id]=record.seq.upper()

    return(recs_dict)

def calculate_nucleotide_non_duplicates(parsed_clusters: list, ref_nucl: str, clust_nucl: str) -> list:
    """
        Reads the list with parsed clusters and nucleotide sequences for references and for within-cluster sequences
        Calculates the number of unique sequences in the cluster and returnes the annotated list
    """
 
    refs_dict=replace_patent_names(parsed_clusters)
    
    recs_dict=get_seq_dict(ref_nucl, clust_nucl)
    ret_rows= get_intra_cluster_identity (refs_dict, recs_dict)

    return(ret_rows)


def get_ipg_pat_list(annot_list: list, bt_dict: dict, tox_list: list, out_dir: str) -> list:
    """
       Reads annotated cluster table and Bt nomenclature accessions and searches patent dbsource with IPG
       Returns an annotated table and writes IPG fetching results 
    """  
    raw_dict=dict()
    filt_dict=dict()
    ret_list=annot_list
    acc_list=list()
    full_IPG_list=list()
    
    for name in tox_list:
        
        raw_dict[name]=list()
        if 'BT' in name:
            acc=bt_dict[name]
        else:
            acc=re.split("\.\d_", name)[1]

        try:
            handle = Entrez.efetch(db="protein",rettype='ipg',retmode='text',id=acc)

        except:
            print('No HTTP acess')
            continue
        try:
            handle_list=[el.split('\t') for el in handle.read().split('\n')][1:-1]

            for row in handle_list:
               raw_dict[name].append(row[1])
               full_IPG_list.append([name]+ row)
            print(name, acc, len(handle_list))
            time.sleep(1)
        except:
            try:
                record = Entrez.read(Entrez.elink(dbfrom="nucleotide", id=acc, db='protein'))
                prim_id = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]][0]
                time.sleep(1)
                handle = Entrez.efetch(db="protein",rettype='ipg',retmode='text',id=prim_id) 
                handle_list=[el.split('\t') for el in handle.read().split('\n')][1:-1]
                
                for row in handle_list:
                   raw_dict[name].append(row[1])
                   full_IPG_list.append([name]+ row)
                print(name, acc, len(handle_list))
            except:
                print('Bad accession!', name, acc)
        time.sleep(1)

    for name in raw_dict:
        if len(raw_dict[name])>0:
            if sum(['PAT' in db_source for db_source in raw_dict[name]])<len(raw_dict[name]):
                filt_dict[name]='not_PAT'
            else:
                filt_dict[name]='PAT'
        else:
            filt_dict[name]='not_PAT'

    for cl_ind in range(len(ret_list)):
        ref_flag=filt_dict[ret_list[cl_ind][0]]
        pat_clust=[]
        for tox in ret_list[cl_ind][9].split(';'):
            pat_clust.append(filt_dict[tox])
  
        ret_list[cl_ind].extend([ref_flag, len([el for el in pat_clust if el=='PAT'])])

    write_csv(full_IPG_list,os.path.join(os.path.realpath(out_dir),'IPG_table_for_clusters.tsv'))

    return(ret_list)

def expand_annotated_table(annot_clust: list, ref_nucl: str, clust_nucl: str, bt_nom_table: str, mail_adr: str, out_dir: str) -> list:
    """
        Reads the list with annotated clusters and expands it with amino acid identity and the number of patent sequences
    """

    refs_dict=replace_patent_names(annot_clust)
    recs_dict=get_seq_dict(ref_nucl, clust_nucl)
    annot_list=get_intra_cluster_identity (refs_dict, recs_dict,mode = 'prot')
    Entrez.email = mail_adr

    
    bt_parse_list=read_csv_to_list(bt_nom_table, headless=False,delim='\t')
    bt_dict=dict()
    tox_list = list(chain.from_iterable([cl_list[9].split(';') for cl_list in annot_list]))
    bt_list=[tox for tox in tox_list if 'BT' in tox]
     
    for row in bt_parse_list:
        if row[1]!='\xa0':
            bt_name='BT_' + row[0]
            if bt_name in bt_list:
                bt_dict[bt_name]=row[1].strip().replace(' ','').replace('\r','')

    return(get_ipg_pat_list(annot_list, bt_dict, tox_list, out_dir)) 
                           

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses CD-HT clusterins results and returns a table with sequences\' lenghts')
    parser.add_argument('-cf', '--cl_f', dest='clust_file', help='path to the CD-HIT cluster file',
                        type=str)
    parser.add_argument('-o', '--out', dest='out_dir', help='the path to the output directory',
                        type=str)
    parser.add_argument('-r', '--rec', dest='rec_file', help='the path to the file with recombination events',
                        type=str)
    parser.add_argument('-rn', '--ref', dest='ref_nucl', help='the path to the reference nucleotide sequences',
                        type=str)
    parser.add_argument('-cn', '--clust', dest='clust_nucl', help='the path to the cluster nucleotide sequences',
                        type=str)
    parser.add_argument('-ct', '--cl_tab', dest='clust_tab', help='the path to the anotated table',
                        type=str)
    parser.add_argument('-m', '--mail_a', dest='mail_adr', help='e-mail adress',
                        type=str)
    parser.add_argument('-bt', '--bt_t', dest='bt_table', help='path to BT nomenclature table',
                        type=str)
    parser.add_argument('-ap', '--ap_k', dest='api_key', help='NCBI API key',
                        type=str)
    args = parser.parse_args()

    parsed_clusters=parse_cd_hit_clusters(args.clust_file)
    recs=read_csv_to_list(args.rec_file, headless=True,delim='\t')
    Entrez.api_key = args.api_key

    if not args.clust_tab:
        #write_csv(parsed_clusters,os.path.join(os.path.realpath(args.out_dir),'cd_hit_clusters_length.tsv'),'new_ref_name','init_ref_name','length_list',
        #                                                                                                               'length_range','num_clusters','dedup_clusters','seq_names_old','seq_names_new')
        #write_csv(filter_rec_clusters(parsed_clusters, recs),os.path.join(os.path.realpath(args.out_dir),'recomb_cd_hit_clusters.tsv'),'new_ref_name','init_ref_name','length_list',
        #                                                                                                               'length_range','num_clusters','dedup_clusters','seq_names_old','seq_names_new')
        calculated_clust=calculate_nucleotide_non_duplicates(parsed_clusters,args.ref_nucl,args.clust_nucl)
        write_csv(calculated_clust,os.path.join(os.path.realpath(args.out_dir),'no_pat_cd_hit_clusters.tsv'),'new_ref_name','init_ref_name','length_list',
                  'length_range','num_clusters','dedup_clusters','seq_names_old','seq_names_new','unique_nucl_full','unique_nucl_names')
        write_csv(filter_rec_clusters(calculated_clust, recs),os.path.join(os.path.realpath(args.out_dir),'no_pat_recs_cd_hit_clusters.tsv'),'new_ref_name','init_ref_name','length_list',
                  'length_range','num_clusters','dedup_clusters','seq_names_old','seq_names_new','unique_nucl_full','unique_nucl_names')
    else:
        annot_clust=read_csv_to_list(args.clust_tab, headless=True,delim='\t')
        write_csv(expand_annotated_table(annot_clust, args.ref_nucl,args.clust_nucl, args.bt_table, args.mail_adr,args.out_dir),os.path.join(os.path.realpath(args.out_dir),'clusters_PAT_annotation.tsv'),'new_ref_name','init_ref_name','length_list',
                  'length_range','num_clusters','dedup_clusters','seq_names_old','seq_names_new','unique_nucl_full','unique_nucl_names','unique_prot_full','unique_prot_names', 'ref_PAT', 'PAT_num')


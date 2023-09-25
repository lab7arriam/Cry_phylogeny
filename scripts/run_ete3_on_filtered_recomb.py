from collections import defaultdict
import argparse
import pandas as pd
from typing import List
import os.path
import subprocess
from ete3 import Tree
import random

def return_event_index(row: pd.Series):
    ev_num = str(row['ID'])
    ret_list=list()
    if not 'partial' in row['Type']:
        if int(row['Min_evol'])>5 and int(row['Min_evol'])<70:
            ret_list.append('minor_event_'+ev_num)
        if int(row['Maj1_evol'])>5 and int(row['Maj1_evol'])<70:
            ret_list.append('major1_event_'+ev_num)
        if int(row['Maj2_evol'])>5 and int(row['Maj2_evol'])<70:
            ret_list.append('major2_event_'+ev_num)
    return ret_list

def mark_nodes(recs: str):
    parsed=recs.split(';')
    if len(parsed)==1:
        return(parsed[0])
    elif len(parsed)>1:
        return(parsed[0]+',,,'+parsed[-1])
    else:
        return(None)


def get_process_row(evol_dir):
    event_set=set()
    def process_row(row: pd.Series):
        result = row.copy()
        events=return_event_index(result)
        if len(events)>0:
            for event in events:
                if event in event_set:
                    break
                event_set.add(event)

                print('Analyzing event {} .........'.format(event))
                event_subdir = os.path.join(os.path.realpath(evol_dir),event)
                names_tup = tuple([os.path.join(event_subdir,x) for x in  ['site_model', 'branch_model_rec', 'branch_site_model_rec']])

                if 'minor' in event:
                    rec_mark=mark_nodes(result['Rec_min'])
                    if 'unknown' not in result['Min_par']:
                        names_tup=tuple([os.path.join(event_subdir,x) for x in  ['site_model', 'branch_model_rec', 'branch_site_model_rec','branch_model_par', 'branch_site_model_par']])
                        par_mark=mark_nodes(result['Min_par'])
                        dir_call = subprocess.call("mkdir {0} {1} {2} {3} {4}".format(*names_tup),shell=True)
                        par_bs_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta','branch_site_model_par', par_mark),shell=True)
                        par_b_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'branch_model_par', par_mark),shell=True)
                    else:
                        dir_call = subprocess.call("mkdir {0} {1} {2}".format(*names_tup),shell=True)

                    rec_bs_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta','branch_site_model_rec', rec_mark),shell=True)
                    rec_b_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'branch_model_rec', rec_mark),shell=True)
                    site_model = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 M1 M2 M3 M7 M8 -i {3}.pdf -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'site_model'),shell=True)


                elif 'major1' in event:
                    rec_mark=mark_nodes(result['Rec_maj1'])
                    if 'unknown' not in result['Maj1_par']:
                        names_tup=tuple([os.path.join(event_subdir,x) for x in  ['site_model', 'branch_model_rec', 'branch_site_model_rec','branch_model_par', 'branch_site_model_par']])
                        par_mark=mark_nodes(result['Maj1_par'])
                        dir_call = subprocess.call("mkdir {0} {1} {2} {3} {4}".format(*names_tup),shell=True)
                        par_bs_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta','branch_site_model_par', par_mark),shell=True)
                        par_b_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'branch_model_par', par_mark),shell=True)
                    else:
                        dir_call = subprocess.call("mkdir {0} {1} {2}".format(*names_tup),shell=True)

                    rec_bs_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta','branch_site_model_rec', rec_mark),shell=True)
                    rec_b_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'branch_model_rec', rec_mark),shell=True)
                    site_model = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 M1 M2 M3 M7 M8 -i {3}.pdf -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'site_model'),shell=True)



                elif 'major2' in event:
                    rec_mark=mark_nodes(result['Rec_maj2'])
                    if 'unknown' not in result['Maj2_par']:
                        names_tup=tuple([os.path.join(event_subdir,x) for x in  ['site_model', 'branch_model_rec', 'branch_site_model_rec','branch_model_par', 'branch_site_model_par']])
                        par_mark=mark_nodes(result['Maj2_par'])
                        dir_call = subprocess.call("mkdir {0} {1} {2} {3} {4}".format(*names_tup),shell=True)
                        par_bs_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta','branch_site_model_par', par_mark),shell=True)
                        par_b_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'branch_model_par', par_mark),shell=True)
                    else:
                        dir_call = subprocess.call("mkdir {0} {1} {2}".format(*names_tup),shell=True)

                    rec_bs_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta','branch_site_model_rec', rec_mark),shell=True)
                    rec_b_call = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'branch_model_rec', rec_mark),shell=True)
                    site_model = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 M1 M2 M3 M7 M8 -i {3}.pdf -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,'ML_'+ event +'.raxml.support',event+'_nucl_reported.fasta', 'site_model'),shell=True)
   

                #dir_call = subprocess.call("mkdir {0} {1} {2}".format(*names_tup), 
                #                                  shell = True)
                #branch_site_model = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 bsA bsA1 bsB -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,event+'.nw',event+'_nucl_reported.fasta','branch_site_model', marking), 
                                                 # shell = True)
                #branch_model = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --model M0 b_free b_neut -i {3}.pdf --mark {4} -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,event+'.nw',event+'_nucl_reported.fasta', 'branch_model', marking), 
                                                  #shell = True)
                #site_model = subprocess.call('sudo docker run -v {0}:/data  --rm ete3/3.1.1 evol -t /data/{1} --alg /data/{2} --models M0 M1 M2 M3 M7 M8 -i {3}.pdf -o {3} --cpu 20 --clear_all > {0}/{3}.log'.format(event_subdir,event+'.nw',event+'_nucl_reported.fasta', 'site_model'), 
                                                  #shell = True)


    return process_row

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze selection on recombination events')
    parser.add_argument('-r', '--recomb', dest='recomb_file', help='path to the recombination table',
                        type=argparse.FileType('r'))
    parser.add_argument('-e', '--evdir', dest='evdir', help='evol_directory', type=str)
    args = parser.parse_args()
    recomb = pd.read_csv(args.recomb_file, sep='\t')
    nrecomb = recomb.apply(get_process_row(args.evdir), axis=1)

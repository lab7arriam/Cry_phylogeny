from collections import defaultdict
import argparse
import pandas as pd
from typing import List
import os.path
import subprocess
from ete3 import Tree
import random
from annotate_clusters_consistency import read_csv_to_list
import logging

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

def get_process_row(evol_dir, mod_ng, rax_ng):
    event_set=set()
    logger = logging.getLogger("evol_preparation")
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler("evol_preparation.log")
    logger.addHandler(fh)

    def process_row(row: pd.Series):
        result = row.copy()
        events=return_event_index(result)

        for event in events: 
            if event in event_set:
                break
            event_set.add(event)
            print('analyzing event ', event)
            logger.info('analyzing event ' + event)

            event_subdir = os.path.join(os.path.realpath(evol_dir),event)
            subprocess.call('{0} -i {1} -t ml --force > /dev/null'.format(os.path.realpath(mod_ng),event_subdir+'/'+event+'_nucl_reported.fasta'), 
                                                      shell = True)
            subprocess.call('cd {0}; rm *.ckp *.topos *.tree *.log'.format(event_subdir), 
                                                      shell = True)
            analyzed_list=read_csv_to_list(event_subdir+'/'+event+'_nucl_reported.fasta.out', headless=True,delim='\t')

            i=0
            for row in analyzed_list:
                if 'Best model according to BIC' in row:
                    break
                i+=1
            model = analyzed_list[i+2][0].split(':')[1].strip().replace('G4','G') 
            print('Model choosen: ', model)
            logger.info('Model choosen: ' + model)
            subprocess.call('{0} --all --msa {1} --model {2} --prefix  {3} --seed 42 --threads 3  --redo --log ERROR'.format(os.path.realpath(rax_ng),event_subdir+'/'+event+'_nucl_reported.fasta',model, event_subdir+'/ML_'+ event), 
                                                      shell = True)
            subprocess.call('cd {0}; rm -f *.bestModel *.bestTree *.bestTreeCollapsed *.bootstraps *.log *.mlTrees *.rba *.phy *.startTree'.format(event_subdir), 
                                                     shell = True)
            if sum(['support' in el for el in os.listdir(event_subdir)]) <1:
                print('tree corrupted')
                logger.info('tree corrupted ' + model)

        return result

    return process_row

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze selection on recombination events')
    parser.add_argument('-r', '--recomb', dest='recomb_file', help='path to the recombination table',
                        type=argparse.FileType('r'))
    parser.add_argument('-e', '--evdir', dest='evdir', help='evol_directory', type=str)
    parser.add_argument('-mg', '--mod_ng', dest='mod_ng', help='path to modeltest_ng file', type=str)
    parser.add_argument('-rg', '--rax_ng', dest='rax_ng', help='path to raxml_ng file', type=str)
    args = parser.parse_args()
    recomb = pd.read_csv(args.recomb_file, sep='\t')
    nrecomb = recomb.apply(get_process_row(args.evdir, args.mod_ng, args.rax_ng), axis=1)


from collections import defaultdict
from ete3 import Tree
import argparse
import pandas as pd
from typing import List
import random
import os.path
import subprocess
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from annotate_clusters_consistency import write_csv, read_csv_to_list
from statistics import mean


rename_dict={"NOV_Cry2Ah2_92.7_AAQ52362.1":"NOV_Cry2Ah2_92.7_AX098631.1", "NOV_Cry40Aa1_67.9_ACC34441.1":"NOV_Cry40Aa1_67.9_CQ868314.1", "NOV_Cry8Bb1_98.4_ABJ38812.1":"NOV_Cry8Bb1_98.4_CS131028.1", "NOV_Cry2Ak1_76.8_ABU37568.1":"NOV_Cry2Ak1_76.8_AX513524.1", "NOV_Cry1Ga1_77.9_AAQ52381.1":"NOV_Cry1Ga1_77.9_AX098669.1","NOV_Cry1Ab18_92.4_AAE49246.1":"NOV_Cry1Ab18_92.4_AX383797.1", "BT_Cry1Jc1":"NOV_Cry1Jc1_99.1_AAS93799.1", "NOV_Cry24Aa1_98.3_WP_086403584.1":"NOV_Cry24Aa1_98.3_WP_086403584.1" }

def gca(tree: Tree, nodes: List[str]) -> Tree:
    """
        Reads a ete3 Tree object and list of leaves to find the common ancestor
        Returns a subtree representing the last common ancestor of a specified list  
    """
    no = [tree & n for n in nodes]
    if len(no) == 1:
        return no[0]
    return tree.get_common_ancestor(no)

def select_tree(row: pd.Series, trees: List[Tree]) -> Tree:
    """
        Returns a tree based on a domain swapped during a recombination event 
    """
    rtype = str(row['Type'])
    if rtype.startswith('domain1'):
        return trees[0]
    elif rtype.startswith('domain2'):
        return trees[1]
    elif rtype.startswith('domain3'):
        return trees[2]
    else:
        return None

def select_oppos_tree(row: pd.Series, trees: List[Tree]) -> list:
    """
        Returns a list with trees for domains not involved in recombination
    """
    rtype = str(row['Type'])
    if rtype.startswith('domain1'):
        return [trees[1],trees[2]]
    elif rtype.startswith('domain2'):
        return [trees[0],trees[2]]
    elif rtype.startswith('domain3'):
        return [trees[0],trees[1]]
    else:
        return None

def return_domain_names(row: pd.Series, mode) -> list:
    """
        Returns domain names for a recombination event either for domains exposed to recombination or the intact ones
    """
    rtype = str(row['Type'])
    if mode == 'str':
        if rtype.startswith('domain1'):
            return ['domain1']
        elif rtype.startswith('domain2'):
            return ['domain2']
        elif rtype.startswith('domain3'):
            return ['domain3']
        else:
            return None
    else:
        if rtype.startswith('domain1'):
            return ['domain2', 'domain3']
        elif rtype.startswith('domain2'):
            return ['domain1','domain3']
        elif rtype.startswith('domain3'):
            return ['domain1','domain2']
        else:
            return None


def get_node_len(node: Tree) -> int:
    """
        Reads a ete3 Tree object and returns the number of leaves 
    """
    return(len([leaf.name for leaf in node.iter_leaves()]))

def get_node_names(node: Tree) -> list:
    """
        Reads a ete3 Tree object and returns the list with respective names of leaves
    """
    return([leaf.name for leaf in node.iter_leaves()])



def make_clusters_dict(clust_rows: list) -> defaultdict(dict):
    """
        Reads parsed rows of clustered file and creates a dictionary with sequence names for each cluster 
    """
    clust_dict=defaultdict(dict)

    for row in clust_rows:
        clust_dict[row[0]][row[1]]=row[2].split(';')

    return(clust_dict)


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
            if len(clust_dict[toxin][domain])<=4:
                return(clust_dict[toxin][domain])
            else:
                return(list(random.sample(clust_dict[toxin][domain],4)))

def extend_node_with_cluster(node_list: list, domain: str, clust_dict: defaultdict(dict), flag_group=None) -> list:
    """
       Reads a list of toxins and returns extended list with disclosed clusters
    """   
    ret_list=node_list

    for tox in node_list:
        if flag_group:

            if tox in flag_group:

                disc_list=disclose_cluster(tox,domain, clust_dict) 

                if len(disc_list)>0:

                    #flag_group=list(set(disc_list+flag_group))
                    flag_group.extend(disc_list)
                    ret_list.extend(disc_list)
        else:
            disc_list=disclose_cluster(tox,domain, clust_dict) 

            if len(disc_list)>0:
                ret_list.extend(disc_list)   

    return(list(set(ret_list)))

def upper_node_search(in_node: Tree, in_tree: Tree, len_lim: int, domain: str,clust_dict: defaultdict(dict)):
    """
       Reads a subtree and a reference tree, goes to the upper levels of the respective tree until exceeding a threshold specified
       Returns a list with names of the leaves with unfolded clusters
    """   
    tree=in_tree.copy("deepcopy")
    node=gca(tree, get_node_names(in_node))
    nodes_list=[]
    parsed_node=extend_node_with_cluster(get_node_names(node), domain, clust_dict)

    if len(parsed_node)>=len_lim:
        return parsed_node 

    else:
        while len(extend_node_with_cluster(get_node_names(node), domain, clust_dict))<=len_lim:
            node=node.up
            nodes_list.append(node)

        if len(nodes_list)==1:
            return(extend_node_with_cluster(get_node_names(nodes_list[0]), domain, clust_dict))
        else:
            return(extend_node_with_cluster(get_node_names(nodes_list[-1]), domain, clust_dict))

def get_identity_heatmap(map_file: str ) -> defaultdict(dict):
    """
    Reads paths to the heatmap file and returns a dictionary with identities for each pair of toxins
    """
    heatmap_dict=defaultdict(dict)
    for row in map_file:

        if row[0] in rename_dict:
            row[0]=rename_dict[row[0]]
        if row[1] in rename_dict:
            row[1]=rename_dict[row[1]]

        heatmap_dict['|'.join(sorted([row[0],row[1]]))]['domain1']=float(row[5])
        heatmap_dict['|'.join(sorted([row[0],row[1]]))]['domain2']=float(row[11])
        heatmap_dict['|'.join(sorted([row[0],row[1]]))]['domain3']=float(row[17])

    return(heatmap_dict)

def check_sets(set1,set2):
    """
    Calculates intersection between two sets and returns False if the intersection is not empty
    """
    if len(set1.intersection(set2))>0:
        return False
    else:
        return True


def tree_iterate(in_node: Tree, in_tree: Tree, len_lim: int, domain: str,clust_dict: defaultdict(dict), heatmap_dict: defaultdict(dict), recs: list):
    """
    Iterates over a tree and creates a subtree for a node specified, adds clades within a tree based on mean sequence identity with recombinants
    Returns a list with toxin names
    """
    tree=in_tree.copy("deepcopy")
    node=gca(tree, [leaf.name for leaf in in_node.iter_leaves()])
    parsed_node=list(set(extend_node_with_cluster(get_node_names(node), domain, clust_dict)))

    #set subtree threshold
    if len(parsed_node)<len_lim:
        len_limit=len_lim
    else:
        len_limit=len(node)*2

    #exctract clades from the tree for each leaf
    parsed_leaves=[]
    for leaf in tree:
        while leaf:
            if len(leaf)>1 and len(leaf)<=5 and check_sets(set(get_node_names(leaf)),set(get_node_names(node))) and get_node_names(leaf) not in parsed_leaves:
                parsed_leaves.append(get_node_names(leaf))
            leaf = leaf.up

    #make a dictionary based on mean sequence identity with recombinants for each clade
    parsed_dict=defaultdict(list)
    node_count=0
    for p_leaf in parsed_leaves:
        passed=list()
        ids=list()
        for p_tox in p_leaf:
            for rec in recs:
                if rec!=p_tox:
                    pair='|'.join(sorted([p_tox,rec]))
                    if pair not in passed:
                        ids.append(float(heatmap_dict[pair][domain]))
                        passed.append(pair)
        parsed_dict[mean(ids)+0.04*node_count]=p_leaf
        node_count+=1
     
    #add subclades to the initial nodes until threshold is reached
    for parse_key in sorted(parsed_dict.keys(), reverse=True):

        if len(parsed_node)>=len_limit:
            break

        parsed_node=list(set(parsed_node + extend_node_with_cluster(parsed_dict[parse_key], domain, clust_dict)))

    return(parsed_node)



def make_evol_trees(rec: list, rec_mask: list, par: list, row: pd.Series, tree: Tree, result: pd.Series, domain: str, clust_dict: defaultdict(dict), heatmap_dict: defaultdict(dict), tree_dir: str, prot_seqs: defaultdict(list), nucl_seqs: defaultdict(list), event_pref: str, unknown=False):
    """
    
    """
    nrec=rec.copy()
    npar=par.copy()
    #print(len(nrec))
    if par[0] == 'none':
        npar=[]
        upper_par_node=[]

    elif par[0]!='none':
        if unknown:
            npar=[]
        if par[0] == 'unknown':
            npar=[]
            upper_par_node=[]
            
        #calculate the initial lentgh of a node
        par_node=gca(tree, rec + npar)
        par_thr_list=get_node_names(par_node)
        par_thr_list=extend_node_with_cluster(par_thr_list, domain, clust_dict, nrec)
        par_thr_list=extend_node_with_cluster(par_thr_list, domain, clust_dict, npar)

        #set a threshold for a subtreee
            
        if len(par_thr_list)<7:
            len_thr = 7

        elif len(par_thr_list)<100:

            len_thr = len(par_thr_list)*2

        else:

            len_thr = len(par_thr_list) + 1
            
        if len(par_thr_list) <=60:
            #go to the upper level 
            upper_par_node=upper_node_search(par_node,tree, len_thr, domain, clust_dict)

            if len(upper_par_node)<=60:
                pass

            #create a tree from sets of leaves
            else:
                upper_par_node=tree_iterate(par_node, tree, len_thr,  domain, clust_dict, heatmap_dict, rec_mask)
        else:
            upper_par_node=get_node_names(par_node)

    if len(upper_par_node) > 5 and len(upper_par_node)<70 and  'partial' not in row['Type']:
        name=os.path.join(os.path.realpath(tree_dir),event_pref +str(row['ID']))
        subprocess.call('mkdir {}'.format(name),shell=True)
        create_fasta_files(tree_dir, event_pref+str(row['ID']), domain, upper_par_node, prot_seqs, nucl_seqs)
    if par[0] == 'none':
        npar=['none']
    if par[0] == 'unknown':
        npar=['unknown']
    if unknown:
        npar=['unknown']
    
    return(upper_par_node, list(set(nrec)), list(set(npar)))
    
def create_fasta_files(tree_dir, event, domain, cry_list, prot_seqs, nucl_seqs):
    p_name=os.path.join(os.path.realpath(tree_dir),event)+'/' +event +'_prot.fasta'
    p_seqs=[cry for cry in prot_seqs[domain] if cry.id in cry_list]
    SeqIO.write(p_seqs,p_name,"fasta")
    n_name=os.path.join(os.path.realpath(tree_dir),event)+'/' +event +'_nucl.fasta'
    n_seqs=[cry for cry in nucl_seqs[domain] if cry.id in cry_list]
    SeqIO.write(n_seqs,n_name,"fasta")
    mafft_call = subprocess.call("mafft --quiet --maxiterate 1000 --localpair --op 3.53 --thread 8  {0} > {1} ".format(p_name,os.path.join(os.path.realpath(tree_dir),event)+'/' +event +'_msa.fasta'), 
                                                      shell = True)
    report_gaps = subprocess.call("python3 report_gaps.py {0} {1}".format(os.path.join(os.path.realpath(tree_dir),event)+'/' +event +'_msa.fasta', n_name), 
                                                      shell = True)
    clean_call = subprocess.call("rm  {0} {1} {2}".format(p_name, n_name,os.path.join(os.path.realpath(tree_dir),event)+'/' +event +'_msa.fasta'), 
                                                      shell = True)


def get_parent_names(proc_str):
    try:
        rec_str = proc_str.split('(')[1].replace(')','').split(';')
    except:
        rec_str = proc_str.replace(')','').split(';')
    return(rec_str)

def get_process_row(trees: List[Tree], tree_dir: str, prot_seqs: defaultdict(list), nucl_seqs: defaultdict(list),clust_dict: defaultdict(dict),heatmap_dict: defaultdict(dict)):
    def process_row(row: pd.Series):
        rec = row['Rec'].split(';')
        rec = [r.strip() for r in rec]
        row['Type']=row['Type'].replace('_partial','')
        if row['Type']=='domain1':
            min=get_parent_names(row['parent_1'])
            min = [m.strip() for m in min]

            maj1=get_parent_names(row['parent_2'])
            maj1 = [m.strip() for m in maj1]

            maj2=get_parent_names(row['parent_3'])
            maj2 = [m.strip() for m in maj2]

        if row['Type']=='domain2':
            min=get_parent_names(row['parent_2'])
            min = [m.strip() for m in min]

            maj1=get_parent_names(row['parent_1'])
            maj1 = [m.strip() for m in maj1]

            maj2=get_parent_names(row['parent_3'])
            maj2 = [m.strip() for m in maj2]

        if row['Type']=='domain3':
            min=get_parent_names(row['parent_3'])
            min = [m.strip() for m in min]

            maj1=get_parent_names(row['parent_1'])
            maj1 = [m.strip() for m in maj1]

            maj2=get_parent_names(row['parent_2'])
            maj2 = [m.strip() for m in maj2]

        #old paradigm for minor/major events
        #min = row['Min_par'].split(';')
        #min = [m.strip() for m in min]
        #maj = row['Maj_par'].split(';')
        #maj = [m.strip() for m in maj]
        
        print(row['ID'],'________________________________________________________')
        rec_mask, min_mask, maj1_mask, maj2_mask = rec.copy(), min.copy(), maj1.copy(), maj2.copy()
        tree = select_tree(row, trees)
        op_tree = select_oppos_tree(row, trees)
        result = row.copy()

        #get lengths of trees
        rec_node=get_node_len(gca(tree, rec))
        try:
            min_node=get_node_len(gca(tree, rec+min))
        except:
            min_node=0
        try:
            maj1_node=get_node_len(gca(op_tree[0], rec+maj1))
        except:
            maj1_node=0
        try:
            maj2_node=get_node_len(gca(op_tree[1], rec+maj2))
        except:
            maj2_node=0
       
        if row['Unknown_flag']=='minor':
            unknown=True
        else:
            unknown=False
        upper_min_res=make_evol_trees(rec, rec_mask, min, row, tree, result, return_domain_names(row, 'str')[0], clust_dict, heatmap_dict, tree_dir, prot_seqs, nucl_seqs, 'minor_event_', unknown)

        if row['Unknown_flag']=='major':
            unknown=True
        else:
            unknown=False
        upper_maj_res1=make_evol_trees(rec, rec_mask, maj1, row, op_tree[0], result, return_domain_names(row, 'no')[0], clust_dict, heatmap_dict, tree_dir, prot_seqs, nucl_seqs,'major1_event_', unknown)
        upper_maj_res2=make_evol_trees(rec, rec_mask, maj2, row, op_tree[1], result, return_domain_names(row, 'no')[1], clust_dict, heatmap_dict, tree_dir, prot_seqs, nucl_seqs,'major2_event_', unknown)

        rec=upper_min_res[1]
        min=upper_min_res[2]
        maj1=upper_maj_res1[2]
        maj2=upper_maj_res2[2]

        result=result.iloc[[0]]
        result['Rec_min']=';'.join(sorted(rec))
        result['Min_par']=';'.join(sorted(min))    
        result['Maj1_par']=';'.join(sorted(maj1))
        result['Rec_maj1']=';'.join(sorted(upper_maj_res1[1]))
        result['Maj2_par']=';'.join(sorted(maj2))
        result['Rec_maj2']=';'.join(sorted(upper_maj_res2[1]))
        result['Type']=row['Type']
        result['Unknown_flag']= row['Unknown_flag']      
        
        for x, y in zip(['Rec_len','Min_len','Maj1_len','Maj2_len','Min_evol', 'Maj1_evol', 'Maj2_evol'],[rec_node,min_node, maj1_node, maj2_node,len(upper_min_res[0]),len(upper_maj_res1[0]),len(upper_maj_res2[0])]):
            result[x]=y
        return result

    return process_row


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create subdirectories for further analysis of evolutionary selection')
    parser.add_argument('-r', '--recomb', dest='recomb_file', help='path to the recombination table',
                        type=argparse.FileType('r'))
    parser.add_argument('-t', '--tree', dest='tree_file', help='trees for domains', type=argparse.FileType('r'),
                        action='append')
    parser.add_argument('-o', '--ofile', dest='ofile', help='output file', type=argparse.FileType('w'))
    parser.add_argument('-e', '--evdir', dest='evdir', help='evol_directory', type=str)
    parser.add_argument('-n', '--nucl', dest='nucl_seqs', help='nucl sequences', type=argparse.FileType('r'),
                        action='append')
    parser.add_argument('-p', '--prot', dest='prot_seqs', help='protein sequences', type=argparse.FileType('r'),
                        action='append')
    parser.add_argument('-c', '--clust', dest='clust_dat', help='table with reference clusters', type=str)
    parser.add_argument('-hf', '--heat_f', dest='heat_file', help='heatmap with paiwise identity of toxins', type=str)


    args = parser.parse_args()

    trees = []
    for tf in args.tree_file:
        t = Tree(tf.name)
        trees.append(t)

    nucl_seqs=defaultdict(list)
    prot_seqs=defaultdict(list)

    i=1
    for nuc in args.nucl_seqs:
        nucl_seqs['domain' + str(i)]=list(SeqIO.parse(nuc, "fasta"))
        i+=1

    i=1
    for prot in args.prot_seqs:
        prot_seqs['domain' + str(i)]=list(SeqIO.parse(prot, "fasta"))
        i+=1

    clust_dict = make_clusters_dict(read_csv_to_list(args.clust_dat, headless=True,delim='\t'))
    heatmap_dict=get_identity_heatmap(read_csv_to_list(args.heat_file, headless=True,delim='\t'))
    recomb = pd.read_csv(args.recomb_file, sep='\t')
    #recomb = recomb.drop(recomb.columns[0], axis=1)
    nrecomb = recomb.apply(get_process_row(trees,args.evdir,prot_seqs,nucl_seqs,clust_dict,heatmap_dict), axis=1)
    
    nrecomb.to_csv(args.ofile.name, sep='\t', index=False)

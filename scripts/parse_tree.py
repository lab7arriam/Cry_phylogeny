import sys
import csv
from collections import defaultdict
import re

tree_str=''
with open (sys.argv[1], 'r',newline='') as csvfile2:
    my_reader3 = csv.reader(csvfile2, delimiter='\t')
    for row in my_reader3:
        tree_str=row[0]

start, stop, node = 0, 0, 0

def remove_last(stack):
    return(stack[0:len(stack)-1])

def get_last(stack):
    return(stack[len(stack)-1])

node_dict=defaultdict(list)
close_node_dict=defaultdict(list)

def prepare_tree(tree_str,close_node_dict):
    replaced_tree=str()
    replaced_tree=tree_str
    ind_stack=list()
    pairs=list()
    for i in range(len(tree_str)):
        if tree_str[i]=='(':
            ind_stack.append(i)
        if tree_str[i]==')':
            pairs.append([get_last(ind_stack),i])
            ind_stack=remove_last(ind_stack)
    pair_lens=sorted(list(set(sorted([int(pair[1])-int(pair[0]) for pair in pairs]))))
    i=len(close_node_dict)+1
    for leng in pair_lens:
        name_set=set()
        for pair in pairs:
            if int(pair[1])-int(pair[0]) == leng:
                if tree_str[int(pair[0]):int(pair[1])+1].count('(')==1:
                    if tree_str[int(pair[0]):int(pair[1])+1].replace(')','').replace('(','') not in name_set:
                        close_node_dict[i].append(tree_str[int(pair[0]):int(pair[1])+1].replace(')','').replace('(',''))
                        replaced_tree=replaced_tree.replace(tree_str[int(pair[0]):int(pair[1])+1],str(i))
                        i+=1
                    name_set.add(tree_str[int(pair[0]):int(pair[1])+1].replace(')','').replace('(',''))
    return(replaced_tree)

for iter_num in range(50):
    tree_str=prepare_tree(tree_str,close_node_dict)


def check_cry(cry_list):
    flag=0
    for el in cry_list:
        if 'Cry' not in el:
            flag=1
    if flag==0:
        return(True)
    else:
        return(False)



with open ('col_3_nucl_branch_stat_extended.tsv', 'w',newline='') as csvfile2:
    my_writer = csv.writer(csvfile2, delimiter='\t')
    my_writer.writerow(['depth','cry_list','cry_num','node','init_structure'])

    for node in close_node_dict:
        #print(node, close_node_dict)
        node_str=[]
        cry_list=close_node_dict[node][0].split(',')
        #print('list',cry_list)
        for i in range(1,50):
            if check_cry(cry_list):
                break
            else:
                passed=list()
                new_list=list()
                for el in cry_list:
                    if 'Cry' not in el:
                        node_str.append(el)
                        new_el=close_node_dict[int(el)][0].split(',')
                        new_list.extend(new_el)

                    else:
                        new_list.append(el)
                cry_list=new_list
            #print(i,cry_list)
        print(node_str)
        my_writer.writerow(["#"+str(i),','.join(cry_list),len(cry_list),node,'; '.join(close_node_dict[node][0].split(','))])
    #print(node, "#"+str(i),cry_list,close_node_dict[node][0].split(','))
    #print(close_node_dict[node],check_cry(close_node_dict[node][0].split(',')))
	



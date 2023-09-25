#/usr/bin/python3.7

import csv
import argparse
from statistics import mean


def get_mean_support_for_tree(tree_file, percent=True, sup_type='r'):
    """
        Calculates mean support for the tree in newick format
        If type is regular (r, by deafult), integer support numbers are expected, otherwise parses tree according to the expected posision of the supporting values
        If percent flag is enabled calculates raw values, otherwize multiplies by 100
    """

    tree_str=''
    with open (tree_file, 'r',newline='') as csvfile2:
        my_reader3 = csv.reader(csvfile2, delimiter='\t')
        for row in my_reader3:
            tree_str=row[0]
    if sup_type=='r':
        tree_list=tree_str.replace('(',',').replace(':',',').replace(')',',').split(',')
        parsed_list=[]

        for el in tree_list:
            if '_' not in el and '.' not in el and len(el)>0 and ';' not in el:
                try:
                    parsed_list.append(int(el)) 
                except:
                    pass


        if percent:
            return(mean(parsed_list))
        else:
            return(mean([el*100 for el in parsed_list]))
    else:
        parsed_list=[]
        tree_list=tree_str.replace('(','').replace(',',')').split(')')
        for el in tree_list:
            if '_' not in el and len(el)>0 and ';' not in el and el.split(':')[0]!='0':
                parsed_list.append(float(el.split(':')[0]))
        if percent:
            return(mean(parsed_list))
        else:
            return(mean([el*100 for el in parsed_list]))

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculates mean supporting value for the tree in the newick format')
    parser.add_argument('-t', '--tr_f', dest='tree_file', help='path to the analyzed tree file',
                        type=str)
    parser.add_argument('-p', '--per', dest='per_flag', type=str2bool, nargs='?',
                        const=True, default=False,
                        help="type of analysis (show percentage if True)")
    args = parser.parse_args()
    try:
        print(get_mean_support_for_tree(args.tree_file, args.per_flag, sup_type='i'))
    except:
        print(get_mean_support_for_tree(args.tree_file, args.per_flag, sup_type='r'))

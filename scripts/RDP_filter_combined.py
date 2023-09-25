import sys
import csv
from collections import defaultdict
import os
import os.path
import subprocess 
from statistics import mean
from Bio import SeqIO
import argparse
from annotate_clusters_consistency import write_csv, read_csv_to_list, create_dict_from_list


rename_dict={"NOV_Cry2Ah2_92.7_AAQ52362.1":"NOV_Cry2Ah2_92.7_AX098631.1", "NOV_Cry40Aa1_67.9_ACC34441.1":"NOV_Cry40Aa1_67.9_CQ868314.1", "NOV_Cry8Bb1_98.4_ABJ38812.1":"NOV_Cry8Bb1_98.4_CS131028.1", "NOV_Cry2Ak1_76.8_ABU37568.1":"NOV_Cry2Ak1_76.8_AX513524.1", "NOV_Cry1Ga1_77.9_AAQ52381.1":"NOV_Cry1Ga1_77.9_AX098669.1","NOV_Cry1Ab18_92.4_AAE49246.1":"NOV_Cry1Ab18_92.4_AX383797.1", "BT_Cry1Jc1":"NOV_Cry1Jc1_99.1_AAS93799.1", "NOV_Cry24Aa1_98.3_WP_086403584.1":"NOV_Cry24Aa1_98.3_WP_086403584.1" }


def clean_name(name_str: str) -> str:
    """
     Cleans RDP generic name by removing special symbols and extra spaces, returns cleaned name
    """
    return(name_str.replace('^','').replace('~','').replace('*','').replace('Unknown (','').replace(')','').replace('[P]','').replace('Unknown(','').replace('[T]',''))

def getOverlap(a, b):
    """
    Calculates the length of an intersected interval for two intervals
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def init_RDP_filter(rdp_file : list) -> list:
    """
    Reads RDP-generated table and filters events marked with swing dash
    Returns a list with parsed rows for further annotaion
    """
    ret_list=[]
    event_dict=defaultdict(dict)
    tilda_names=set()

    for row in rdp_file:

        if len(row)>1:

            if len(row[0])>0:

                key=row[0].replace(' ','')
 
                if '~' in row[1]:
                    tilda_names.add(key)

                if len(row)==21:

                    event_dict[key]['recs']=[clean_name(row[8])]
                    event_dict[key]['pvals']=row[11:20]
                    event_dict[key]['coords']=sorted([int(clean_name(x)) for x in  [row[5],row[4]]])
                    event_dict[key]['major_par']=[clean_name(row[10])]
                    event_dict[key]['minor_par']=[clean_name(row[9])]

                    if 'Unknown' in row[10]:    
                        event_dict[key]['unknown']='major'
                    if 'Unknown' in row[9]:
                        event_dict[key]['unknown']='minor'
                    if 'Unknown' not in row[9] and 'Unknown' not in row[10]:
                        event_dict[key]['unknown']='no'

                elif len(row)==11:

                    event_dict[key]['recs'].append(clean_name(row[8]))
                    event_dict[key]['major_par'].append(clean_name(row[10]))
                    event_dict[key]['minor_par'].append(clean_name(row[9]))


    for event in event_dict:
        if event not in tilda_names: 
            recs = ';'.join([rec for rec in event_dict[event]['recs'] if len(rec)>0])
            min_p = ';'.join([rec for rec in event_dict[event]['minor_par'] if len(rec)>0])
            maj_p = ';'.join([rec for rec in event_dict[event]['major_par'] if len(rec)>0])
            ret_list.append([recs,min_p,maj_p,event_dict[event]['unknown']]+event_dict[event]['coords']+[min([float(pval) for pval in event_dict[event]['pvals'] if pval!='NS'])])

    return(ret_list)



def read_annotation_files(tree_parsed_dir : str, sim_file : str , map_file : str ):
    """
    Reads paths to files with annotation infomation and creates respective global dictionaries
    """

    global parsed_tree_dict 
    global similarity_dict
    global dom_mappings_dict
    
    parsed_tree_dict = defaultdict(list)
    similarity_dict=defaultdict(list)
    dom_mappings_dict=defaultdict(list)

    prep_list = ['full_nucl','1_nucl','2_nucl','3_nucl'] 

    for prep in prep_list:

        parsed_tree_dict[prep]=list()
        rel_name = os.path.join(tree_parsed_dir,prep+'_branch_stat_extended.tsv')

        for row in read_csv_to_list(rel_name, headless=True,delim='\t'):
            parsed_tree_dict[prep].append(row[0:2])

    for row in read_csv_to_list(sim_file, headless=True,delim='\t'):

        if row[0] in rename_dict:
            row[0]=rename_dict[row[0]]
        if row[1] in rename_dict:
            row[1]=rename_dict[row[1]]

        similarity_dict['|'.join(sorted([row[0],row[1]]))]=[row[5],row[11], row[17]]
    
    for row in read_csv_to_list(map_file, headless=True,delim='\t'):
        dom_mappings_dict[row[0]]=[int(x)*3+1 for x in row[1:len(row)]]


def find_dom_overlaps(start:int, stop:int, map_list:list, filt_mode=False) -> list:
    """
    Calculates overlaps between recombination breakpoints and domain coordinates
    Returs a list of overlaps for 1,2, and 3 domains
    """
    over1=getOverlap([start,stop],map_list[0:2])/(int(map_list[1])-int(map_list[0]))
    over2=getOverlap([start,stop],map_list[1:3])/(int(map_list[2])-int(map_list[1]))
    over3=getOverlap([start,stop],map_list[2:4])/(int(map_list[3])-int(map_list[2]))

    if filt_mode:
        if over1>over2 and over1>over3 and over1>0.6 and over3<0.3 and over2<0.3:
            dom='domain1'
        elif over2>over1 and over2>over3 and over2>0.6 and over3<0.3 and over1<0.3:
            dom='domain2'
        elif over3>over2 and over3>over1 and over3>0.6 and over1<0.3 and over2<0.3:
            dom='domain3'
        else:
            dom='other'

    else:
        if over1>over2 and over1>over3:
            dom='domain1'
        elif over2>over1 and over2>over3:
            dom='domain2'
        elif over3>over2 and over3>over1:
            dom='domain3'

    if max([over1,over2,over3])<0.70:
        dom+='_partial'

    return([over1,over2,over3,dom])


def reformat_coords_for_manhattan(overlap_list: list) -> list:
    """
    Reads domain mappings results and creates coordinates for the Manhattan plot
    """
    if 'domain3' in overlap_list[3]:
        if overlap_list[2]==1:
            ret_list=[[(1-overlap_list[1])*100+100,overlap_list[4],'domain3'],[298,overlap_list[4],'domain3']]
        else:
            if sum([x>0 for x in overlap_list[0:3]])>1:
                ret_list=[[(1-overlap_list[1])*100+100,overlap_list[4],'domain3'],[overlap_list[2]*100+200,overlap_list[4],'domain3']]
            else:  
                ret_list=[[overlap_list[2]*100+200,overlap_list[4],'domain3'],[199,overlap_list[4],'domain3']]

    elif 'domain2' in overlap_list[3]:
        if overlap_list[0]==0 and overlap_list[2]==0:
            ret_list=[[100,overlap_list[4],'domain2'],[overlap_list[1]*100+100,overlap_list[4],'domain2']]
        elif overlap_list[1]==1:
            ret_list=[[(1-overlap_list[0])*100,overlap_list[4],'domain2'],[overlap_list[2]*100+200,overlap_list[4],'domain2']]
        elif overlap_list[0]>0 and overlap_list[1]>0 and overlap_list[2]==0:
            ret_list=[[(1-overlap_list[0])*100,overlap_list[4],'domain2'],[overlap_list[1]*100+100,overlap_list[4],'domain2']]
        else:
            ret_list=[[(1-overlap_list[1])*100+100,overlap_list[4],'domain2'],[overlap_list[2]*100+200,overlap_list[4],'domain2']]

    else:
        if overlap_list[0]==1:
            ret_list=[[overlap_list[1]*100+100,overlap_list[4],'domain1'],[2,overlap_list[4],'domain1']]
        else:
            if sum([x>0 for x in overlap_list[0:3]])>1:
                ret_list=[[(1-overlap_list[0])*100,overlap_list[4],'domain1'],[overlap_list[1]*100+100,overlap_list[4],'domain1']]
            else:   
                ret_list=[[overlap_list[0]*100,overlap_list[4],'domain1'],[1,overlap_list[4],'domain1']]
    return(ret_list)


def search_tree(input_list : list,tree_parse_list : list) -> list:
    """
    Finds a subtree encompassing all the items specified in a list
    Returns a list of leaf names for a subtree
    """
    cry_list=[]
    for branch in tree_parse_list:
        if sum([x in branch[1].split(',') for x in input_list])==len(input_list):
            cry_list=branch[1].split(',')
            break
    return(cry_list)   


def inter_group_pair(group1 : list,group2 : list, similarity_dict : defaultdict(list), pair_flag=False) -> list:
    """
    Calculates mean identity for groups of toxins in each domain
    Returns a list with respective mean similarities or the maximal similarity pair
    """

    dom1=list()
    dom1_names=list()
    dom2=list()
    dom2_names=list()
    dom3=list()
    dom3_names=list()
    passed=list()

    for cry1 in group1.replace(' ','').split(';'):
        for cry2 in group2.replace(' ','').split(';'):
            if cry1!=cry2:
                pair='|'.join(sorted([cry1,cry2]))

                if pair not in passed:
                    dom1.append(float(similarity_dict[pair][0]))
                    dom1_names.append(pair)
                    dom2.append(float(similarity_dict[pair][1]))
                    dom2_names.append(pair)
                    dom3.append(float(similarity_dict[pair][2]))
                    dom3_names.append(pair)
                    passed.append(pair)

    if len(dom1)==0:
        dom1.append(100)
    if len(dom2)==0:
        dom2.append(100)
    if len(dom3)==0:
        dom3.append(100)

    if not pair_flag:
        return([mean(dom1+dom2+dom3),mean(dom1), mean(dom2), mean(dom3)])

    else:
        return([str(round(max(dom1),2))+ ' (' +dom1_names[dom1.index(max(dom1))]+ ')',str(round(max(dom2),2))+ ' (' +dom2_names[dom2.index(max(dom2))]+ ')',str(round(max(dom3),2))+ ' (' +dom3_names[dom3.index(max(dom3))]+ ')'])


def getlistpercent(a, b):
    """
    Calculates the percentage of elements from one list in another list
    """
    return sum([x in b for x in a])/len(b)

def RDP_annotation(rdp_file : list, parsed_tree_dict : defaultdict(list), similarity_dict : defaultdict(list)) -> list:
    """
    Reads RDP table and the data for annotation (parsed tree tables and similarity )
    Returns a list with annotated RDP events
    """
    annot_list=[]
    for row in rdp_file:

        tree_percents=list()
        similarity_inf=list()

        for cry_list in row[0:3]:
            for tree_list in parsed_tree_dict.keys():
                if len(cry_list.split(';'))==1:
                    tree_percents.append(str(100))
                else:
                    tree_percents.append(str(round(getlistpercent(row[0].split(';'),search_tree(row[0].split(';'), parsed_tree_dict[tree_list]))*100,2)))
            similarity_inf.extend(inter_group_pair(cry_list,cry_list,similarity_dict))

        for ind_pair in [(0,1),(0,2),(1,2)]:
            similarity_inf.extend(inter_group_pair(row[ind_pair[0]],row[ind_pair[1]],similarity_dict) + inter_group_pair(row[ind_pair[0]],row[ind_pair[1]],similarity_dict,pair_flag=True))

        annot_list.append(row[0:14] + [';'.join(tree_percents[0:4]),';'.join(tree_percents[4:8]),';'.join(tree_percents[8:12])] + similarity_inf)

    return(annot_list)
    
        
def RDP_filter_pipeline(rdp_file : list, tree_parsed_dir : str, sim_file : str , map_file : str, filt_mode=False) -> list:
    """
    Reads raw RDP table and the data for annotation
    Returns an filtered list with parsed RDP events
    """
    raw_rdp_file = init_RDP_filter(rdp_file)
    read_annotation_files(tree_parsed_dir,sim_file,map_file)

    filt_RDP_file = list()
    manhat_file = list()
    for row in raw_rdp_file:

        rec_key=row[0].split(';')[0]
       
        map_bonds=[1,dom_mappings_dict[rec_key][2]] +dom_mappings_dict[rec_key][4:]
        overlaps_filt=find_dom_overlaps(row[4],row[5], map_bonds, filt_mode=True)

        if 'other' in overlaps_filt[3]:
            over1=find_dom_overlaps(1,row[4], map_bonds, filt_mode=True)
            over3=find_dom_overlaps(row[5],dom_mappings_dict[rec_key][5], map_bonds, filt_mode=True)
            if over1[0]>0.6 and  over1[1]<0.2 and over1[2]<0.2:
                row[5]=row[4]
                row[4]=1
            elif over3[2]>0.6 and overlaps_filt[0]>0.6 and overlaps_filt[1]>0.6 and overlaps_filt[2]<0.2:
                row[4]=row[5]
                row[5]=dom_mappings_dict[rec_key][5]
            elif (overlaps_filt[1]>0.7 and overlaps_filt[2]>0.7 and overlaps_filt[0]<0.5):
                row[5]=row[4]
                row[4]=1
            elif over3[2]>0.3 and overlaps_filt[0]>0.7 and overlaps_filt[1]>0.7:
                row[4]=row[5]
                row[5]=dom_mappings_dict[rec_key][5]

        overlaps=find_dom_overlaps(row[4],row[5], map_bonds)
        if max(overlaps[0:3])>0.35:
            filt_RDP_file.append(row[0:4]+[overlaps[3]]+row[4:6]+dom_mappings_dict[rec_key] + [row[6]])
            if 'partial' not in overlaps[3]:
                manhat_file.extend(reformat_coords_for_manhattan(overlaps+[row[6]]))

    annot_RDP=RDP_annotation(filt_RDP_file, parsed_tree_dict,similarity_dict)
    write_csv(manhat_file,'recomb_manhattan_filtered.tsv')
    write_csv(annot_RDP,'RDP_pre_filtered_fixed_unpivot.tsv','Rec','Min_par','Maj_par','Unknown_flag','Type','start','stop','d1_start','d1_stop','d2_start','d2_stop','d3_start','d3_stop','Pval', 
           'rec_tree_cons:f|1|2|3','min_tree_cons','maj_tree_cons','rec_a','rec_1','rec_2','rec_3','min_a','min_1','min_2','min_3','maj_a','maj_1','maj_2','maj_3',
           'r_min_a','r_min_1','r_min_2','r_min_3','r_min_1_pair','r_min_2_pair','r_min_3_pair','r_maj_a','r_maj_1','r_maj_2','r_maj_3','r_maj_1_pair','r_maj_2_pair','r_maj_3_pair',
           'max_min_a','max_min_1','max_min_2','max_min_3','max_min_1_pair','max_min_2_pair','max_min_3_pair')
    
def make_manhattan_only(rdp_file : list, tree_parsed_dir : str, sim_file : str , map_file : str, filt_mode=False,new_annot=True) -> list:
    """
    Reads raw RDP table and the data for annotation
    Returns manhattan mappings
    """
    if new_annot:
        raw_rdp_file=rdp_file
        read_annotation_files(tree_parsed_dir,sim_file,map_file)
        manhat_file = list()
        for row in raw_rdp_file:
            rec_key=row[1].split(';')[0]
            map_bonds=[1,dom_mappings_dict[rec_key][2]] +dom_mappings_dict[rec_key][4:]
            overlaps=find_dom_overlaps(int(row[6]),int(row[7]), map_bonds)
            if max(overlaps[0:3])>0.35:
                if 'partial' not in overlaps[3]:
                    coords_lists=reformat_coords_for_manhattan(overlaps+[row[14]])
                    coords_lists[0]=coords_lists[0] + [row[0]]
                    coords_lists[1]=coords_lists[1] + [row[0]]
                    manhat_file.extend(coords_lists)
        write_csv(manhat_file,'recomb_manhattan_filtered_unpivot.tsv')
    else:
        raw_rdp_file = init_RDP_filter(rdp_file)
        read_annotation_files(tree_parsed_dir,sim_file,map_file)

        manhat_file = list()
        for row in raw_rdp_file:
            rec_key=row[0].split(';')[0]
      
            map_bonds=[1,dom_mappings_dict[rec_key][2]] +dom_mappings_dict[rec_key][4:]
            overlaps_filt=find_dom_overlaps(row[4],row[5], map_bonds, filt_mode=True)

            if 'other' in overlaps_filt[3]:
                over1=find_dom_overlaps(1,row[4], map_bonds, filt_mode=True)
                over3=find_dom_overlaps(row[5],dom_mappings_dict[rec_key][5], map_bonds, filt_mode=True)
                if over1[0]>0.6 and  over1[1]<0.2 and over1[2]<0.2:
                    row[5]=row[4]
                    row[4]=1
                elif over3[2]>0.6 and overlaps_filt[0]>0.6 and overlaps_filt[1]>0.6 and overlaps_filt[2]<0.2:
                    row[4]=row[5]
                    row[5]=dom_mappings_dict[rec_key][5]
                elif (overlaps_filt[1]>0.7 and overlaps_filt[2]>0.7 and overlaps_filt[0]<0.5):
                    row[5]=row[4]
                    row[4]=1
                elif over3[2]>0.3 and overlaps_filt[0]>0.7 and overlaps_filt[1]>0.7:
                    row[4]=row[5]
                    row[5]=dom_mappings_dict[rec_key][5]

            overlaps=find_dom_overlaps(row[4],row[5], map_bonds)
            if max(overlaps[0:3])>0.35:
                #if 'partial' not in overlaps[3]:
                manhat_file.extend(reformat_coords_for_manhattan(overlaps+[row[6]]))
        write_csv(manhat_file,'recomb_manhattan_with_partials.tsv')
  
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reads raw RDP table and cleans detected events')
    parser.add_argument('-o', '--out', dest='out_dir', help='the path to the output directory',
                        type=str)
    parser.add_argument('-r', '--rec', dest='rec_tab', help='the path to the raw RDP table',
                        type=str)
    parser.add_argument('-t', '--tree', dest='tree_path', help='the path to the nucleotide trees',
                        type=str)
    parser.add_argument('-p', '--parsed_t', dest='parsed_trees', help='the path to the parsed tree tables',
                        type=str)
    parser.add_argument('-s', '--sim_t', dest='sim_tab', help='the path to the table with similarities',
                        type=str)
    parser.add_argument('-m', '--map_f', dest='map_file', help='the path to the file with domain mappings',
                        type=str)
    args = parser.parse_args()

    RDP_filter_pipeline(read_csv_to_list(args.rec_tab, headless=True,delim=','), args.parsed_trees, args.sim_tab, args.map_file)
    make_manhattan_only(read_csv_to_list(args.rec_tab, headless=True,delim=','), args.parsed_trees, args.sim_tab, args.map_file, new_annot=False)


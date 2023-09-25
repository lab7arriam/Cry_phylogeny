import click
import os
import os.path
import pandas as pd
from collections import defaultdict
from overall_toxicity_summary import make_tox_summary_dict, make_strains_dictionary, make_asmbls_dictionary
from annotate_clusters_consistency import write_csv, read_csv_to_list
from get_breakpoints_from_alignments import  get_major_and_minor_parents, unpacking_list_of_lists
from count_host_changes import clean_tox_dict, crop_cry_name, get_parent_names


def make_stat_for_strains(strain_to_tox_dict):
    strain_res = [['Strain', 'Num', 'Rank']]
    new_strains_stat = [['Strain', 'Cry','Num']]

    for strain in strain_to_tox_dict:
       
        third_rank_tox_list = [crop_cry_name(cry) for cry in strain_to_tox_dict[strain]]
        strain_res.append([strain, len(strain_to_tox_dict[strain]),'Rank4'])
        strain_res.append([strain, len(set(third_rank_tox_list)), 'Rank3'])

        if strain == 'HS18-1':
            print(strain, len(strain_to_tox_dict[strain]),len(set(third_rank_tox_list)), list(set(third_rank_tox_list)))

        for cry in set(third_rank_tox_list):
            new_strains_stat.append([strain, cry, len(set(third_rank_tox_list))])

    write_csv(strain_res, 'Num_tox_for_strains_per_rank.csv')
    write_csv(new_strains_stat, 'Cry_strains_rank3.csv')



def make_stat_for_asmbls(asmbls_to_tox_dict):
    strain_res = [['Asmbl', 'Num', 'Rank']]

    for assembly in asmbls_to_tox_dict:
       
        third_rank_tox_list = [crop_cry_name(cry) for cry in asmbls_to_tox_dict[assembly]]
        strain_res.append([assembly, len(asmbls_to_tox_dict[assembly]),'Rank4'])
        strain_res.append([assembly, len(set(third_rank_tox_list)), 'Rank3'])

    write_csv(strain_res, 'Num_tox_for_assemblies_per_rank.csv')




def make_num_stat(res_to_tox_dict, mode):
    results_list = []
    for res in res_to_tox_dict:
        results_list.append([res, len(res_to_tox_dict[res]), mode])

    return(results_list)


def make_summary_of_cry_nums(res_list):
    final_res = [['Accession','Cry_num',  'Dataset']]
    for res_dict, mode in zip(res_list, ['Strains', 'Assemblies', 'Contigs', 'Plasmids', 'Chromosomes']):
        for inter_list in make_num_stat(res_dict, mode):
            final_res.append(inter_list)
    
    write_csv(final_res, 'summary_of_cry_nums_per_dataset.csv')



def get_all_toxins_from_by_dict(tox, tox_to_res_dict, res_to_tox_dict):
    tox_set = set()

    for res in tox_to_res_dict[tox]:
        for new_tox in res_to_tox_dict[res]:
            tox_set.add(new_tox)

    if tox!='unknown':
       tox_set.add(tox) 
	
    return(tox_set)



def get_jaccard_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/len(set(list1).union(set(list2))))

def get_simpson_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/min(len(set(list1)),len(set(list2))))

def get_braun_coeff(list1,list2):
    return(len(set(list1).intersection(set(list2)))/max(len(set(list1)),len(set(list2))))

def get_rec_coeff(recs,pars):
    return(len(set(recs).intersection(set(pars)))/len(recs))


def prepare_make_all_coeffs_calculations(rec_id, recs, minor_pars,major_inter_pars,major_all_pars, min_maj_strict_pars, min_maj_all_pars, dict_list, dataset, mode='for_coeffs'):
    rec_list = []
    par_strict_list = []
    par_all_list = []
    common_tox_res = []

    min_list = []
    maj_list_strict = []
    maj_list_all = []

    recs_total_tox = set()
    min_maj_strict_total_tox = set()
    min_maj_all_total_tox = set()

    min_total_tox = set()
    maj_strict_total_tox = set()
    maj_all_total_tox = set()

    for tox in minor_pars:
        crop_tox = crop_cry_name(tox)
        min_list.append(crop_tox)
        min_total_tox = min_total_tox.union(get_all_toxins_from_by_dict(crop_tox, dict_list[0], dict_list[1]))

    for tox in major_inter_pars:
        crop_tox = crop_cry_name(tox)
        maj_list_strict.append(crop_tox)
        maj_strict_total_tox = maj_strict_total_tox.union(get_all_toxins_from_by_dict(crop_tox, dict_list[0], dict_list[1]))

    for tox in major_all_pars:
        crop_tox = crop_cry_name(tox)
        maj_list_all.append(crop_tox)
        maj_all_total_tox = maj_all_total_tox.union(get_all_toxins_from_by_dict(crop_tox, dict_list[0], dict_list[1]))


    for tox in recs:
        crop_tox = crop_cry_name(tox)
        recs_total_tox = recs_total_tox.union(get_all_toxins_from_by_dict(crop_tox, dict_list[0], dict_list[1]))

        rec_list.append(crop_tox)

    for tox in min_maj_strict_pars:
        crop_tox = crop_cry_name(tox)

        min_maj_strict_total_tox = min_maj_strict_total_tox.union(get_all_toxins_from_by_dict(crop_tox, dict_list[0], dict_list[1]))

        par_strict_list.append(crop_tox)


    for tox in min_maj_all_pars:
        crop_tox = crop_cry_name(tox)

        min_maj_all_total_tox = min_maj_all_total_tox.union(get_all_toxins_from_by_dict(crop_tox, dict_list[0], dict_list[1]))
        par_all_list.append(crop_tox)

    if mode == 'for_combinations':
        return([rec_list, min_list, maj_list_strict, par_strict_list, maj_list_all, par_all_list])


    if len(min_maj_strict_total_tox)>0:
        jac = get_jaccard_coeff(recs_total_tox, min_maj_strict_total_tox)
        sim = get_simpson_coeff(recs_total_tox, min_maj_strict_total_tox)
        braun = get_braun_coeff(recs_total_tox, min_maj_strict_total_tox)
        rec_coeff = get_rec_coeff(recs_total_tox,min_maj_strict_total_tox)
        par_coeff = get_jaccard_coeff(min_total_tox, maj_strict_total_tox)

        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), jac, 'Jaccard', dataset, 'Strict'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), sim, 'Simpson', dataset, 'Strict'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), braun, 'Braun', dataset, 'Strict'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), rec_coeff, 'Rec_coeff', dataset, 'Strict'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), par_coeff, 'Par_coeff', dataset, 'Strict'])


    if len(min_maj_all_total_tox)>0:
        jac = get_jaccard_coeff(recs_total_tox, min_maj_all_total_tox)
        sim = get_simpson_coeff(recs_total_tox, min_maj_all_total_tox)
        braun = get_braun_coeff(recs_total_tox, min_maj_all_total_tox)
        rec_coeff = get_rec_coeff(recs_total_tox,min_maj_all_total_tox)
        par_coeff = get_jaccard_coeff(min_total_tox, maj_all_total_tox)

        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_all_total_tox), jac, 'Jaccard', dataset, 'All'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_all_total_tox), sim, 'Simpson', dataset, 'All'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), braun, 'Braun', dataset, 'All'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), rec_coeff, 'Rec_coeff', dataset, 'All'])
        common_tox_res.append([rec_id, len(recs_total_tox), len(min_maj_strict_total_tox), par_coeff, 'Par_coeff', dataset, 'All'])    

    return((rec_list, par_strict_list, par_all_list, common_tox_res))


def calculate_num_of_common_tox(rec_table, dict_list, dataset):

    file_reader = read_csv_to_list(os.path.realpath(rec_table), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    common_tox_res = []

    rec_list = []
    par_strict_list = []
    par_all_list = []

    for line in rec_df.iterrows():
       
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']
        
        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        min_maj_strict_pars = minor_pars+major_inter_pars
        min_maj_all_pars = minor_pars+major_all_pars

        res = prepare_make_all_coeffs_calculations(rec_id, recs, minor_pars,major_inter_pars,major_all_pars,min_maj_strict_pars, min_maj_all_pars, dict_list, dataset)
        rec_list.extend(res[0])
        par_strict_list.extend(res[1])
        par_all_list.extend(res[2])
        common_tox_res.extend(res[3])
        
    tox_sets = [list(set([el for el in rec_list if el!='unknown'])), list(set([el for el in par_strict_list if el!='unknown'])), list(set([el for el in par_all_list if el!='unknown']))]

    return((common_tox_res, tox_sets))


def get_tox_from_strain_to_compare(tox, tox_to_res_dict, res_to_tox_dict, tox_type, dataset, other_sets, mode = 'soft'):
    tox_num = []
    if tox_type == 'par_all':
        other_tox_sets = set(unpacking_list_of_lists(other_sets[:-2]))

    else:
        other_tox_sets = set(unpacking_list_of_lists(other_sets))

    for res in tox_to_res_dict[tox]:
        if mode =='strict':
            if len(set(res_to_tox_dict[res]).intersection(other_tox_sets))==0:
                tox_num.append((tox_type, res, len(res_to_tox_dict[res]), dataset, mode))
        else:
            tox_num.append((tox_type, res, len(res_to_tox_dict[res]), dataset, mode))


    return(tox_num)



def make_comparision_for_rec_and_all_tox(tox_to_strain_dict, strain_to_tox_dict, tox_to_asmbls_dict, asmbls_to_tox_dict, rec_tab, ref_list):

    tox_sets = calculate_num_of_common_tox(rec_tab, [tox_to_strain_dict, strain_to_tox_dict],'Strains')[1]
    ref_res = list(set([crop_cry_name(cry) for cry in unpacking_list_of_lists(read_csv_to_list(os.path.realpath(ref_list), headless=False))]))
    ref_no_rec = []
    
    for ref_tox in ref_res:
        if ref_tox not in unpacking_list_of_lists(tox_sets[:2]):
            ref_no_rec.append(ref_tox)

    
    
    num_tox_per_set = [['Group', 'Accession', 'Num', 'Dataset', 'Mode']]
    tox_sets = [ref_no_rec]+tox_sets[:3]

    for tox_set_ind, tox_type in zip([0,1,2,3, 4], ['ref_no_rec','rec', 'par_strict', 'par_all']): #ref_no_rec_strict, ref_no_rec_all

        other_sets = [tox_sets[ind] for ind in range(4) if ind!= tox_set_ind]

        num_pars_per_tox_set = []
        print(tox_type)
        
        for tox in tox_sets[tox_set_ind]:
            num_pars_per_tox_set.extend(get_tox_from_strain_to_compare(tox, tox_to_strain_dict, strain_to_tox_dict, tox_type, 'Strains', other_sets, 'strict'))
            num_pars_per_tox_set.extend(get_tox_from_strain_to_compare(tox, tox_to_asmbls_dict, asmbls_to_tox_dict, tox_type, 'Assemblies', other_sets, 'strict'))


            num_pars_per_tox_set.extend(get_tox_from_strain_to_compare(tox, tox_to_strain_dict, strain_to_tox_dict, tox_type, 'Strains', other_sets, 'soft'))
            num_pars_per_tox_set.extend(get_tox_from_strain_to_compare(tox, tox_to_asmbls_dict, asmbls_to_tox_dict, tox_type, 'Assemblies', other_sets, 'soft'))

        for tox_tuple in num_pars_per_tox_set:
            num_tox_per_set.append(list(tox_tuple))

    write_csv(num_tox_per_set, 'Tox_for_recs_and_refs.csv')

def check_tox_group(tox, tox_sets_dict):
    check = []
    for tox_set in tox_sets_dict:
        if tox in tox_sets_dict[tox_set]:
            check.append(tox_set)

    if len(check)==0:
        check.append('ref_no_rec')

    if len(check) >1 and 'rec' in check:
        check.remove('rec')
        check.append('both')

    return(check)

def classify_tox_from_dataset(tox_to_res_dict, res_to_tox_dict, rec_tab, ref_list, mode='Strains'):
       
    tox_sets = calculate_num_of_common_tox(rec_tab, [tox_to_res_dict, res_to_tox_dict],'Strains')[1]
    ref_res = list(set([crop_cry_name(cry) for cry in unpacking_list_of_lists(read_csv_to_list(os.path.realpath(ref_list), headless=False))]))
    ref_no_rec = []
    
    for ref_tox in ref_res:
        if ref_tox not in unpacking_list_of_lists(tox_sets):
            ref_no_rec.append(ref_tox)
    
    tox_sets = [ref_no_rec]+tox_sets
    
    classified_tox_types_res = []
    tox_sets_dict = {'ref_no_rec': tox_sets[0], 'rec': tox_sets[1], 'par_strict':tox_sets[2], 'par_all':tox_sets[3]}

    for res in res_to_tox_dict:
        for tox in res_to_tox_dict[res]:
            for group in check_tox_group(tox, tox_sets_dict):
                classified_tox_types_res.append([res, mode, group, len(res_to_tox_dict[res])])

    return(classified_tox_types_res)



def find_pars_and_recs_in_asmbls_strain(tox_list1, tox_list2, tox_to_res_dict, res_to_tox_dict, pair_type, dataset, rec_id, rec_type):
    tox_res1 = []
    tox_res2 = []

    for tox in tox_list1:
        if tox in tox_to_res_dict:
            tox_res1.extend(tox_to_res_dict[tox])

    for tox in tox_list2:
        if tox in tox_to_res_dict:
            tox_res2.extend(tox_to_res_dict[tox])

    tox_set_intersections = list(set(tox_res1).intersection(set(tox_res2)))

    res_list = []
    if len(tox_set_intersections) >0:
        for res in tox_set_intersections:
            for tox in res_to_tox_dict[res]:
                if tox in tox_list1:
                    res_list.append([rec_id, rec_type, 'rec',tox,res,dataset])
                if tox in tox_list2:
                    res_list.append([rec_id, rec_type, pair_type,tox,res,dataset])

    return(res_list)


def check_presence_of_parents_and_recs_in_events(rec_table, dict_list, dataset):

    file_reader = read_csv_to_list(os.path.realpath(rec_table), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])


    combinations_results = []

    for line in rec_df.iterrows():
       
        rec_id = line[1]['ID']
        recs = get_parent_names(line[1]['Rec'])
        rec_type = line[1]['Type']
        
        minor_pars, major_all_pars, major_inter_pars = get_major_and_minor_parents(line, mode='intersect')
        min_maj_strict_pars = minor_pars+major_inter_pars
        min_maj_all_pars = minor_pars+major_all_pars

        res = prepare_make_all_coeffs_calculations(rec_id, recs, minor_pars,major_inter_pars,major_all_pars,min_maj_strict_pars, min_maj_all_pars, dict_list, dataset, 'for_combinations')
        
        extra_pars = set(res[4]) - set(res[2])
        combinations_results.extend(find_pars_and_recs_in_asmbls_strain(res[0], res[1], dict_list[0],dict_list[1], 'min', dataset, rec_id, rec_type))
        combinations_results.extend(find_pars_and_recs_in_asmbls_strain(res[0], res[2], dict_list[0],dict_list[1], 'maj', dataset, rec_id, rec_type))
        combinations_results.extend(find_pars_and_recs_in_asmbls_strain(res[0], list(extra_pars), dict_list[0],dict_list[1], 'extra_pars', dataset, rec_id, rec_type))


        #4 Rank
        #extra_pars = set(min_maj_all_pars) - set(min_maj_strict_pars)
        #combinations_results.extend(find_pars_and_recs_in_asmbls_strain(recs, minor_pars, dict_list[0],dict_list[1], 'min', dataset, rec_id, rec_type))
        #combinations_results.extend(find_pars_and_recs_in_asmbls_strain(recs, major_inter_pars, dict_list[0],dict_list[1], 'maj', dataset, rec_id, rec_type))
        #combinations_results.extend(find_pars_and_recs_in_asmbls_strain(recs, list(extra_pars), dict_list[0],dict_list[1], 'extra_pars', dataset, rec_id, rec_type))


    return(combinations_results)



def mearged_combinations_stat(asmbl_stat, plasm_stat, strains_stat, rec_tab, tox_summary_dict, ref_list):
    #tox_to_strain_dict, strain_to_tox_dict =  reduce_dicts_to_third_rank(*make_strains_dictionary(strains_stat))

    tox_to_strain_dict, strain_to_tox_dict =  make_strains_dictionary(strains_stat)


    tox_to_asmbls_dict, asmbls_to_tox_dict =  make_asmbls_dictionary(asmbl_stat, plasm_stat, mode = 'asmbl')
    tox_to_cont_dict, cont_to_tox_dict =  make_asmbls_dictionary(asmbl_stat, plasm_stat, mode = 'contig')
    tox_to_plasm_dict, plasm_to_tox_dict =  make_asmbls_dictionary(asmbl_stat, plasm_stat, mode = 'plasm')
    tox_to_chrom_dict, chrom_to_tox_dict =  make_asmbls_dictionary(asmbl_stat, plasm_stat, mode = 'chrom')
 
    #make_summary_of_cry_nums([strain_to_tox_dict, asmbls_to_tox_dict, cont_to_tox_dict, plasm_to_tox_dict, chrom_to_tox_dict])
    print(len(tox_to_plasm_dict))
    print(len(tox_to_chrom_dict))

    common_tox_res = [['ID', 'Rec_num', 'Par_num', 'Comon_num', 'Coeff', 'Dataset', 'Par_type']]
    common_tox_res.extend(calculate_num_of_common_tox(rec_tab, [tox_to_strain_dict, strain_to_tox_dict],'Strains')[0])
    common_tox_res.extend(calculate_num_of_common_tox(rec_tab, [tox_to_asmbls_dict, asmbls_to_tox_dict],'Assemblies')[0])
    
    write_csv(common_tox_res, 'Num_common_tox_pars_recs.csv')

    make_comparision_for_rec_and_all_tox(tox_to_strain_dict, strain_to_tox_dict, tox_to_asmbls_dict, asmbls_to_tox_dict, rec_tab, ref_list)


    classified_tox_types_res = [['Accession', 'Dataset', 'Tox_type', 'Total_tox_num']]
    combinations_parents = [['ID', 'Domain', 'Tox_type', 'Tox', 'Common_accession','Dataset_type']]
    
    for tox_set_pair, dataset_type in zip([(tox_to_strain_dict, strain_to_tox_dict),(tox_to_asmbls_dict, asmbls_to_tox_dict ),(tox_to_plasm_dict, plasm_to_tox_dict),(tox_to_chrom_dict, chrom_to_tox_dict)], ['Strains', 'Assemblies', 'Plasmids', 'Chromosomes']):
        classified_tox_types_res.extend(classify_tox_from_dataset(tox_set_pair[0], tox_set_pair[1], rec_tab, ref_list, dataset_type))
        combinations_parents.extend(check_presence_of_parents_and_recs_in_events(rec_tab, [tox_set_pair[0], tox_set_pair[1]],dataset_type))


    combinations_parents.extend(check_presence_of_parents_and_recs_in_events(rec_tab, [tox_to_cont_dict, cont_to_tox_dict],'Contig'))


    write_csv(classified_tox_types_res, 'Num_tox_types_per_dataset_with_both.csv')
    write_csv(combinations_parents, 'Recs_and_parents_combinations_Rank3.csv')
    write_csv(combinations_parents, 'Recs_and_parents_combinations_Rank4.csv')


@click.command()           
@click.option('--rec_tab', '-r', help="the path with recombination results", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--spec_tab', '-s', help="the path to the table denoting host specificity", 
              type=click.Path(exists=True), metavar='<PATH>')  
@click.option('--strains_stat', '-t', help="the path to the table denoting strains attributions", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--ref_list', '-f', help="the path to the list of reference proteins", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--plasm_stat', '-p', help="the path to attributions of nucleotide acessions in the assemblies", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--asmbl_stat', '-a', help="the path to cry proteins disctribution on assemblies", 
              type=click.Path(exists=True), metavar='<PATH>')



def main(rec_tab, spec_tab, strains_stat, ref_list, plasm_stat, asmbl_stat): 

    tox_df = pd.read_csv(os.path.realpath(spec_tab), sep='\t')
    tox_summary_dict = clean_tox_dict(make_tox_summary_dict(tox_df))

    mearged_combinations_stat(asmbl_stat, plasm_stat, strains_stat, rec_tab, tox_summary_dict, ref_list)

if __name__ == '__main__':
   main()





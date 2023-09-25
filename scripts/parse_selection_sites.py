
import pandas as pd
import os
import os.path
import click
import random
from Bio import SeqIO
from collections import defaultdict
from annotate_clusters_consistency import write_csv, read_csv_to_list


def get_best_model(row, pvals_inds, ln_inds):
    if row[ln_inds[0]]=='-':
       return(None)

    ln_dict = {ind:float(row[ind]) for ind in ln_inds}
    pval_dict = {ind:float(row[ind]) for ind in pvals_inds}
    susp_model = [key for key in ln_dict if ln_dict[key]==max(list(ln_dict.values()))][0]
    pval_tests = []
    evol_type = ''
    
    for pval_key in pval_dict:
        model_ind = susp_model.split('_')[1].lower().replace('bsa1','mA1').replace('bsa','mA').replace('bsb','mB')
        if model_ind in pval_key:
            pass
            #if 'M1' in susp_model: 
            #    print(susp_model,ln_dict[susp_model], pval_key, pval_dict[pval_key])


    model_ind = susp_model.split('_')[1].lower().replace('bsa1','mA1').replace('bsa','mA').replace('bsb','mB')

    if 'M1' in susp_model: 
        if pval_dict[ln_inds[0].split('_')[0]+'_'+'m0_vs_m1'] < 0.05: 
            evol_type='relaxation'

    if 'M8' in susp_model: 
        if pval_dict[ln_inds[0].split('_')[0]+'_'+'m0_vs_m8'] < 0.05:
            if pval_dict[ln_inds[0].split('_')[0]+'_'+'m7_vs_m8'] < 0.05:
                evol_type='positive'
            else:
                evol_type='non-signif'
 
    if 'M2' in susp_model: 
        if pval_dict[ln_inds[0].split('_')[0]+'_'+'m0_vs_m2'] < 0.05:
            if pval_dict[ln_inds[0].split('_')[0]+'_'+'m1_vs_m2'] < 0.05:
                evol_type='positive'
            else:
                evol_type='non-signif'
        else:
            evol_type='non-signif' 

    if 'M3' in susp_model: 
        if pval_dict[ln_inds[0].split('_')[0]+'_'+'m3_vs_m0'] < 0.05:
            if pval_dict[ln_inds[0].split('_')[0]+'_'+'m7_vs_m8']  < 0.05:
                evol_type='positive'
                susp_model='rec_M8_lnL'
            else:
                evol_type='relaxation'

    if 'mB' in model_ind: 
        if pval_dict[ln_inds[0].split('_')[0]+'_'+'mA1_vs_mA'] < 0.05:
            evol_type='positive'

        elif pval_dict[ln_inds[0].split('_')[0]+'_'+'mA1_vs_mB'] < 0.05:
            evol_type='positive'
        else:
            evol_type='non-signif'

    if 'mA' in model_ind: 
        if pval_dict[ln_inds[0].split('_')[0]+'_'+'mA1_vs_mA'] < 0.05:
            evol_type='positive'

        elif pval_dict[ln_inds[0].split('_')[0]+'_'+'mA1_vs_mB'] < 0.05:
            evol_type='positive'
            susp_model='rec_bsB_lnL'
        else:
            evol_type='non-signif'

    return([susp_model.split('_')[1], evol_type])


def check_branch_model(row, pref):
    if row[1][pref+'omega_b']=='-':
        return(None)

    background = float(row[1][pref+'omega_b'])
    foreground = float(row[1][pref+'omega_1'])

    m0_check = float(row[1][pref+'M0_vs_bfree'])
    bfree_check = float(row[1][pref+'bneut_vs_bfree'])
   
    evol_type, ground_check = '', ''

    if bfree_check > 0.05:
        evol_type = 'non-signif'
    else:
        if foreground > 1:
            evol_type = 'positive'
        else:
            evol_type = 'purifying'

    if m0_check > 0.05:
        ground_check = 'non-signif'
    else:
        ground_check = 'signif'

    return([evol_type, ground_check, background, foreground])

def get_coords_percent(model_col, fasta_link, domain, site_type, tox_type):
    coords_dict={}	
    return_dict=defaultdict(list)
    for item in model_col.split(';'):
        if '-' not in item:
            coords_dict[int(item.split('(')[0])]=item.split('(')[1].replace(')','')
        else:
            if item!='-':
                prob = item.split('(')[1].replace(')','')
                items = [x for x in range(int(item.split('(')[0].split('-')[0]),int(item.split('(')[0].split('-')[1]))]
                
                for i in items:
                   coords_dict[i] = prob

    for record in SeqIO.parse(fasta_link, "fasta"):
        fasta_length = len(record.seq)
        break

    fasta_length /= 3

    if len(coords_dict)  < fasta_length * 0.4:
        for i in coords_dict:
            coord_percent = i
            coord_percent /= fasta_length
            coord_percent *= 100

            if domain=='domain2':
                coord_percent += 100
            elif domain=='domain3':
                coord_percent += 200

            return_dict[coord_percent] = [coords_dict[i], 'passed', site_type, tox_type]

    elif len(coords_dict)  >= fasta_length * 0.4:
        for i in coords_dict:
            coord_percent = i
            coord_percent /= fasta_length
            coord_percent *= 100

            if domain=='domain2':
                coord_percent += 100
            elif domain=='domain3':
                coord_percent += 200

            return_dict[coord_percent] = [coords_dict[i], 'not_passed', site_type, tox_type]

    return(return_dict)


def parse_sites_data(evol_df, evol_dir, simp_tab, ord_tab):
    branch_inds = ['M0_vs_bfree', 'bneut_vs_bfree', 'M0_lnL', 'b_bfree_L', 'b_neut_L']
    branch_site_inds = ['m0_vs_mA', 'm0_vs_mA1', 'm0_vs_mB', 'mA1_vs_mA', 'mA_vs_mB', 'mA1_vs_mB','M0_lnL', 'bsA_lnL', 'bsA1_lnL', 'bsB_lnL']
    site_inds = ['m0_vs_m1','m0_vs_m2', 'm3_vs_m0','m0_vs_m7', 'm0_vs_m8','m1_vs_m2','m1_vs_m3','m1_vs_m8','m2_vs_m3','m7_vs_m2', 'm7_vs_m3','m7_vs_m8', 'm8_vs_m3',
                  'M0_lnL', 'M1_lnL', 'M2_lnL', 'M3_lnL', 'M7_lnL', 'M8_lnL']

    branch_results = [['ID','Rec_domain', 'Curr_domain','Tox_type','evol_type', 'signif_flag',  'omega_b','omega_f','Num_tox','Num_rec_pars']]

    #branch_results_extended = [['ID','Rec_domain', 'Curr_domain','Tox_type', 'unknown_flag','evol_type', 'signif_flag',  'omega_b','omega_f','Num_tox','Num_rec_pars', 'Num_recs_pars_all', 
    #                   'Sites_pos', 'Sites_cons', 'Sites_rel','Best_branc_site', 'Brsites_pos', 'Brsites_cons', 'Brsites_rel','Passed_flag', 'Orders_flag','Simpson_species']]

    branch_results_extended = [['ID','Rec_domain', 'Curr_domain','Tox_type', 'unknown_flag','evol_type', 'signif_flag',  'omega_b','omega_f','Num_tox','Num_rec_pars', 'Num_recs_pars_all', 'best_site_model','Sites_pos', 'Sites_cons', 'Sites_rel','Best_branc_site', 'Brsites_pos', 'Brsites_cons', 'Brsites_rel','Passed_flag', 'Orders_flag','Simpson_species']]

    site_results = [['ID','Rec_domain', 'Curr_domain','Tox_type', 'Model', 'signif_flag','Coord', 'Site_prob','Site_flag','Site_type']]

    orders_dict = dict()
    simp_dict = dict()

    orders_table = read_csv_to_list(os.path.realpath(ord_tab), headless=True)
    for row in orders_table:
        if row[5]=='Strains_add':
            orders_dict[row[0]] = row[1]

    simp_table = read_csv_to_list(os.path.realpath(simp_tab), headless=True)
    for row in simp_table:
        if row[5]=='Strains_add':
            if row[2]=='all':
                simp_dict[row[0]] = row[6]
 
    for row in evol_df.iterrows():
        ID=row[1]['Id']
        dom_type=row[1]['Rec_flag'].replace('Rec','Min').replace('Maj1','Maj').replace('Maj2','Maj')
        curr_domain=row[1]['Current_domain'] 
        domain_type=row[1]['Type']

        num_recs = len(row[1]['Rec'].split(';'))
        num_tox = row[1]['len']
        unknown_par_flag = 'known'

        if str(ID) in orders_dict:
            diff_flag = orders_dict[str(ID)]
        else:
            diff_flag = 'NA'

        if str(ID) in simp_dict:
            simp_res = simp_dict[str(ID)]
        else:
            simp_res = 'NA'

        if row[1]['Rec_flag']=='Rec':
            event_type='minor'
            num_pars = len(row[1]['Min_par'].split(';'))
            if row[1]['Min_par']=='unknown':
                num_pars= 0
                unknown_par_flag = 'minor'
                
        if row[1]['Rec_flag']=='Maj1':
            event_type='major1'
            num_pars = len(row[1]['Maj_par1'].split(';'))
            if row[1]['Maj_par1']=='unknown':
                num_pars= 0
                unknown_par_flag = 'major'
                
        if row[1]['Rec_flag']=='Maj2':
            event_type='major2'
            num_pars = len(row[1]['Maj_par2'].split(';'))
            if row[1]['Maj_par2']=='unknown':
                num_pars= 0
                unknown_par_flag = 'major'

        num_recs_pars_all = num_pars+num_recs

        event_dir=os.path.join(os.path.realpath(evol_dir),event_type+'_event_'+str(ID),event_type+'_event_'+str(ID)+'_nucl_reported.fasta')

        branch_site_par = get_best_model(row[1], ['par_' +el for el in branch_site_inds[:6]], ['par_' +el for el in branch_site_inds[6:]])
        branch_site_rec = get_best_model(row[1], ['rec_' +el for el in branch_site_inds[:6]], ['rec_' +el for el in branch_site_inds[6:]])
        best_site_model = get_best_model(row[1], ['rec_' +el for el in site_inds[:13]], ['rec_' +el for el in site_inds[13:]])

        branch_res_par = check_branch_model(row, 'par_')
        branch_res_rec = check_branch_model(row, 'rec_')
        
        site_model = best_site_model[0].lower()
        
        c_sites = get_coords_percent(row[1]['rec_'+site_model+'_sites_c'], event_dir, curr_domain, 'cons','rec')
        s_sites = get_coords_percent(row[1]['rec_'+site_model+'_sites_s'], event_dir, curr_domain, 'pos','rec')
        r_sites = get_coords_percent(row[1]['rec_'+site_model+'_sites_r'], event_dir, curr_domain, 'relaxed','rec')

        for site_set in [c_sites, s_sites,r_sites]:
            if len(site_set)>0:
                for site_coord in site_set:
                    site_results.append([ID, domain_type,curr_domain, 'Rec', site_model, best_site_model[1], site_coord] + site_set[site_coord][:-1])


        if branch_site_par:
            bs_par_model =  branch_site_par[0].replace('bs', 'm')

            c_brsites_par = get_coords_percent(row[1]['par_'+bs_par_model+'_sites_c'], event_dir, curr_domain, 'cons','par')
            s_brsites_par = get_coords_percent(row[1]['par_'+bs_par_model+'_sites_s'], event_dir, curr_domain, 'pos','par')
            r_brsites_par = get_coords_percent(row[1]['par_'+bs_par_model+'_sites_r'], event_dir, curr_domain, 'relaxed','par')


            for site_set in [c_brsites_par, s_brsites_par,r_brsites_par]:
                if len(site_set)>0:
                    for site_coord in site_set:
                        site_results.append([ID, domain_type,curr_domain, 'Par', bs_par_model, branch_site_par[1], site_coord] + site_set[site_coord][:-1])


        if branch_site_rec:
            bs_rec_model =  branch_site_rec[0].replace('bs', 'm')

            c_brsites_rec = get_coords_percent(row[1]['rec_'+bs_rec_model+'_sites_c'], event_dir, curr_domain, 'cons','rec')
            s_brsites_rec = get_coords_percent(row[1]['rec_'+bs_rec_model+'_sites_s'], event_dir, curr_domain, 'pos','rec')
            r_brsites_rec = get_coords_percent(row[1]['rec_'+bs_rec_model+'_sites_r'], event_dir, curr_domain, 'relaxed','rec')


            for site_set in [c_brsites_rec, s_brsites_rec,r_brsites_rec]:
                if len(site_set)>0:
                    for site_coord in site_set:
                        site_results.append([ID, domain_type,curr_domain, 'Rec', bs_rec_model, branch_site_rec[1], site_coord] + site_set[site_coord][:-1])

        if branch_res_rec:
            branch_results.append([ID, domain_type,curr_domain, 'Rec']+ branch_res_rec +[num_tox,num_recs])
            rec_row = [ID, domain_type,curr_domain, 'Rec', unknown_par_flag]+ branch_res_rec +[num_tox,num_recs, num_recs_pars_all]
            sites_row = [best_site_model[1],site_model, len(s_sites), len(c_sites), len(r_sites)]
            branch_site_row = ['NA','NA','NA','NA']

            if branch_site_rec:

                passed_flag = 'passed'
                for site_set in [c_brsites_rec, s_brsites_rec,r_brsites_rec]:
                    if len(site_set)>0:
   
                        if random.choice(list(site_set.values()))[1]=='not_passed':
                            passed_flag='not_passed'


                branch_site_row = [bs_rec_model, len(s_brsites_rec), len(c_brsites_rec), len(r_brsites_rec), passed_flag]

            branch_results_extended.append(rec_row+sites_row+branch_site_row+[diff_flag, simp_res])



        if branch_res_par:
            branch_results.append([ID, domain_type,curr_domain, 'Par']+ branch_res_par + [num_tox,num_pars])

            par_row = [ID, domain_type,curr_domain, 'Par', unknown_par_flag]+ branch_res_par +[num_tox,num_pars, num_recs_pars_all]
            sites_row = [best_site_model[1], site_model, len(s_sites), len(c_sites), len(r_sites)]

            branch_site_row = ['NA','NA','NA','NA']

            if branch_site_par:
                passed_flag = 'passed'
                for site_set in [c_brsites_par, s_brsites_par,r_brsites_par]:
                    if len(site_set)>0:
   
                        if random.choice(list(site_set.values()))[1]=='not_passed':
                            passed_flag='not_passed'

                branch_site_row = [bs_par_model, len(s_brsites_par), len(c_brsites_par), len(r_brsites_par), passed_flag]

            branch_results_extended.append(par_row+sites_row+branch_site_row+[diff_flag, simp_res])
  
    write_csv(branch_results,'Branch_models_classified_res.csv')
    write_csv(branch_results_extended,'Branch_models_with_site_models.csv')
    write_csv(site_results,'Site_models_classified_res.csv')


@click.command()           
@click.option('--evol_dir', '-e', help="the path with evolutionary_selection results", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--summ_tab', '-s', help="the path to file with selection results", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--simp_tab', '-t', help="the path to file with simpson coefficient based on affected species", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--ord_tab', '-o', help="the path to file with orders change flag", 
              type=click.Path(exists=True), metavar='<PATH>')


def main(evol_dir, summ_tab, simp_tab, ord_tab): 

    evol_df = pd.read_csv(os.path.realpath(summ_tab), sep='\t')
    parse_sites_data(evol_df, evol_dir, simp_tab, ord_tab)

if __name__ == '__main__':
   main()

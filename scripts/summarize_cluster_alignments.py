import sys
from Bio import SeqIO, Entrez
import argparse
import os.path
import os
from agregate_selection_results import read_site_model, read_branch_model, read_branch_site_model
from annotate_clusters_consistency import write_csv, read_csv_to_list, create_dict_from_list
from collections import defaultdict, Counter
import pandas as pd
from itertools import combinations
from make_data_for_heatmap_universal import find_sequence_identity


def get_coords_percent(model_col, fasta_link, domain):
    coords_list=[]
    for item in model_col.split(';'):
        if '-' not in item:
            coords_list.append(int(item))
        else:
            if item!='-':
                items = [int(x) for x in range(int(item.split('-')[0]),int(item.split('-')[1]))]
                for i in items:
                   coords_list.append(int(i))

    for record in SeqIO.parse(fasta_link, "fasta"):
        fasta_length = len(record.seq)
        break

    fasta_length /= 3

    if len(coords_list)  < fasta_length * 0.4:
        for i in range(len(coords_list)):
            coords_list[i] /= fasta_length
            coords_list[i] *= 100
            if domain=='domain2':
                coords_list[i] += 100
            elif domain=='domain3':
                coords_list[i] += 200

        return(coords_list)

    else:
        return([])


def get_coords_percent_full(model_col, fasta_link, coord_dict, cluster):
    coords_list=[]
    for item in model_col.split(';'):
        if '-' not in item:
            coords_list.append(int(item))
        else:
            if item!='-':
                items = [int(x) for x in range(int(item.split('-')[0]),int(item.split('-')[1]))]
                for i in items:
                   coords_list.append(int(i))

    for record in SeqIO.parse(fasta_link, "fasta"):
       fasta_length = len(record.seq) // 3
       break

    if len(coords_list) * 2.5 < fasta_length:
        for domain in ['domain1', 'domain2','domain3']:
            coords_domain=coord_dict[cluster][int(domain[-1])-1]

            for i in range(len(coords_list)):
                if domain=='domain1' and coords_list[i]<=coords_domain[1] and coords_list[i]>=coords_domain[0]:

                    coords_list[i] /= coords_domain[1]
                    coords_list[i] *= 100

                elif domain=='domain2' and coords_list[i]<=coords_domain[1] and coords_list[i]>=coords_domain[0]:
                    coords_list[i] /= coords_domain[1]
                    coords_list[i] *= 100
                    coords_list[i] += 100

                elif domain=='domain3' and coords_list[i]<=coords_domain[1] and coords_list[i]>=coords_domain[0]:
                    coords_list[i] /= coords_domain[1]
                    coords_list[i] *= 100
                    coords_list[i] += 200
        return(coords_list)

    else:
        return([])

def ananlyze_selection(evol_dir: str) -> list:
    """
        Summarizes sites and branch models for multi-clusters
    """
    results=[]
    for seq_dir in ['domain1', 'domain2','domain3','full']:
        type_dir=os.path.join(os.path.realpath(evol_dir),seq_dir)
        for clust in os.listdir(type_dir):
            clust_dir=os.path.join(os.path.realpath(type_dir),clust)
            sites=os.path.join(os.path.realpath(clust_dir),'site_model.log')
            branches=os.path.join(os.path.realpath(clust_dir),'branch_model.log')
            branch_sites=os.path.join(os.path.realpath(clust_dir),'branch_site_model.log')

            site_res=read_site_model(sites)
            branch_res=read_branch_model(branches)
            branch_site_res=read_branch_site_model(branch_sites)
            results.append([seq_dir, clust]+branch_res+site_res+branch_site_res)
    return(results)

def dislose_sites(site_list, site_name, parse_row, domain, dom_type, cluster):
    ret_list=[]
    if len(site_list)>0:
        if 'full' not in dom_type: 
            mark_type='single'
        else:
            mark_type='full'

        if 'm8' in site_name and 'par' not in site_name:
            pval=parse_row[1]['m7_vs_m8']
        elif 'm2' in site_name and 'par' not in site_name:
            pval=parse_row[1]['m1_vs_m2']

        elif 'mA' in site_name and 'par' not in site_name:
            pval=parse_row[1]['mA1_vs_mA']

        elif 'mB' in site_name and 'par' not in site_name:
            pval=parse_row[1]['mA1_vs_mB']


        if site_name.split('_')[2]=='s': 
            site_type='pos'
        else:
            site_type='cons'

        for site in site_list:
            ret_list.append([cluster, domain, site, pval, site_name.split('_')[0], mark_type, site_type])

        return(ret_list)
    else:
        return([])

def aggregate_sites(evol_dir, sum_evol):
    write_list=[]
    write_list.append(['ID','domain','coord','pval', 'model','mark_type', 'site_type'])

    summary_csv=read_csv_to_list(sum_evol, headless=False)
    summary_df=pd.DataFrame(summary_csv[1:], columns=summary_csv[0])
    
    coord_dict=defaultdict(list)
    for row in summary_df.iterrows():

        if row[1]['Type']!='full':

            dom_type='single'
            domain=row[1]['Type']
            dom_dir=os.path.join(os.path.realpath(evol_dir),domain)
            clust_dir=os.path.join(os.path.realpath(dom_dir),row[1]['Cluster'])
            fasta_file=os.path.join(os.path.realpath(clust_dir),row[1]['Cluster']+'_extracted_reported.fasta')

            for record in SeqIO.parse(fasta_file, "fasta"):
               coord_len = len(record.seq)//3
               break

            if row[1]['Cluster'] not in coord_dict:
                coord_dict[row[1]['Cluster']]=[[0, coord_len]]
            else:
                coord_dict[row[1]['Cluster']].append([0, coord_len])

            for model in ['m8_sites_s', 'm8_sites_c', 'm2_sites_s','m2_sites_c','mA_sites_s','mA_sites_c','mB_sites_s','mB_sites_c']:

                model_sites=get_coords_percent(row[1][model], fasta_file, domain)
                model_res=dislose_sites(model_sites, model, row, domain, dom_type, row[1]['Cluster'])
                
                if len(model_res)>0:
                    for res in model_res:
                        write_list.append(res)

        else:
            dom_type='full'
            domain=row[1]['Type']

            for ind in range(1,3):
                coord_dict[row[1]['Cluster']][ind][0]+=coord_dict[row[1]['Cluster']][ind-1][1]+1
                coord_dict[row[1]['Cluster']][ind][1]+=coord_dict[row[1]['Cluster']][ind][0]

                dom_dir=os.path.join(os.path.realpath(evol_dir),domain)

            clust_dir=os.path.join(os.path.realpath(dom_dir),row[1]['Cluster'])
            fasta_file=os.path.join(os.path.realpath(clust_dir),row[1]['Cluster']+'_extracted_reported.fasta')

            for model in ['m8_sites_s', 'm8_sites_c', 'm2_sites_s','m2_sites_c','mA_sites_s','mA_sites_c','mB_sites_s','mB_sites_c']:
                full_coords = get_coords_percent_full(row[1][model], fasta_file, coord_dict, row[1]['Cluster'])
                if len(full_coords)>0:
                    for coord in full_coords:
                        if coord>=0 and coord<=100:
                            model_res=dislose_sites([coord], model, row, 'domain1', dom_type, row[1]['Cluster'])
                            write_list.append(model_res[0])

                        if coord>100 and coord<=200:
                            model_res=dislose_sites([coord], model, row, 'domain2', dom_type, row[1]['Cluster'])
                            write_list.append(model_res[0])

                        if coord>200 and coord<=300:
                            model_res=dislose_sites([coord], model, row, 'domain3', dom_type, row[1]['Cluster'])
                            write_list.append(model_res[0])

    write_csv(write_list,'multiclusters_selection_sites_summary.csv')


def make_id_distr_for_multi(alignment_dict_nucs):

    ids_distr = [['cluster', 'domain', 'id']]
    for domain in alignment_dict_nucs:
        for cluster in alignment_dict_nucs[domain]:
            print(cluster, domain)
            seqs_list = []
            for seq in alignment_dict_nucs[domain][cluster]:
                seqs_list.append(seq.replace('-',''))
            for seq_pair in list(combinations(seqs_list,2)):
                ident = find_sequence_identity(seq_pair, ret_mode='ident')[0]
                ids_distr.append([cluster,domain, ident])
    write_csv(ids_distr,'multiclusters_id_per_domain.csv')


def expand_annotated_table(align_dir: str) -> list:
    """
        Aggregates sequences witthin multiclusters
    """
    alignment_dict_prots=defaultdict(dict)
    alignment_dict_nucs=defaultdict(dict)


    for domain in ['domain1','domain2','domain3']:

        alignment_dict_prots[domain]=dict()
        alignment_dict_nucs[domain]=dict()
        domain_dir=os.path.join(os.path.realpath(align_dir),domain)
        msa_prot_dir=os.path.join(os.path.realpath(domain_dir),'msa_prot')
        msa_nucl_dir=os.path.join(os.path.realpath(domain_dir),'msa_nucl_reported')

        for cluster_file in os.listdir(msa_prot_dir):

            cluster=cluster_file.replace('_extracted_prot.msa','')
            records=[]
            aligment_file=os.path.join(os.path.realpath(msa_prot_dir),cluster_file)

            for record in SeqIO.parse(aligment_file, "fasta"):
                records.append(str(record.seq))

            if domain=='domain1':
                alignment_dict_prots['domain1'][cluster]=records
            if domain=='domain2':
                alignment_dict_prots['domain2'][cluster]=records
            if domain=='domain3':
                alignment_dict_prots['domain3'][cluster]=records

        for cluster_file in os.listdir(msa_nucl_dir):

            cluster=cluster_file.replace('_extracted_reported.fasta','')
            records=[]
            aligment_file=os.path.join(os.path.realpath(msa_nucl_dir),cluster_file)

            for record in SeqIO.parse(aligment_file, "fasta"):
                records.append(str(record.seq))

            if domain=='domain1':
                alignment_dict_nucs['domain1'][cluster]=records
            if domain=='domain2':
                alignment_dict_nucs['domain2'][cluster]=records
            if domain=='domain3':
                alignment_dict_nucs['domain3'][cluster]=records


    mis_dict_prots=defaultdict(dict)
    mis_dict_nucs=defaultdict(dict)

    for domain in alignment_dict_nucs:
        mis_dict_nucs[domain]=dict()
        for clust in alignment_dict_nucs[domain]:
            mis_dict_nucs[domain][clust]=dict()
            for seq in alignment_dict_nucs[domain][clust]:
                for ind in range(0, len(seq)-2,3):
                    if ind//3 not in mis_dict_nucs[domain][clust]:
                        mis_dict_nucs[domain][clust][ind//3]=[seq[ind:ind+3]]
                    else:
                        mis_dict_nucs[domain][clust][ind//3].append(seq[ind:ind+3])
            for ind in mis_dict_nucs[domain][clust]:
                mis_dict_nucs[domain][clust][ind]=dict(Counter(mis_dict_nucs[domain][clust][ind]))


    for domain in alignment_dict_prots:
        mis_dict_prots[domain]=dict()
        for clust in alignment_dict_prots[domain]:
            mis_dict_prots[domain][clust]=dict()
            for seq in alignment_dict_prots[domain][clust]:
                for ind in range(len(seq)):
                    if ind not in mis_dict_prots[domain][clust]:
                        mis_dict_prots[domain][clust][ind]=[seq[ind]]
                    else:
                        mis_dict_prots[domain][clust][ind].append(seq[ind])

            for ind in mis_dict_prots[domain][clust]:
                mis_dict_prots[domain][clust][ind]=dict(Counter(mis_dict_prots[domain][clust][ind]))


    mismatch_ready_dict=dict()
    passed_inds=dict()

    for domain in mis_dict_prots:
        mismatch_ready_dict[domain]=dict()
        passed_inds[domain]=dict()

        for clust in mis_dict_prots[domain]:
            mismatch_ready_dict[domain][clust]=[]
            passed_inds[domain][clust]=list()
            align_len=len(mis_dict_prots[domain][clust])

            for ind in mis_dict_prots[domain][clust]:
                if len(mis_dict_prots[domain][clust][ind])==2 and '-' in mis_dict_prots[domain][clust][ind]:
                    if ind>=3 and align_len-ind>=4:
                        mismatch_ready_dict[domain][clust].append([align_len, ['indel', ind]])
                        print([align_len, ['indel', ind]], domain, clust)
                        passed_inds[domain][clust].append(ind)

                elif len(mis_dict_prots[domain][clust][ind])>=2 and '-' not in mis_dict_prots[domain][clust][ind]:
                    mismatch_ready_dict[domain][clust].append([align_len, ['non_syn', ind]])
                    passed_inds[domain][clust].append(ind)

                elif len(mis_dict_prots[domain][clust][ind])>=2 and '-'  in mis_dict_prots[domain][clust][ind]:
                    if ind>=3 and align_len-ind>=4:
                        mismatch_ready_dict[domain][clust].append([align_len, ['non_syn|indel', ind]])
                        print([align_len, ['non_syn|indel', ind]], domain, clust)
                        passed_inds[domain][clust].append(ind)

    for domain in mis_dict_nucs:
        for clust in mis_dict_nucs[domain]:
            align_len=len(mis_dict_nucs[domain][clust])
            for ind in mis_dict_nucs[domain][clust]:
                if len(mis_dict_nucs[domain][clust][ind])>=2 and '---' not in mis_dict_nucs[domain][clust][ind] and ind not in passed_inds[domain][clust]:
                    mismatch_ready_dict[domain][clust].append([align_len, ['syn', ind]])

    write_list=list()
    write_list.append(['cluster', 'domain','coord','site_type','align_length', 'clust_abundance'])


    for domain in mismatch_ready_dict:
        for cluster in mismatch_ready_dict[domain]:
            num_alns=len(alignment_dict_prots[domain][cluster])
            for site in mismatch_ready_dict[domain][cluster]:
                length =site[0]
                rel_site=site[1][1]/site[0]*100
                if domain=='domain2':
                    rel_site += 100
                elif domain=='domain3':
                    rel_site+=200
                write_list.append([cluster, domain,rel_site,site[1][0], length, num_alns])


    write_csv(write_list,'multiclusters_mismacth_sites.csv')
    #make_id_distr_for_multi(alignment_dict_nucs)


def count_conserve_sites(align_dir: str) -> list:
    """
        Summarizes per-site conservative scrore
    """
    inds_gap=defaultdict(dict)
    inds_no_gap=defaultdict(dict)

    gap_alignments = ['domain1.msa', 'domain2.msa', 'domain3.msa']
    #min_gap_aln = ['min_gaps_domain1.msa', 'min_gaps_domain2.msa', 'min_gaps_domain3.msa']
    mismatch_ready_dict=dict()

    cons_scores=[]
    cons_scores.append(['coord', 'domain','site_type','cons_score'])
    len_ind=0

    for fasta_file in gap_alignments:
        domain=fasta_file.split('.')[0]
        inds_gap[domain]=dict()

        for alseq in SeqIO.parse(os.path.join(os.path.realpath(align_dir),fasta_file),"fasta"): 
            for ind in range(len(str(alseq.seq))):

                if ind+len_ind not in  inds_gap[domain]:

                    inds_gap[domain][ind+len_ind]=[str(alseq.seq)[ind]]

                else:
                    inds_gap[domain][ind+len_ind].append(str(alseq.seq)[ind])

            new_ind=len(str(alseq.seq))

        len_ind+=new_ind

    for domain in inds_gap:
        for ind in inds_gap[domain]:
            inds_gap[domain][ind]=dict(Counter(inds_gap[domain][ind]))

    for domain in inds_gap:
        for ind in inds_gap[domain]: 
            max_count=0
            for symb in inds_gap[domain][ind]:
                if symb!='-':
                    if inds_gap[domain][ind][symb]>max_count:
                        max_count=inds_gap[domain][ind][symb]

            cons_scores.append([ind, domain, 'gap_aln', max_count/367])

    #len_ind=0
    #for fasta_file in min_gap_aln:
    #    domain=fasta_file.split('.')[0].split('_')[2]
    #    inds_no_gap[domain]=dict()

    #    for alseq in SeqIO.parse(os.path.join(os.path.realpath(align_dir),fasta_file),"fasta"): 
    #        for ind in range(len(str(alseq.seq))):
    #            if ind+len_ind not in  inds_no_gap[domain]:

    #                inds_no_gap[domain][ind+len_ind]=[str(alseq.seq)[ind]]

    #            else:
    #                inds_no_gap[domain][ind+len_ind].append(str(alseq.seq)[ind])

    #    new_ind=len(str(alseq.seq))

       
    #    len_ind+=new_ind

    #for domain in inds_no_gap:
    #    for ind in inds_no_gap[domain]:
    #        inds_no_gap[domain][ind]=dict(Counter(inds_no_gap[domain][ind]))

    #for domain in inds_no_gap:
    #    for ind in inds_no_gap[domain]: 
    #        max_count=0

    #        for symb in inds_no_gap[domain][ind]:
    #            if symb!='-':
    #                if inds_no_gap[domain][ind][symb]>max_count:
    #                    max_count=inds_no_gap[domain][ind][symb]

    #        cons_scores.append([ind, domain, 'min_gap_aln', max_count/367])
    
    write_csv(cons_scores,'conservative_scores.csv')

def count_conserve_sites_no_boundaries(align_dir: str) -> list:
    """
        Summarizes per-site conservative scrore for the whole Cry alignment
    """
    gap_alignment = 'all_sequences_nucl.aln.fasta'
    #min_gap_aln = '733_names_extracted.aln.fasta'
    mismatch_ready_dict=dict()

    cons_scores=[]
    cons_scores.append(['coord', 'site_type','cons_score'])
    inds_gap=defaultdict(list)

    for record in SeqIO.parse(os.path.join(os.path.realpath(align_dir),gap_alignment),"fasta"):
        for ind in range(len(str(record.seq))):
            if ind not in  inds_gap:
                inds_gap[ind]=[str(record.seq)[ind]]
            else:
                inds_gap[ind].append(str(record.seq)[ind])

    for ind in inds_gap:
        inds_gap[ind]=dict(Counter(inds_gap[ind]))


    for ind in inds_gap: 
        max_count=0
        for symb in inds_gap[ind]:
            if symb!='-':
                if inds_gap[ind][symb]>max_count:
                    max_count=inds_gap[ind][symb]

        cons_scores.append([ind, 'gap_aln', max_count/732])

    #inds_gap_min=defaultdict(list)
    #for record in SeqIO.parse(os.path.join(os.path.realpath(align_dir),min_gap_aln),"fasta"):
    #    for ind in range(len(str(record.seq))):
    #        if ind not in  inds_gap:
    #            inds_gap_min[ind]=[str(record.seq)[ind]]
    #        else:
    #            inds_gap_min[ind].append(str(record.seq)[ind])

    #for ind in inds_gap_min:
    #    inds_gap_min[ind]=dict(Counter(inds_gap_min[ind]))


    #for ind in inds_gap_min: 
    #    max_count=0
    #    for symb in inds_gap_min[ind]:
    #        if symb!='-':
    #            if inds_gap_min[ind][symb]>max_count:
    #                max_count=inds_gap_min[ind][symb]

    #    cons_scores.append([ind, 'min_gap_aln', max_count/732]) #367


    write_csv(cons_scores,'733_conservative_scores_full_aln.csv')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='agregates multiple sequence alignments from clusters with more than 10 sequences')
    parser.add_argument('-al', '--al_d', dest='align_dir', help='path to the alignments\' directory',
                        type=str)
    #parser.add_argument('-e', '--evol', dest='evol_dir', help='directory with selection results',
    #                    type=str)
    #parser.add_argument('-s', '--sum', dest='sum_evol', help='the path to the summary table with selection results',
    #                    type=str)
    parser.add_argument('-n', '--nuc', dest='nuc_evol', help='the path to nucleotide alignments',
                        type=str)
    args = parser.parse_args()
    
    #aggregate_sites(args.evol_dir, args.sum_evol)
    expand_annotated_table(args.align_dir)
    #count_conserve_sites(args.nuc_evol)
    count_conserve_sites_no_boundaries(args.nuc_evol)

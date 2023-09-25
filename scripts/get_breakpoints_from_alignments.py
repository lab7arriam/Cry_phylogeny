#/usr/bin/python3.7
import click
import os
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from make_data_for_heatmap_universal import find_sequence_identity, make_sequence_dict
from annotate_clusters_consistency import write_csv, read_csv_to_list
import pandas as pd

def make_beds(bed_tab, asmbl_tab):

    asmbl_tab_raw=read_csv_to_list(os.path.realpath(asmbl_tab),delim='\t' ,headless=False)
    asmbl_dict=defaultdict(list)

    for asmbl_row in asmbl_tab_raw:
        asmbl_dict[asmbl_row[4]]=asmbl_row[5].split('/')[9]

    coord_tab_raw=read_csv_to_list(os.path.realpath(bed_tab),delim='\t' ,headless=False)
    bed_dict=defaultdict(list)

    for row in coord_tab_raw:
        bed_dict[row[4]].extend([[row[5], row[7], row[8],row[3]]])

    for assembly in bed_dict:
        write_csv(bed_dict[assembly], os.path.realpath(asmbl_dict[assembly]+'.bed'))


def get_parent_names(proc_str):
    try:
        rec_str = proc_str.split('(')[1].replace(')','').split(';')
    except:
        rec_str = proc_str.replace(')','').split(';')

    return([r.strip() for r in rec_str])



def fasta_file_to_dict(fasta_file):
    ret_dict=dict()

    for record in SeqIO.parse(os.path.realpath(fasta_file),"fasta"):
        record.seq=Seq(str(record.seq.upper()))
        ret_dict[record.id]=record

    return(ret_dict)


def make_beds_dict_from_bed(bed_file, bed_type='regular'):
    bed_dict=defaultdict(dict)
    for row in read_csv_to_list(bed_file, headless=False,delim='\t'):
        if 'N_term' not in row and 'C_term' not in row:
           if bed_type=='regular':
               bed_dict[row[0]][row[3]]=[int(x) for x in row[1:3]]
           else:
               bed_dict[row[0]]['domain1']=[int(x)*3+1 for x in row[1:3]]
               bed_dict[row[0]]['domain2']=[int(x)*3+1 for x in row[3:5]]
               bed_dict[row[0]]['domain3']=[int(x)*3+1 for x in row[5:7]]

    return(bed_dict)

def get_major_and_minor_parents(row, mode='mearged'):

    if row[1]['Type']=='domain1':
        min=get_parent_names(row[1]['parent_1'])
        min = [m.strip() for m in min]

        maj1=get_parent_names(row[1]['parent_2'])
        maj1 = [m.strip() for m in maj1]

        maj2=get_parent_names(row[1]['parent_3'])
        maj2 = [m.strip() for m in maj2]
            


    if row[1]['Type']=='domain2':
        min=get_parent_names(row[1]['parent_2'])
        min = [m.strip() for m in min]

        maj1=get_parent_names(row[1]['parent_1'])
        maj1 = [m.strip() for m in maj1]

        maj2=get_parent_names(row[1]['parent_3'])
        maj2 = [m.strip() for m in maj2]

    if row[1]['Type']=='domain3':
        min=get_parent_names(row[1]['parent_3'])
        min = [m.strip() for m in min]

        maj1=get_parent_names(row[1]['parent_1'])
        maj1 = [m.strip() for m in maj1]

        maj2=get_parent_names(row[1]['parent_2'])
        maj2 = [m.strip() for m in maj2]

    if mode=='mearged':
        return([list(set(min)), list(set(maj1+maj2))])

    elif mode=='intersect':
        return([list(set(min)), list(set(maj1+maj2)), list(set(maj1).intersection(set(maj2)))])


def get_start_stop_easy(part_seq, full_seq, aln):
    start_in_full_seq = full_seq.find(part_seq)

    full_seq_before = full_seq[0:start_in_full_seq]
    seq_before_find = ''
    aln_start = 0
    aln_start_no_gaps = 0
    seq_within_find = ''

    if start_in_full_seq!=0: #new function is subseq starts with 0 index
        for symb in aln:
            aln_start+=1
            if symb=='-':
                continue
            else:
                seq_before_find+=symb
                aln_start_no_gaps+=1
            if seq_before_find==full_seq_before:
                break

    aln_stop = aln_start
    aln_stop_no_gaps = aln_start_no_gaps


    for symb in aln[aln_start:]:
        aln_stop+=1
        if symb=='-':
            continue
        else:
            seq_within_find+=symb
            aln_stop_no_gaps+=1
        if seq_within_find==part_seq:
            break

    return((aln_start, aln_stop))


def analyze_tox_list_for_blocks(tox_list, rec_seq, recomb_coords_full, seq_dist, ID):
    ret_list=[]
    for tox in tox_list:

        alns = find_sequence_identity((rec_seq.seq, seq_dist[tox].seq), gap_mode='local', get_alns='yes')

        rec_seq_part = rec_seq.seq[recomb_coords_full[0]:recomb_coords_full[1]]
        start, stop = get_start_stop_easy(rec_seq_part, rec_seq.seq, alns[0])
        ret_list.append([tox,alns[1], start, stop])
    return(ret_list)

def extract_flanking_regions(aln_seq, start, step, length):

    flank_seq=''
    if step==1:
        for coord in range(start, len(aln_seq), step):
            if len(flank_seq)==length:
                break
            if aln_seq[coord]=='-':
                pass
            else:
                flank_seq+=aln_seq[coord]

    else:
        for coord in range(start, 0, step):
            if len(flank_seq)==length:
                break
            if aln_seq[coord]=='-':
                pass
            else:
                flank_seq+=aln_seq[coord]

    return(flank_seq)

def make_parents_dicts(line):
    minor_dom_dicts = defaultdict(list)
    major_dom_dicts = defaultdict(list)

    domain=line[1]['Type'].replace('_partial', '')
    domain_mappings_dict={'parent_1':'domain1','parent_2':'domain2','parent_3':'domain3'}

    for parent_type in ['parent_1','parent_2','parent_3']:
        for parent in get_parent_names(line[1][parent_type]):

            if domain_mappings_dict[parent_type]==domain:
                minor_dom_dicts.setdefault(parent, []).append(domain_mappings_dict[parent_type])
            else:
                major_dom_dicts.setdefault(parent, []).append(domain_mappings_dict[parent_type])

    return(minor_dom_dicts, major_dom_dicts)

def unpacking_list_of_lists(list_of_lists):
    total_list = []
    for list_sub in list_of_lists:
        for el in list_sub:
            total_list.append(el)
    return(total_list)

def make_row_with_domain_mappings(ID, domain, recs, rec_seq, recomb_coords_full, recomb_coords_proc, coord_pair, full_beds, proc_beds, seq_dist_proc, minor_dom_dicts, major_dom_dicts, seq_dist, par_type):
    ref_rec = recs[0]
    recs_name=';'.join(recs)
    Par_name = coord_pair[0]

    #Get coordinates and domain data for recombinants
    rec_start, rec_stop = recomb_coords_proc
    doms_cord_rec_proc = unpacking_list_of_lists([proc_beds[ref_rec][dom] for dom in proc_beds[ref_rec].keys()])

    rec_start_full, rec_stop_full = recomb_coords_full
    doms_coord_rec_full = unpacking_list_of_lists([full_beds[ref_rec][dom] for dom in full_beds[ref_rec].keys()])

    num_gaps_before = coord_pair[1][:coord_pair[2]].count('-')
    num_gaps_after = coord_pair[1][coord_pair[2]:coord_pair[3]].count('-')

    #print(coord_pair)

    coords_parent_fixed = [coord_pair[2]-num_gaps_before, coord_pair[3]-num_gaps_after-num_gaps_before]
    if coords_parent_fixed[1] > len(str(seq_dist[Par_name].seq)):
        coords_parent_fixed[1] = len(str(seq_dist[Par_name].seq))

    coords_parent_dom_only = [coords_parent_fixed[0]-full_beds[Par_name]['domain1'][0],coords_parent_fixed[1]-full_beds[Par_name]['domain1'][0]]


    #Calculate potential offset
    d3_dif_between_full_proc = coords_parent_fixed[1] - full_beds[Par_name]['domain3'][1]
    d1_dif_between_full_proc = full_beds[Par_name]['domain1'][0] - coords_parent_fixed[0] 
  
    #Adjust coordinates
    if d1_dif_between_full_proc>0:
        coords_parent_dom_only[0]+=d1_dif_between_full_proc
        coords_parent_dom_only[1]+=d1_dif_between_full_proc

    if d3_dif_between_full_proc>0:
        coords_parent_dom_only[0]-=d3_dif_between_full_proc
        coords_parent_dom_only[1]-=d3_dif_between_full_proc

    #Get coordinates and domain data for parents
    par_start, par_stop = coords_parent_dom_only
    doms_cord_par_proc = unpacking_list_of_lists([proc_beds[Par_name][dom] for dom in proc_beds[Par_name].keys()])

    par_start_full, par_stop_full = coords_parent_fixed
    doms_coord_par_full = unpacking_list_of_lists([full_beds[Par_name][dom] for dom in full_beds[Par_name].keys()])

    rec_proc = seq_dist_proc[recs[0]][recomb_coords_proc[0]:recomb_coords_proc[1]].seq
    rec_full = rec_seq.seq[ rec_start_full:rec_stop_full]
    
    par_proc = seq_dist_proc[Par_name][par_start:par_stop].seq
    par_full = seq_dist[Par_name][par_start_full:par_stop_full].seq

    if d1_dif_between_full_proc<0:
        d1_dif_between_full_proc=0

    if d3_dif_between_full_proc<0:
        d3_dif_between_full_proc=0

    ret_list=[]

    if par_type=='minor':
        ret_list.append([ID, recs_name, Par_name, domain, minor_dom_dicts[Par_name][0],par_type, rec_start, rec_stop, rec_start_full, rec_stop_full]+ doms_cord_rec_proc+doms_coord_rec_full+ [par_start, par_stop, par_start_full, par_stop_full] +doms_cord_par_proc +doms_coord_par_full)
    else:
        for par_domain in major_dom_dicts[Par_name]:
            ret_list.append([ID, recs_name, Par_name, domain, par_domain,par_type, rec_start, rec_stop, rec_start_full, rec_stop_full]+ doms_cord_rec_proc+doms_coord_rec_full+[par_start, par_stop, par_start_full, par_stop_full] +doms_cord_par_proc +doms_coord_par_full)

    return(ret_list)


def get_conservative_block_for_rec_parents(rec_tab, seq_dist, dom_seq_dict_all, full_beds, proc_beds, seq_dist_proc):

    blocks_stat=[]
    blocks_stat.append(['ID', 'Recs', 'Mins', 'Majs','Par1', 'Par2' ,'Dom1_id', 'Dom2_id', 'Dom3_id', 'Full_id', 'Proc_id', 'Block1_id', 'Block2_id'])

    file_reader = read_csv_to_list(os.path.realpath(rec_tab), headless=False)
    rec_df = pd.DataFrame(file_reader[1:], columns=file_reader[0])

    flanking_seqs = []

    coord_results=[['ID', 'Recs', 'Parent', 'Domain', 'Curr_Domain','Parent_type', 'r_start_proc', 'r_stop_proc', 'r_start_full', 'r_stop_full',
                        'd1_r_start_proc',  'd1_r_stop_proct','d2_r_start_proc',  'd2_r_stop_proc','d3_r_start_proc',  'd3_r_stop_proc', 
                        'd1_r_start_full',  'd1_r_stop_full','d2_r_start_full',  'd2_r_stop_full','d3_r_start_full',  'd3_r_stop_full', 
                        'p_start_proc', 'p_stop_proc', 'p_start_full', 'p_stop_full',
                        'd1_p_start_proc',  'd1_p_stop_proct','d2_p_start_proc',  'd2_p_stop_proc','d3_p_start_proc',  'd3_p_stop_proc', 
                        'd1_p_start_full',  'd1_p_stop_full','d2_p_start_full',  'd2_p_stop_full','d3_p_start_full',  'd3_p_stop_full']]


    for line in rec_df.iterrows(): 
        print(coord_results)

        line[1]['Type']= line[1]['Type'].replace('_partial', '')
        domain=line[1]['Type'].replace('_partial', '')
        ID=line[1]['ID']
        recs=get_parent_names(line[1]['Rec'])
        parents_lists = get_major_and_minor_parents(line)

        ref_rec = recs[0]
        recomb_coords_proc = [int(line[1]['start']),int(line[1]['stop'])]
        rec_proc_break = seq_dist_proc[ref_rec][recomb_coords_proc[0]:recomb_coords_proc[1]]

        rec_seq=seq_dist[ref_rec]

        rec_full_break_start = str(rec_seq.seq).find(str(rec_proc_break.seq))
        rec_full_break_stop = rec_full_break_start+len(rec_proc_break)
        recomb_coords_full = [rec_full_break_start, rec_full_break_stop]
        

        min_parent_list=[tox for tox in parents_lists[0] if tox!='unknown']
        maj_parent_list=[tox for tox in parents_lists[1] if tox!='unknown']


        event_type ='Full'
        event_type='Partial'
        minor_dom_dicts, major_dom_dicts = make_parents_dicts(line)


        for rec in recs:
            curr_rec_seq = seq_dist[rec]
            rec_left_flank=extract_flanking_regions(curr_rec_seq, recomb_coords_full[0], -1, 30)
            rec_right_flank=extract_flanking_regions(curr_rec_seq, recomb_coords_full[1], 1, 30)
            rec_ID_left=rec+'|'+ID+'|'+domain+'|left|rec'
            rec_ID_right=rec+'|'+ID+'|'+domain+'|right|rec'

            if len(min_parent_list)==0 or len(maj_parent_list)==0:
                rec_ID_left+='|unknown'
                rec_ID_right+='|unknown'

            left_seq_req = SeqRecord(Seq(str(rec_left_flank)),
                               id=rec_ID_left,
                               description='')
            right_seq_req = SeqRecord(Seq(str(rec_right_flank)),
                               id=rec_ID_right,
                               description='')

            flanking_seqs.append(left_seq_req)
            flanking_seqs.append(right_seq_req)

        if len(min_parent_list)>0:

            coords_list = analyze_tox_list_for_blocks(min_parent_list, rec_seq, recomb_coords_full, seq_dist, ID)
            for coord_pair in coords_list:
                domain_map_res = make_row_with_domain_mappings(ID, domain, recs, rec_seq, recomb_coords_full, recomb_coords_proc, 
                                             coord_pair, full_beds, proc_beds, seq_dist_proc, minor_dom_dicts, major_dom_dicts,seq_dist, 'minor')

                coord_results.extend(domain_map_res)

                left_flank = extract_flanking_regions(coord_pair[1], coord_pair[2], -1, 30)
                right_flank = extract_flanking_regions(coord_pair[1], coord_pair[3], 1, 30)

                par_ID_left=coord_pair[0]+'|'+ID+'|'+domain+'|left|min'
                par_ID_right=coord_pair[0]+'|'+ID+'|'+domain+'|right|min'

                if len(maj_parent_list)==0:
                    par_ID_left+='|unknown'
                    par_ID_right+='|unknown'

                left_seq_req_par = SeqRecord(Seq(str(left_flank)),
                                   id=par_ID_left,
                                   description='')
                right_seq_req_par = SeqRecord(Seq(str(right_flank)),
                                   id=par_ID_right,
                                   description='')
                flanking_seqs.append(left_seq_req_par)
                flanking_seqs.append(right_seq_req_par)
        
        
        if len(maj_parent_list)>0:
            #continue
            coords_list = analyze_tox_list_for_blocks(maj_parent_list, rec_seq, recomb_coords_full, seq_dist, ID)
            for coord_pair in coords_list:
                domain_map_res = make_row_with_domain_mappings(ID, domain, recs, rec_seq, recomb_coords_full, recomb_coords_proc, 
                                             coord_pair, full_beds, proc_beds, seq_dist_proc, minor_dom_dicts, major_dom_dicts,seq_dist, 'major')
                coord_results.extend(domain_map_res)

                left_flank = extract_flanking_regions(coord_pair[1], coord_pair[2], -1, 30)
                right_flank = extract_flanking_regions(coord_pair[1], coord_pair[3], 1, 30)
                par_ID_left=coord_pair[0]+'|'+ID+'|'+domain+'|left|maj'
                par_ID_right=coord_pair[0]+'|'+ID+'|'+domain+'|right|maj'

                if len(min_parent_list)==0:
                    par_ID_left+='|unknown'
                    par_ID_right+='|unknown'

                left_seq_req_par = SeqRecord(Seq(str(left_flank)),
                                   id=par_ID_left,
                                   description='')
                right_seq_req_par = SeqRecord(Seq(str(right_flank)),
                                   id=par_ID_right,
                                   description='')
                flanking_seqs.append(left_seq_req_par)
                flanking_seqs.append(right_seq_req_par)
        
    write_csv(coord_results, 'Mappings_for_breakpoints_de_novo_search_partials.csv')

    left_seqs_all=[rec for rec in flanking_seqs if 'left' in rec.id]
    right_seqs_all=[rec for rec in flanking_seqs if 'right' in rec.id]

    left_seqs_pars=[rec for rec in flanking_seqs if 'left' in rec.id and 'unknown' not in rec.id]
    right_seqs_pars=[rec for rec in flanking_seqs if 'right' in rec.id and 'unknown' not in rec.id]
        
    #SeqIO.write(left_seqs_all, 'left_with_uknowns.fasta', 'fasta')
    #SeqIO.write(right_seqs_all, 'right_with_uknowns.fasta', 'fasta')
    #SeqIO.write(left_seqs_pars, 'left_no_uknowns.fasta', 'fasta')
    #SeqIO.write(right_seqs_pars, 'right_no_uknowns.fasta', 'fasta')

    

@click.command()           
@click.option('--cry_tab', '-c', help="the file with assembiles' statistics", 
              type=click.Path(exists=True), metavar='<PATH>')   
@click.option('--asbml_tab', '-a', help="the file with assembiles' statistics", 
              type=click.Path(exists=True), metavar='<PATH>')
#@click.option('--fasta_file', '-f', help="the fasta file", 
#              type=click.Path(exists=True), metavar='<PATH>')  #for deleting ambiguous nucleotides
@click.option('--all_nucl', '-an', help="the path to the mearged nucleotide sequences", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--proc_nucl', '-pn', help="the path to the processed nucleotide sequences", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--rec_tab', '-rt', help="'the path to the raw RDP table", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--all_doms', '-ad', help="the path to the domains of nucleotide sequences", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--full_bed', '-fb', help="the path to the domains mappings for full sequences", 
              type=click.Path(exists=True), metavar='<PATH>')
@click.option('--proc_bed', '-pb', help="the path to the domains mappings for processed sequences", 
              type=click.Path(exists=True), metavar='<PATH>')



def main(cry_tab, asbml_tab, all_nucl, rec_tab, all_doms, full_bed, proc_bed, proc_nucl): #,fasta_file
    #make_beds(cry_tab, asbml_tab)
    #clean_from_ambig_sites(fasta_file)
    seq_dist = fasta_file_to_dict(all_nucl)
    seq_dist_proc = fasta_file_to_dict(proc_nucl)
    dom_seq_dict_all = make_sequence_dict(all_doms)
    full_beds = make_beds_dict_from_bed(full_bed, bed_type='regular')
    proc_beds = make_beds_dict_from_bed(proc_bed, bed_type='regular')

    get_conservative_block_for_rec_parents(rec_tab, seq_dist, dom_seq_dict_all, full_beds, proc_beds, seq_dist_proc)

if __name__ == '__main__':
   main()

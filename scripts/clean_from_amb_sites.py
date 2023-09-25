
import click
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter


global ambig_sites   
ambig_sites=['M','K','R','Y','W','S','B','H','D','V', 'N','U',
             'm','k','r','y','w','s','b','h','d','v', 'n','u']

def make_symbols_dict(fasta: str) -> defaultdict(dict):
    ##Parses fasta sequences
    ##Returns a dictionary with the dictionary of symbols' frequency per coordinate in the alignment
   

    #Make lists of sequences in string format for each sample
    seqs_lists = list(map(lambda x: list(str(x.seq)), SeqIO.parse(os.path.realpath(fasta), "fasta")))

    #Transform each sequence into a Counter object with sybmols' frequency
    Counters = list(map(lambda x:Counter(x), zip(*seqs_lists)))

    #Transform into dictionary with coordinates
    per_symbol_dict = {i:dict(Counters[i]) for i in range(len(Counters))}

    return(per_symbol_dict)


def get_best_symbol(symb_counter):
    ret_code = 1 
    symb_counter = dict(sorted(symb_counter.items(), key=lambda item: item[1],  reverse=True))
    for symb_key in symb_counter:
        if symb_key not in ambig_sites and symb_key!='-':
            ret_code = 0
            return(str(symb_key))
    if ret_code!=0:
        #print(symb_counter)
        return('N')

def clean_from_ambig_sites(fasta_file, fasta_symbols_dict):
    clean_recs=[]

    for record in SeqIO.parse(os.path.realpath(fasta_file),"fasta"):
        new_rec=record
        new_seq=list(str(record.seq))
        for base_ind in range(len(new_seq)):
            if new_seq[base_ind] in ambig_sites:
                new_seq[base_ind] = get_best_symbol(fasta_symbols_dict[base_ind])
        new_rec.seq=Seq(''.join(new_seq).upper())
        clean_recs.append(new_rec)
    
            
    SeqIO.write(clean_recs,os.path.realpath(fasta_file).replace('.fasta','.fixed.fasta').replace('.aln','.fixed.fasta'), "fasta")



@click.command()
@click.option('--fasta_file', '-f', help="the fasta file", 
              type=click.Path(exists=True), metavar='<PATH>') 


def main(fasta_file):
    print('Creating symbols dictionary')
    fasta_symbols_dict=make_symbols_dict(os.path.realpath(fasta_file))
    print('Cleaning from ambiguous sites')
    clean_from_ambig_sites(fasta_file, fasta_symbols_dict)

if __name__ == '__main__':
   main()

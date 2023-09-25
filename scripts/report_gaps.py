from Bio import pairwise2 as pw2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein, generic_dna
from Bio.SeqRecord import SeqRecord
import sys
import csv
from collections import defaultdict


prot_al=defaultdict(list)
nucl_fasta=dict()

for record in SeqIO.parse(sys.argv[1],"fasta"):
    prot_al[record.id]=[str(record.seq).replace('-',''),str(record.seq)]
for record in SeqIO.parse(sys.argv[2],"fasta"):
    nucl_fasta[record.id]=str(record.seq)

#print(Seq(nucl_fasta['BT_Cry1Jc1'],generic_dna).translate())


#global_align1 = pw2.align.globalxx(Seq(nucl_fasta['BT_Cry1Jc1'],generic_dna).translate(), prot_al['BT_Cry1Jc1'][0])
#for i in range(len(global_align1[0][1])):
#    if global_align1[0][1][i]!='-':
#        if global_align1[0][1][i]!=Seq(nucl_fasta['BT_Cry1Jc1'],generic_dna).translate()[i]:
#            print(global_align1[0][1][i],Seq(nucl_fasta['BT_Cry1Jc1'],generic_dna).translate()[i])

reported_gaps=[]
print('reporting')

for tox in nucl_fasta:
    global_align1 = pw2.align.globalxx(Seq(nucl_fasta[tox],generic_dna).translate(), prot_al[tox][0])
    codes_dict=dict()
    code_ind=0
    for i in range(len(global_align1[0][1])):
        if global_align1[0][1][i]!='-':
            codes_dict[code_ind]=str(Seq(nucl_fasta[tox],generic_dna)[i*3:i*3+3])
            code_ind+=1
            #print(tox, i, global_align1[0][1][i],Seq(nucl_fasta[tox],generic_dna)[i*3:i*3+3],Seq(nucl_fasta[tox],generic_dna).translate()[i],Seq(nucl_fasta[tox],generic_dna)[i*3:i*3+3].translate() )
            #if global_align1[0][1][i]!=Seq(nucl_fasta[tox],generic_dna).translate()[i]:
             #   print(global_align1[0][1][i],Seq(nucl_fasta[tox],generic_dna).translate()[i])
    new_code_ind=0
    reported_str=''
    for symb in prot_al[tox][1]:
        if symb!='-':
            reported_str=reported_str+codes_dict[new_code_ind]
            new_code_ind+=1
        else: 
            #print('rep')
            reported_str=reported_str+'---'
    reported_gaps.append(SeqRecord(Seq(str(reported_str),generic_dna),id=tox,description=''))
    #print(tox, reported_str,prot_al[tox][1],prot_al[tox][0],nucl_fasta[tox])
    #print(tox, len(Seq(nucl_fasta[tox],generic_dna).translate()), len(prot_al[tox][0]), len(nucl_fasta[tox]), len(reported_str), Seq(nucl_fasta[tox],generic_dna).translate(), prot_al[tox][1] , global_align1[0][1])
    if  len(Seq(nucl_fasta[tox],generic_dna).translate()) != len(prot_al[tox][0]):
        print('unequal_lengths', tox, len(Seq(nucl_fasta[tox],generic_dna).translate()), len(prot_al[tox][0]))
    if '-' in global_align1[0][1]:
        print(tox, 'non perfect match for translation')
    if 'N 'in Seq(nucl_fasta[tox],generic_dna).translate():
        print(tox, 'N in align')
SeqIO.write(reported_gaps,sys.argv[2].split('.')[0]+'_reported.fasta',"fasta")

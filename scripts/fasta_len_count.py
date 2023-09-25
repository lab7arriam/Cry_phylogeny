import sys
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord


inputfile = sys.argv[1]

inlist=list()

for record in SeqIO.parse(inputfile,"fasta"):
    print(str(record.id) + '\t'+ str(len(record.seq)))
    #print(str(inputfile.split('/')[0]).replace('.fasta','') + '\t'+ str(len(record.seq)))
    #break

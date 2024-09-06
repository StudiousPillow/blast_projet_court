from Bio import Entrez
from Bio import SeqIO

Entrez.email = "pouet"
with Entrez.esearch(db='nucleotide', term='SARS-CoV-2', retmax=15) as handle:
    record = Entrez.read(handle)
id_list = record['IdList']
with open('data/Sars-CoV-2_sequences.fasta','w') as output:
    for seq_id in id_list:
        with Entrez.efetch(db='nucleotide', id=seq_id, rettype='fasta',retmode='text') as handle:
            seq_record = SeqIO.read(handle, 'fasta')
            print(seq_record)
            print("-----")
            SeqIO.write(seq_record, output, 'fasta')

with Entrez.esearch(db='nucleotide', term='SARS-CoV-2', retmax=15) as handle:
    record = Entrez.read(handle)
id_list = record['IdList']
with open('data/Sars-CoV-2_sequences.fasta','w') as output:
    for seq_id in id_list:
        with Entrez.efetch(db='nucleotide', id=seq_id, rettype='fasta',retmode='text') as handle:
            seq_record = SeqIO.read(handle, 'fasta')
            print(seq_record)
            print("-----")
            SeqIO.write(seq_record, output, 'fasta')
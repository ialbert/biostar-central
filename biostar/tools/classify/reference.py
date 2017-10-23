# This script takes accession list as input and produce a multifasta reference file and accession-taxid table.
# Input : accession list file.
# output1 : sequence.fa (multifasta reference file corresponding to accessions.)
# output2 : acc2taxa.txt (tab delimited file where accession is first column and taxid is second column).


import sys,os
from Bio import Entrez
from Bio import SeqIO


def get_fasta(accessions):
    """
    Takes accession list as input and produce a multifasta reference file and accession-taxid table.
    Input : accession list
    outputs: sequence.fa, acc2taxa.txt
    """
    # create a dictionary of list where 1st element is sequence and second element is taxid.
    genome_info = dict()
    Entrez.email = "name@abc.com"
    for accession in accessions:
        print (f'working with {accession}')
        try:
            with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="fasta", id=accession) as handle:
                seq_record = SeqIO.read(handle, "fasta")
                species_name = ' '.join([seq_record.description.split()[1],seq_record.description.split()[2]])
                taxon_id = os.popen("echo {0} | taxonkit name2taxid".format(species_name)).read().split('\t')[1]
                genome_info[accession]=[seq_record.seq,taxon_id.strip()]
        except Exception:
            print(f'{accession} gave error.{accession} dicarded.')

        for key, val in genome_info.items():
            id = key
            seq = genome_info[key][0]
            taxid = genome_info[key][1]
            print(id, taxid)

    return


if __name__ =="__main__":
    #accessions = sys.argv[1]
    accessions=['AP013046','JN024756','NC_0008600000']
    get_fasta(accessions)

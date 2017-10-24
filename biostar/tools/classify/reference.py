# This script takes accession list as input and produce a multifasta reference file and accession-taxid table.
# Input : accession list file.
# output1 : sequence.fa (multifasta reference file corresponding to accessions.)
# output2 : acc2taxa.txt (tab delimited file where accession is first column and taxid is second column).


import sys,os
import argparse
from Bio import Entrez
from Bio import SeqIO


def accession_details(accessions):
    """
    Takes accessions as a list and produce a multifasta reference file and accession-taxid table.
    Input : accession list
    outputs: sequence.fa, acc2taxa.txt
    """
    # create a dictionary of list of the form
    #d = { 'acc' : [species_name, sequence, taxon_id] }

    accession_info = dict()
    Entrez.email = "name@abc.com"
    for accession in accessions:
        accession = accession.strip()
        try:
            with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="fasta", id=accession) as handle:
                seq_record = SeqIO.read(handle, "fasta")
                species_name = ' '.join([seq_record.description.split()[1],seq_record.description.split()[2]])
                accession_info[accession]=[species_name,seq_record.seq]

        except Exception:
            print(f'{accession} gave error.{accession} discarded.')

    if args.sequence:
        for k, v in accession_info.items():
            header = ">" + k
            seq = accession_info[k][1]
            print(header)
            print(seq)

    if args.taxid:
        for k, v in accession_info.items():
            species = accession_info[k][0]
            taxon_id = os.popen("echo {0} | taxonkit name2taxid".format(species)).read().split('\t')[1]
            accession_info[k].append(taxon_id.strip())
        for k, v in accession_info.items():
            taxon_id = accession_info[k][2]
            out = "\t".join([k,taxon_id])
            print(out)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--ids", help="File with accessions.")
    parser.add_argument("--taxid", action='store_true',help="Print accession - taxid table.")
    parser.add_argument("--sequence", action='store_true', help="Print fasta sequences of accessions.")
    args = parser.parse_args()

    if not args.ids:
        print("This requires a file with accessions")

    if args.ids:
        accession_list = open(args.ids).readlines()
        accession_details(accession_list)

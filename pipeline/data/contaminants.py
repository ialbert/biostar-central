import sys,csv


def get_contaminants(fname):
    primer_collect = dict()
    barcode_collect = dict()

    with open(fname) as csvfile:

        stream =csv.DictReader(csvfile,delimiter="\t",restkey='extra',fieldnames =
        ("sample_name", "sample_group","barcode","fwd_primer","rev_primer","target_gene"))
        next(stream)
        for row in stream:
            idx = 1
            fp = row['fwd_primer']
            primer_collect[fp] = "fwd_primer" +str(idx)
            rp = row['rev_primer']
            primer_collect[rp] = "rev_primer" + str(idx)
            bcode = row['barcode']
            barcode_collect[bcode] = "barcode" +str(idx)
            idx +=1

        return primer_collect, barcode_collect


def write_fasta(data, outfile):
    '''
    writes data in a dictionary to a file named fname
    '''
    for k, v in data.items():
        outfile.write(">" + v + "\n")
        outfile.write(k)
        outfile.write("\n")


if __name__ == "__main__":

    samplesheet = sys.argv[1]
    contaminant_type = sys.argv[2]  # contaminant type  can be primer, barcode or all
    outfile = sys.stdout
    primer_seqs, barcode_seqs = get_contaminants(samplesheet)

    # writes primers.fa
    if contaminant_type == "primer" and primer_seqs:
        write_fasta(primer_seqs, outfile)

    # write barcodes.fa
    if contaminant_type == "barcode" and  barcode_seqs:
        write_fasta(barcode_seqs, outfile)

    # write both
    if contaminant_type == "all":
        cont_seqs = {**primer_seqs, **barcode_seqs}
        write_fasta(cont_seqs, outfile)




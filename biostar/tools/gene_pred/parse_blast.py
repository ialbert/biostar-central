import csv,re


def parse_blast(blast_file):
    """
    parses blast tabular output to produce an annotation for each query sequence.
    It discards uncharacterized and hypothetical proteins and selects for hits from Bovine(Bos) genus.

    output is parsed blast file in tabular format.

    """
    # order of blast tabular output
    # qseqid sseqid stitle staxids pident qlen slen length qstart qend sstart send evalue bitscore score

    results=dict()
    header="\t".join(["id","protein","taxid","perc_iden","qlength","slength","aln_length","evalue"])
    print(header)

    with open(blast_file) as csvfile:
        stream = csv.reader(csvfile, delimiter="\t")
        for row in stream:
            if re.search("uncharacterized",row[2]):
                continue
            if re.search("hypothetical",row[2]):
                continue
            if re.search("Bos", row[2]):
                gid=row[0]
                if gid in results.keys():
                    continue
                taxid = row[3]
                protein = row[2]
                pid = row[4]
                qlen = row[5]
                slen = row[6]
                alen = row[7]
                evalue = row[12]
                results[gid] = [gid, protein, taxid, pid, qlen, slen, alen, evalue]

        for k, v in results.items():
            out = "\t".join(v)
            print(out)


if __name__ =="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Parses blast tabular output to produce annotation for query sequence.')

    parser.add_argument('--blast', help='File with blast tabular output')
    args = parser.parse_args()

    # Read the arguments.
    blast_file = args.blast

    parse_blast(blast_file)

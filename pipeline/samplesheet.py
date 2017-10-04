import sys,csv
import os,fnmatch
'''
creates a modified samplesheet by adding filename fields.
'''


class Bunch:
    pass


def update_sampleinfo(fname, datadir):
    '''
    reads sample info as a text file and add columns file1 and file2 to the samplesheet.
    returns the modified samplesheet as a list.
    '''
    data =[]

    with open(fname) as csvfile:

        stream =csv.DictReader(csvfile,delimiter="\t",restkey='extra',fieldnames =
        ("sample_name", "sample_group","barcode","fwd_primer","rev_primer","target_gene"))

        next(stream)

        for row in stream:
            if any(row[key] in (None, "") for key in row):
                sys.stderr.write("Error in {0}: Empty fields in samplesheet.\n\n".format(os.path.basename(__file__)))


            bunch = Bunch()

            bunch.sname = row['sample_name']
            bunch.sgroup = row['sample_group']
            bunch.tgene = row['target_gene']
            bunch.barcode = row['barcode']
            bunch.fwdprimer = row['fwd_primer']
            bunch.revprimer = row['rev_primer']

            file1, file2 = get_fnames(bunch.sname, datadir)
            bunch.file1 = file1
            bunch.file2 = file2
            data.append(bunch)
    return data


def get_fnames(sname , datadir):
    '''
    get R1 and R2 files matching samplename from "datadir"
    '''

    file1, file2 = None,None

    for file in os.listdir(datadir):
        R1patt = sname + "*R1*gz"
        R2patt = sname + "*R2*gz"

        if fnmatch.fnmatch(file, R1patt):
            file1 = file

        if fnmatch.fnmatch(file, R2patt):
            file2 = file

    if file1 is None or file2 is None:
        sys.stderr.write("File not found for {0}\n\n".format(sname))
    return file1, file2


if __name__ == "__main__":

    samplesheet = sys.argv[1]
    datadir =sys.argv[2]

    #samplesheet ="./test/sampleinfo.txt"
    #datadir ="test/data"

    data = update_sampleinfo(samplesheet, datadir)

    outfile = open("updated_sampleinfo.txt", "w")
    header ="\t".join(['sample_name','sample_group','target_gene','barcode','fwd_primer','rev_primer','file1','file2'])
    outfile.write(header + "\n")

    for bunch in data:
        row = "\t".join([bunch.sname,bunch.sgroup,bunch.tgene,bunch.barcode,bunch.fwdprimer,bunch.revprimer,bunch.file1,bunch.file2])
        outfile.write(row + "\n")
        #outfile.write("\n")

    outfile.close

    #if ERR_LIST:
    #    for e in ERR_LIST:
    #        print(e)



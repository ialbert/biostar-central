import csv,fnmatch,os,sys
import tarfile
import argparse
from const import INPUT_FILE_TYPE


class Bunch(object):
    pass


def setup_dir(dir):
    if not os.path.exists(dir):
        cmd='mkdir -p {0}'.format(dir)
        os.system(cmd)


def check_file(fname):
    if not os.path.isfile(fname):
        print("Error: {0} file not found".format(os.path.basename(fname)))
        sys.exit(1)


def read_sampleinfo(fname):
    '''
    reads metadata as a text file and returns a dictionary of metadata.
    '''
    data =[]

    with open(fname) as csvfile:
        stream =csv.DictReader(csvfile,delimiter="\t",restkey='extra',fieldnames =
        ( "sample_name","sample_group","target_gene", "barcode","fwd_primer","rev_primer" ))
        next(stream)
        for row in stream:
            bunch = Bunch()

            bunch.sname = row['sample_name'] if row['sample_name'] else "" #sample name
            bunch.sgroup = row['sample_group'] if row['sample_group'] else "" # tissue mix
            #bunch.primer = row['Primer Pair'] if row['Primer Pair'] else "" # primer pair
            bunch.tgene = row['target_gene'] if row['target_gene'] else "" # target gene
            bunch.barcode = row['barcode'] if row['barcode'] else ""
            bunch.fwdprimer = row['fwd_primer'] if row['fwd_primer'] else ""
            bunch.revprimer = row['rev_primer'] if row ['rev_primer'] else ""
            data.append(bunch)
    return data


def extract_files(fname,destdir):

    sampledir = destdir
    check_file(fname)
    print("extracting files....")
    tar = tarfile.open(fname, "r:gz")
    for tarinfo in tar:
        if tarinfo.isdir():
            sampledir = tarinfo.name
    tar.extractall(destdir)
    tar.close()
    return sampledir


def check_samples(sampleinfo, sampledir, datatype=None):
    found =1
    for bunch in sampleinfo:
        R1patt = "{0}_*R1*gz".format(bunch.sname)
        R2patt = "{0}_*R2*gz".format(bunch.sname)

        matched = find_files(sampledir,R1patt)
        if not matched:
            print("Read1 file for {0} is missing".format(bunch.sname))
            found = 0

        if datatype == "PE":
            matched = find_files(sampledir,R2patt)
            if not matched:
                print("Read2 file for {0} is missing".format(bunch.sname))
                found = 0

    if found == 1:
        print("Found all files. Ready to go.")
    else:
        print("\nError: Input files not found. Check before proceeding.\n")
        sys.exit(1)


def find_files(dir,patt):
    '''Return list of files matching pattern in base folder.'''

    return [n for n in fnmatch.filter(os.listdir(dir), patt) if
    os.path.isfile(os.path.join(dir, n))]


def option_parser():
        parser = argparse.ArgumentParser(description="Checks input files.")
        parser.add_argument("-i","--input", default=None, help="input data. Needs to be a gzipped file.")
        parser.add_argument("-s","--sampleinfo", default=None, help="sampleinfo and metadata. Needs to be a tab "
                                                                    "separated file.")
        parser.add_argument("-w","--workdir", default="./work", help="working directory.")
        return parser


if __name__ == "__main__":

    parser = option_parser()
    args = parser.parse_args()

    # check mandatory inputs
    if not (args.sampleinfo and args.input):
        parser.print_help()
        sys.exit(1)

    # read inputs
    sampleinfo = args.sampleinfo
    inputdata = args.input
    workdir = args.workdir

    # check if input files exist.
    if not os.path.basename(inputdata).endswith(INPUT_FILE_TYPE):
        print("Not tar zipped file; exiting")

    check_file(sampleinfo)
    check_file(inputdata)

    # set up directories.
    setup_dir(workdir)

    sampleinfo = read_sampleinfo(sampleinfo)

    # change sampledir to workdir/data
    data_dir = os.path.join(os.path.realpath(workdir),'data')

    if os.path.exists(data_dir):
        print("data directory exists. Removing..")
        os.system('rm -rf {0}'.format(data_dir))

    sampledir = extract_files(inputdata,workdir)

    if sampledir == workdir:
        setup_dir(os.path.join(sampledir,"data"))
        sampledir = os.path.join(sampledir,"data")
        os.system('mv {0}/*.* {1}'.format(workdir,sampledir))
    elif sampledir == "data":
        pass
    else:
        sampledir = os.path.join(workdir,sampledir)
        os.system('mv {0} {1}'.format(sampledir, data_dir))

    sampledir = os.path.join(workdir, "data")
    print(sampledir)
    check_samples(sampleinfo=sampleinfo, sampledir=sampledir, datatype="PE")

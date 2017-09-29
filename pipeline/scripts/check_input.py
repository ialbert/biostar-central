import csv,fnmatch,os,sys
import tarfile
import argparse


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


def read_metadata(fname):
    '''
    reads metadata as a text file and returns a dictionary of metadata.
    '''
    data =[]

    with open(fname) as csvfile:
        stream =csv.DictReader(csvfile,delimiter="\t")
        for row in stream:
            bunch = Bunch()
            bunch.sname = row['Sample Name']  # sample name
            bunch.tissmix = row['Tissue Mix']  # tissue mix
            bunch.primer = row['Primer Pair']  # primer pair
            bunch.tgene = row['Target gene']    # target gene
            bunch.barcode = row['Forward inline barcode']
            bunch.fwdprimer = row['Forward Primer']
            bunch.revprimer = row['Reverse Primer']
            data.append(bunch)
    return data


def extract_files(fname):

    check_file(fname)
    print("extracting files....")
    tar = tarfile.open(fname, "r:gz")
    for tarinfo in tar:
        if tarinfo.isdir():
            sampledir = tarinfo.name
    tar.extractall()
    tar.close()
    return sampledir


def check_samples(metadata, sampledir, datatype=None):
    found =1
    for bunch in metadata:
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
        parser.add_argument("-s","--sampleinfo", default=None, help="sample info and metadata. Needs to be a tab "
                                                                    "separated file.")
        parser.add_argument("-w","--workdir", default="./work", help="working directory.")
        return parser


if __name__ == "__main__":

    parser = option_parser()
    args = parser.parse_args()

    # check mandentory inputs
    if not (args.sampleinfo and args.input):
        parser.print_help()
        sys.exit(1)

    # read inputs
    metadata = args.sampleinfo
    inputdata = args.input
    workdir = args.workdir

    # check if input files exist.
    check_file(metadata)
    check_file(inputdata)

    # set up directories.
    setup_dir(workdir)

    metadata = read_metadata(metadata)

    if os.path.basename(inputdata).endswith("tar.gz"):   ## TO DO uncompress a gzipped file( not a tar file).
        sampledir = extract_files(inputdata)
    else:
        sampledir = inputdata

    # change sampledir to workdir/data
    data_dir = os.path.join(workdir,'data')

    if os.path.exists(data_dir):
        print("data directory exists. Removing..")
        os.system('rm -rf {0}'.format(data_dir))

    os.system('mv {0} {1}'.format(sampledir, data_dir))
    sampledir = "{0}/data".format(workdir)

    check_samples(metadata=metadata, sampledir=sampledir, datatype="PE")





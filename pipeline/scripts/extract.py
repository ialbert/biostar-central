import tarfile,sys,os

def check_file(fname):
    if not os.path.isfile(fname):
        print("Error: {0} file not found".format(os.path.basename(fname)))
        sys.exit(1)

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

if __name__ == "__main__":
	fname = sys.argv[1]
	check_file(fname)
	sampledir=extract_files(fname)
	print(sampledir)
		

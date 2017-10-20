import sys,csv,glob
import os,fnmatch

'''
creates a modified samplesheet by adding filename fields.
'''


class Bunch:
    def __init__(self, adict):
        for key, value in adict.items():
            setattr(self, key, value)

    def keys(self):
        return self.__dict__.keys()

    def items(self):
        return self.__dict__.items()

    def values(self):
        return self.__dict__.values()


def update_sampleinfo(fname,datadir):
    """
    Reads sample info as a text file and add additional columns to the samplesheet.
    Returns the modified samplesheet as a list.
    """
    data =[]

    with open(fname) as csvfile:

        stream =csv.DictReader(csvfile,delimiter="\t",restkey='extra',fieldnames =
        ("sample_name", "sample_group","barcode","fwd_primer","rev_primer","target_gene"))

        next(stream)

        for row in stream:

            if any(row[key] in (None, "") for key in row):
                sys.stderr.write("Error in {0}: Empty fields in samplesheet.\n\n".format(os.path.basename(__file__)))

            bunch = Bunch(row)

            patt1 = f'data/{bunch.sample_name}*R1*.gz'
            patt2 = f'data/{bunch.sample_name}*R2*.gz'

            def get_one(patt):
                files = glob.glob(patt)
                if len(files) != 1:
                    print(f"error for pattern {patt}")
                    sys.exit()
                return files[0]

            bunch.file1 = get_one(patt1)
            bunch.file2 = get_one(patt2)

            bunch.result1 = f'work/{bunch.sample_name}_result_R1.fq.gz'
            bunch.result2 = f'work/{bunch.sample_name}_result_R2.fq.gz'

            bunch.temp1 = f'work/{bunch.sample_name}_temp_R1.fq.gz'
            bunch.temp2 = f'work/{bunch.sample_name}_temp_R2.fq.gz'

            data.append(bunch)
    return data


if __name__ == "__main__":

    samplesheet = sys.argv[1]
    datadir = sys.argv[2]

    data = update_sampleinfo(samplesheet, datadir)

    first = data[0]

    header = "\t".join(first.keys())
    print(header)

    for bunch in data:
        line = "\t".join(bunch.values())
        print(line)


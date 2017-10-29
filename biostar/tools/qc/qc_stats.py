'''
This script takes an input fastq file and returns the no.of sequences in it.

'''

import sys, os, re

if __name__ == "__main__":
	sample = sys.argv[1]
	fname = sys.argv[2]
	cmd = "bioawk -c fastx \'END {{print NR }}\' {0}".format(fname)
	total = os.popen(cmd).read()
	total_pe = int(total) * 2
	print("\t".join([sample, str(total_pe)]))
'''
This script parses a bam file and produce the no. of reads mapped to each chromosome
after applying filters.
Filter remove alignments that are unmapped, supplementary, secondary and alignment lwngth not within specified
distance of read length.
'''

import pysam, sys
from collections import Counter

LENGTH_FILTER = 6


def parse_bam(fname, len_filter):
    bamfile = pysam.AlignmentFile(fname, "rb")
    stats = Counter()

    # Get chromosomes.
    chroms = []
    arr = bamfile.header['SQ']
    for chr in arr:
        chroms.append(chr['SN'])

    # get stats for each chromosome after applying filters.
    # Filter remove reads that are unmapped,supplementary, secondary and
    # those with alignment length not within query_length+-5.

    print("Chromosome\tMapped\tFiltered")
    for chr in chroms:
        iter = bamfile.fetch(chr)
        chr_total = chr + "_total"

        for read in iter:
            # total alignments in chromosome.
            stats[chr_total] += 1
            # Applying filters.
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
            if read.query_alignment_length in range(read.query_length - int(len_filter),
                                                    read.query_length + int(len_filter) + 1):
                stats[chr] += 1

        outline = "\t".join([chr, str(stats[chr_total]), str(stats[chr])])
        print(outline)


if __name__ == "__main__":
    fname = sys.argv[1]
    length_filter = sys.argv[2]
    # fname = "aligned.bam"
    # length_filter = 6
    parse_bam(fname, length_filter)

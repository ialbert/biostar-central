# Stop on any error.
set -uxe

# The table of contents of all data.
TOC={{reads.toc}}

# The assembly scaffolds.
REF={{reference.path}}

# How many of the largest scaffolds to keep.
TOPN={{topn.value}}

# How many reads to select randomly.
READ_NUM={{readnum.value}}

# The directory that holds the scaffolds.
IDX_DIR={{runtime.work_dir}}/fasta

# The name of the scaffold file.
IDX=$IDX_DIR/scaffolds_${TOPN}.fa

# The URL to the IGV genome.
IGV_GENOME={{runtime.job_url}}/fasta/scaffolds_${TOPN}.fa

# Index the reference so that we can extract the sequences.
mkdir -p $IDX_DIR

echo "Indexing with samtools."
samtools faidx $REF

# Find the largest scaffolds.
cat $REF | bioawk -c fastx  '{ print length($seq), $name }' | sort -k1,1rn | head -${TOPN} | cut -f 2 > largest.txt

# Truncate the file if it exists.
cat /dev/null >| $IDX

# Extract the scaffolds.
for ACC in  $(cat "largest.txt"); do
    echo "Extracting accession $ACC"
    samtools faidx $REF $ACC >> $IDX
done

# Index the selected sequences with samtools.
samtools faidx $IDX

# Build the index for the top N scaffolds.
bwa index ${IDX}

# The directory for the BAM files.
BAM={{runtime.work_dir}}/bam

# The URL for the bam files.
URL={{runtime.job_url}}/bam

# The directory to hold the READ samples.
FASTQ={{runtime.work_dir}}/fq

# Make the BAM directory.
mkdir -p $BAM $FASTQ

# Subselect only the paired end reads. Sample each with seqtk.
cat $TOC | sort | egrep "fq|fastq" | egrep "r1|r2|R1|R2" | parallel seqtk sample -2 {} $READ_NUM '>' $FASTQ/{/}.fq

# Process the input files by two.
ls -1 $FASTQ/*.fq | sort |  parallel -j 5 bwa mem ${IDX} {1} {2} '|' samtools sort '>' $BAM/{1/}.bam :::: - -

# Indexing alignment files.
ls -1 $BAM/*.bam | parallel -j 5 samtools index {}

# File with mapping stats.
MAPPED_STATS={{runtime.work_dir}}/mapping-stats.txt

# Create mapping statistics.
ls -1 $BAM/*.bam | parallel -j 5 "(echo {/} && samtools idxstats {})" >> $MAPPED_STATS

# Create IGV session for the data.
python -m biostar.tools.igv.bams --base $URL --bams $BAM --genome $IGV_GENOME > igv.xml

# File with total read counts.
READ_COUNTS={{runtime.work_dir}}/read-counts.txt

# Get total reads in each file.
cat $TOC | egrep "fq|fastq" | egrep "r1|r2|R1|R2" | parallel "echo {/} && bioawk -c fastx 'END{print NR}' {} " >>$READ_COUNTS

# Create barcharts with normalized mapped reads.
python -m biostar.tools.align.scaffold_plotter --mapped $MAPPED_STATS --total $READ_COUNTS --selected $READ_NUM >index.html

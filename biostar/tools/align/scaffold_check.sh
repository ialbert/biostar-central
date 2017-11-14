# Stop on any error.
set -ue

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

echo "Finding largest $TOPN scaffolds."
# Find the largest scaffolds.
cat $REF | bioawk -c fastx  '{ print length($seq), $name }' | sort -k1,1rn | head -${TOPN} | cut -f 2 > largest.txt

# Truncate the file if it exists.
cat /dev/null >| $IDX

echo "Extracting sequences for each scaffold."
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

echo "Subselecting $READ_NUM reads"
cat $TOC | egrep ".fq|.fastq" | parallel seqtk sample -2 {} $READ_NUM '>' $FASTQ/{/}.fq

echo  "Mapping reads to the genome."
ls -1 $FASTQ/*.fq | egrep ".fq|.fastq" | parallel -j 5 bwa mem ${IDX} {} '|' samtools sort '>' $BAM/{/}.bam

echo "Indexing alignment files."
ls -1 $BAM/*.bam | parallel -j 5 samtools index {}

# Create mapping statistics.
ls -1 $BAM/*.bam | parallel -j 5 "(echo {/} && samtools idxstats {})" >> {{runtime.work_dir}}/mapping-stats.txt

# Create IGV session for the data.
python -m biostar.tools.igv.bams --base $URL --bams $BAM --genome $IGV_GENOME > igv.xml

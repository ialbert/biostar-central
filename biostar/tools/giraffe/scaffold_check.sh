set -ue

# Get the table of contents.
TOC={{reads.toc}}

# The assembly scaffolds.
REF={{reference.path}}

# The index directory
IDX_DIR={{runtime.workdir}}/bwa

# How many scaffolds to keep.
TOPN={{topn.value}}

# How many reads to check
READ_NUM={{readnum.value}}

# The name of the index.
IDX=$IDX_DIR/scaffolds_${TOPN}.fa

# Make the index directory
mkdir -p $IDX_DIR

# Build index .
echo "Building the bwa index for the $TOPN scaffolds."

# Find the largest scaffolds.
cat $REF | bioawk -c fastx  '{ print length($seq), $name }' | sort -k1,1rn | head -${TOPN} | cut -f 2 > largest.txt

# Create a fasta index for the reference.
echo "Indexing with samtools."
samtools faidx $REF

# Extract each sequence in turn.
echo "Extracting sequences"
readarray -t ITEMS < largest.txt
for ACC in ${ITEMS[@]} ; do
    echo "Extracting accession $ACC"
    samtools faidx $REF $ACC >> $IDX
done

# Build the index for the top N scaffolds.
bwa index ${IDX}

# The directory to hold the BAM files.
BAM={{runtime.work_dir}}/bam

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

# Get total reads in each file.
#cat $TOC |egrep ".fq|.fastq" | parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" >>read_counts.txt


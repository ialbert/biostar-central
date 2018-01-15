# Stop on any error.
set -uxe

#
# Parameters.
#

# Assembly scaffolds.
ASSEMBLY={{assembly.value}}

# Species to be used as training set.
SPECIES={{species.value}}

# No. of processors to be used.
NPROC={{processors.value}}

#
# Preparing for Augustus.
#

# Size-sorted assembly fasta file.
ASSM_SORTED={{runtime.work_dir}}/assembly_sorted.fa

# Filtering for scaffolds >1kb length and sorting them by size .
{% verbatim %}
cat $ASSEMBLY | bioawk -c fastx 'length($seq)>1000 { print length($seq),$name,$seq } ' | sort -k1nr,1  | awk '{ print ">"$2"\n"$3"\n"}' | seqtk seq -l 80 - >$ASSM_SORTED
{% endverbatim %}

# Creating samtools index.
samtools faidx $ASSM_SORTED

# Augustus results directory.
AUGUSTUS={{runtime.work_dir}}/augustus

# Create AUGUSTUS directory.
mkdir -p $AUGUSTUS

# Augustus result files.
GENES=${AUGUSTUS}/genes.gtf
PROTEINS=${AUGUSTUS}/proteins.fa

#
# Running Augustus.
#

# Run augustus gene prediction.
mkdir -p tmp
{% verbatim %}
cat $ASSM_SORTED | parallel --j $NPROC --blocksize 5M --recstart '>' --pipe "cat {} > tmp/{%} && augustus --species=$SPECIES  tmp/{%}" > $GENES
{% endverbatim %}
rm -rf tmp

# Make proteins.fa from Augustus predicted protein sequences.
getAnnoFasta_mod.pl $GENES
mv ${AUGUSTUS}/genes.aa ${PROTEINS}

#
# Preparing for diamond blastp.
#

# Diamond NR database.
DIAMOND_NR=/export/refs/diamond-dbs/nr/nr

# Diamond protein-accession2taxon map file.
DIAMOND_TAXON=/export/refs/diamond-dbs/nr/prot.accession2taxid.gz

# Diamond blastp results directory.
DIAMOND={{runtime.work_dir}}/diamond

# Diamond blastp results.
DIAMOND_RES=${DIAMOND}/diamond-blastp.txt

# Create DIAMOND directory.
mkdir -p $DIAMOND

# Writing blastp header into output file.
HEADER="qseqid sseqid stitle staxids pident qlen slen length qstart qend sstart send evalue bitscore score"
echo $HEADER  |tr [:blank:] \\t >$DIAMOND_RES

#
# Running diamond blastp.
#

# Running diamond blastp on predicted proteins.
diamond blastp -f 6 $HEADER -d $DIAMOND_NR --taxonmap $DIAMOND_TAXON --max-target-seqs 15 -q $PROTEINS -p $NPROC >>$DIAMOND_RES

#
# Parsing diamond-blastp results.
#

# Parsed blastp results.
PARSED=${DIAMOND}/diamond-blastp-parsed.txt

# Parsing blastp results for cattle specific hits.
python -m biostar.tools.gene_pred.blast_parse --blast $DIAMOND_RES >$PARSED

# Main result view.
cat $PARSED | head -100 >result-preview.txt




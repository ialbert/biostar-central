# Stop on any error.
set -uxe

# Assembly scaffolds.
ASSEMBLY={{assembly.path}}

# Species to be used as training set.
SPECIES={{species.value}}

# No. of processors to be used.
NPROC={{processors.value}}

# Diamond NR database.
DIAMOND_NR={{runtime.local_root}}/diamond-dbs/nr/nr

# Protein-accession2taxon map file.
TAXON_MAP={{runtime.local_root}}/accession-maps/prot.accession2taxid.gz

# Size-sorted assembly fasta file.
ASSM_SORTED={{runtime.work_dir}}/assembly_sorted.fa

# Sorting assembly based on size.
cat $ASSEMBLY | bioawk -c fastx ' { print length($seq),$name,$seq } '  | sort -k1nr,1  | awk '{ print ">"$2"\n"$3"\n"}' > $ASSM_SORTED

# Augustus results directory.
AUGUSTUS={{runtime.work_dir}}/augustus

# Create AUGUSTUS directory.
mkdir -p $AUGUSTUS

# Augustus predicted genes file.
GENES=${AUGUSTUS}/genes.gff

# Run augustus gene prediction.
echo "Running augustus and predicting genes."
echo "-------------------------------------"
cat $ASSM_SORTED | parallel --j $NPROC --blocksize 5M --recstart '>' --pipe "cat {} > tmp/{%} && augustus --species=$SPECIES  tmp/{%}" > $GENES

# Bed file with predicted transcripts.
TRANS_BED=${AUGUSTUS}/transcripts.bed

# Fasta file with predicted transcripts.
TRANS_FASTA=${AUGUSTUS}/transcripts.fa

# Creating a bed file of augustus predicted transcripts.
cat $GENES | grep -v "#" | awk '$3=="transcript" {print $0}' |awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$1"_"$9,".\t"$7}' > $TRANS_BED

# Creating transcripts fasta file.
bedtools getfasta -fi $ASSM_SORTED -bed $TRANS_BED -name -s -fo tmp.fa

# Formatting transcripts fasta file.
seqtk seq -l 80 tmp.fa >$TRANS_FASTA
rm -f tmp.fa

# Diamond blastx results directory.
DIAMOND={{runtime.work_dir}}/blastx

# Diamond blastx results.
DIAMOND_RES=${DIAMOND}/diamond-blastx.txt

# Create DIAMOND directory.
mkdir -p $DIAMOND

# Writing blastx header into output file.
HEADER="qseqid sseqid stitle staxids pident qlen slen length qstart qend sstart send evalue bitscore score"
echo $HEADER  |tr [:blank:] \\t >$DIAMOND_RES

# Running diamond blastx on predicted transcripts.
diamond blastx -f 6 $HEADER -d $DIAMOND_NR --taxonmap $TAXON_MAP --max-target-seqs 15 -q $TRANS_FASTA -p 10 >>$DIAMOND_RES

# Parsing blastx results.




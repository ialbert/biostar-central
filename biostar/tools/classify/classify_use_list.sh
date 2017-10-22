set -ue
# This script takes an accession list  and builds a multifasta file of all genomes.
#  It also creates centrifuge index and runs centrifuge.

INPUT_DATA={{data.path}}
INPUT_SAMPLE_INFO={{sampleinfo.path}}
ACC_LIST={{taxamap.path}}

DATA_DIR=data
RESULT_DIR=results
UPDATED_SAMPLE_INFO=updated_sampleinfo.txt
REF_FILE=genomes.fa

# Get genome sequence.
for accession in $(cat $ACC_LIST)
do
	echo "getting sequence for $accession"
	efetch -db=nuccore -format=fasta -id=$accession >>REF_FILE

	# Get taxid of the genome.
	taxid=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=${accession}&amp;rettype=fasta&amp;retmode=xml" | grep TSeq_taxid | sed -e s'/<TSeq_taxid>//g' -e s'/\<\/TSeq_taxid\>//g' | tr -d '[:space:]' )

	# Create conversion table for centrifuge.
	echo -e "$accession\t$taxid" >>name2taxa.txt
done 

# Modify fasta header.
sed -i .bak 's/\..*//g' $OUT_FILE

# Download ncbi taxonomy.
centrifuge-download -o taxonomy taxonomy

# Build centrifuge index.
echo "Building centrifuge index."
centrifuge-build -p 4 --conversion-table name2taxa.txt --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp $OUT_FILE $REF_NAME

# Extract data.
tar -xzvf $INPUT_DATA

# Update sample sheet to include file names. The new sample sheet is named $(UPDATED_SAMPLE_INFO)
python -m biostar.tools.data.samplesheet $INPUT_SAMPLE_INFO $DATA_DIR > $UPDATED_SAMPLE_INFO

#  Run classification using centrifuge.
mkdir -p $RESULT_DIR
cat $UPDATED_SAMPLE_INFO |parallel --verbose --progress  --header : --colsep '\t' centrifuge \
-x $INDEX -q -1 {file1} -2 {file2} --report-file ${RESULT_DIR}/{sample_name}_report.txt \
-S ${RESULT_DIR}/{sample_name}_classify.txt
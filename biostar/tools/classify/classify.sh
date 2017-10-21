
INPUT_DATA={{{data.path}
INPUT_SAMPLE_INFO={{sampleinfo.path}}
REFERENCE={{reference.path}}
THREADS={{threads.value}}

DATA_DIR=data
RESULT_DIR=results
INDEX_DIR=index
INDEX_BASE=$INDEX_DIR/genomes
UPDATED_SAMPLE_INFO=updated_sampleinfo.txt

# Extract data.
tar -xzvf $INPUT_DATA

# Update sample sheet to include file names. The new sample sheet is named $(UPDATED_SAMPLE_INFO)
python -m biostar.tools.data.samplesheet $INPUT_SAMPLE_INFO $DATA_DIR > $UPDATED_SAMPLE_INFO

# Create index.
mkdir -p $INDEX_DIR
# Download taxonomy
centrifuge-download -o taxonomy taxonomy
# make accession2taxid table -acc2tax.txt
# TO DO
# centrifuge-build -p 6 --conversion-table acc2tax.txt --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp $REFERENCE $INDEX_BASE


#  Run classification using centrifuge.
mkdir -p $RESULT_DIR
cat $UPDATED_SAMPLE_INFO |parallel --verbose --progress -j $THREADS  --header : --colsep '\t' centrifuge \
-x $INDEX_BASE -q -1 ${DATA_DIR}/{file1} -2 ${DATA_DIR}/{file2} --report-file \
${RESULT_DIR}/{sample_name}_report.txt -S ${RESULT_DIR}/{sample_name}_classify.txt

# Register results to result/index.html
# TO DO



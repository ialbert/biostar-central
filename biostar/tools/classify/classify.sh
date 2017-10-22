
INPUT_DATA={{data.path}}
INPUT_SAMPLE_INFO={{sampleinfo.path}}
TAXA_MAP={{taxamap.path}}
REFERENCE={{reference.path}}
INDEX=${REFERENCE%.*}

DATA_DIR=data
RESULT_DIR=results
UPDATED_SAMPLE_INFO=updated_sampleinfo.txt

# Extract data.
tar -xzvf $INPUT_DATA

# Update sample sheet to include file names. The new sample sheet is named $(UPDATED_SAMPLE_INFO)
python -m biostar.tools.data.samplesheet $INPUT_SAMPLE_INFO $DATA_DIR > $UPDATED_SAMPLE_INFO

# Download ncbi taxonomy.
centrifuge-download -o taxonomy taxonomy

# Build centrifuge index.
centrifuge-build -p 4 --conversion-table $TAXA_MAP --taxonomy-tree taxonomy/nodes.dmp \
--name-table taxonomy/names.dmp $REFERENCE $INDEX

#  Run classification using centrifuge.
mkdir -p $RESULT_DIR
cat $UPDATED_SAMPLE_INFO |parallel --verbose --progress  --header : --colsep '\t' centrifuge \
-x $INDEX -q -1 {file1} -2 {file2} --report-file ${RESULT_DIR}/{sample_name}_report.txt \
-S ${RESULT_DIR}/{sample_name}_classify.txt

# Register results to result/index.html
# TO DO




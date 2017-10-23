# Get parameters.
INPUT_DATA={{data.path}}
TAXA_MAP={{taxamap.path}}
REFERENCE={{reference.path}}
TAX_DIR={{taxonomy.dirname}}
LOCAL_DIR={{runtime.local_root}}
RESULT_VIEW={{settings.index}}

# Internal parameters.
INDEX=${REFERENCE%.*}
TAX_NODES="${LOCAL_DIR}/${TAX_DIR}/nodes.dmp"
TAX_NAMES="${LOCAL_DIR}/${TAX_DIR}/names.dmp"

DATA_DIR=data
RESULT_DIR=results

# Build centrifuge index.
centrifuge-build -p 4 --conversion-table $TAXA_MAP --taxonomy-tree $TAX_NODES --name-table $TAX_NAMES $REFERENCE $INDEX

# Run classification using centrifuge.
mkdir -p $RESULT_DIR
centrifuge -x $INDEX -U $INPUT_DATA --report-file ${RESULT_DIR}/report.txt -S ${RESULT_DIR}/classification.txt

# Register results to settings.index.
cp ${RESULT_DIR}/report.txt ${RESULT_VIEW}



